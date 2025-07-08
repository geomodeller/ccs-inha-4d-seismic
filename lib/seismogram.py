import numpy as np
import pandas as pd
import os

class CCS_Seismogram():
    t_step, NZ, NY, NX = 6, 16, 32, 32
    R = 83.14462618 # gas constant (bar*cm^3/mol/K)
    DZ = np.ones((NZ, NY, NX)) * 1.905
    v_above = 6.45 # Unit: km/s
    dt = 0.001   #sampleing interval
    t_max = 3.0   # max time to create time vector
    t = np.arange(0, t_max, dt)
    f=30            # wavelet frequency
    length=0.512    # wavelet vector length
    
    v_above = 6.45 # Unit: km/s
    dt = 0.001   #sampleing interval
    t_max = 3.0   # max time to create time vector
    t = np.arange(0, t_max, dt)
    
    f=30            # wavelet frequency
    length=0.512    # wavelet vector length
    
    def __init__(self, 
                 P = None,
                 T = None,
                 S = None,
                 por = None,
                 Sg = None,
                 TOP_D = None,
                 current_directory = None,
                 heat_capacity_file = None
                 ):
        """
        Initializes the object with the current directory, template directory, and total injection amount.

        Parameters:
            currert_directory (str): The current directory path.
            template_dir (str): The template directory path.
            total_injection_ammount (int or float): The total injection amount.

        Returns:
            None
        """
        if P.dtype in [int, float]:
            self.P = P
            self.Ppr_co2 = self.P/73.8
        else:
            assert False, "only numpy array is allowed and dtype must be int and float"
        if T.dtype in [int, float]:
            self.T = T
            self.Tpr_co2 = self.T/304.25
        else:
            assert False, "only numpy array is allowed and dtype must be int and float"
        if S.dtype in [int, float]:
            self.S = S
        else:
            assert False, "only numpy array is allowed and dtype must be int and float"            
        if por.dtype in [int, float]:
            self.por = por
        else:
            assert False, "only numpy array is allowed and dtype must be int and float"            
        if Sg.dtype in [int, float]:
            self.Sg = Sg
        else:
            assert False, "only numpy array is allowed and dtype must be int and float"            
        if TOP_D.dtype in [int, float]:
            self.TOP_D = TOP_D
        else:
            assert False, "only numpy array is allowed and dtype must be int and float"
        if current_directory is not None:
            self.current_directory = current_directory
        if heat_capacity_file is not None:
            self.heat_capacity_file = os.path.join(self.current_directory,heat_capacity_file)

        assert os.path.isfile(self.heat_capacity_file), f'{heat_capacity_file} does not exist'
    
    def reset_pressure(self, P):
        if P.dtype in [int, float]:
            self.P = P
            self.Ppr_co2 = self.P/73.8
        else:
            assert False, "only numpy array is allowed and dtype must be int and float"
            
    def reset_temperature(self, T):
        if T.dtype in [int, float]:
            self.T = T
            self.Tpr_co2 = self.T/304.25
        else:
            assert False, "only numpy array is allowed and dtype must be int and float"
            
    def reset_salinity(self, S):
        if S.dtype in [int, float]:
            self.S = S
        else:
            assert False, "only numpy array is allowed and dtype must be int and float"

    def _density_of_brine(self):
        P_mpa = self.P * 0.1 # conversion from bars to MPa
        T_c = self.T - 273.15 # conversion from kelvin to celcius
        rho_w = 1+10**(-6)*(-80*T_c-3.3*T_c**2+0.00175*T_c**3+489*P_mpa-2*T_c*P_mpa+0.016*T_c**2*P_mpa-1.3*10**(-5)*T_c**(3)*P_mpa-0.333*P_mpa**(2)-0.002*T_c*P_mpa**2) # Unit: g/cm3
        rho_b = rho_w + self.S*(0.668+0.44*self.S+10**(-6)*(300*P_mpa-2400*P_mpa*self.S+T_c*(80+3*T_c-3300*self.S-13*P_mpa+47*P_mpa*self.S)))
        return rho_b
    
    def _velocity_of_brine(self):
        V_b = -np.ones((self.t_step, self.NZ, self.NY, self.NX))
        P_mpa = self.P * 0.1 # conversion from bars to MPa
        T_c = self.T - 273.15 # conversion from kelvin to celcius
        coef = np.array([[1402.85, 1.524, 0.003437, -1.197*10**(-5)],
                        [4.871, -0.0111, 1.739*10**(-4), -1.628*10**(-6)],
                        [-0.04783, 2.747*10**(-4), -2.135*10**(-6), 1.237*10**(-8)],
                        [1.487*10**(-4), -6.503*10**(-7), -1.455*10**(-8), 1.327*10**(-10)],
                        [-2.197*10**(-7), 7.987*10**(-10), 5.23*10**(-11), -4.614*10**(-13)]])
        for t in range(self.t_step):
            for k in range(self.NZ):
                for j in range(self.NY):
                    for i in range(self.NX):
                        T_matrix = np.array([[1, T_c[t,k,j,i], T_c[t,k,j,i]**2, T_c[t,k,j,i]**3, T_c[t,k,j,i]**4]]).T
                        P_matrix = np.array([[1, P_mpa[t,k,j,i], P_mpa[t,k,j,i]**2, P_mpa[t,k,j,i]**3]])
                        V_w = np.sum(coef*T_matrix*P_matrix)
                        V_b_ = V_w + self.S[t,k,j,i]*(1170-9.6*T_c[t,k,j,i]+0.055*T_c[t,k,j,i]**2-8.5*10**(-5)*T_c[t,k,j,i]**3+ \
                                    2.6*P_mpa[t,k,j,i]-0.0029*T_c[t,k,j,i]*P_mpa[t,k,j,i]-0.0476*P_mpa[t,k,j,i]**2)+ \
                                    self.S[t,k,j,i]*(780-10*P_mpa[t,k,j,i]+0.16*P_mpa[t,k,j,i]**2)-820*self.S[t,k,j,i]**2
                        V_b[t,k,j,i] = V_b_/1000

        return V_b
    
    def _compressibility_of_co2(self):
        E = 0.109*(3.85-self.Tpr_co2)**2*np.exp(-(0.45+8*(0.56-1/self.Tpr_co2)**2)*self.Ppr_co2**(1.2)/self.Tpr_co2)
        Z_co2 = (0.03+0.00527*(3.5-self.Tpr_co2)**3)*self.Ppr_co2+(0.642*self.Tpr_co2-0.007*self.Tpr_co2**4-0.52)+E
        return Z_co2
    
    def _density_of_co2(self):
        rho_co2 = 44.01*self.P/(self.R*self.Z_co2*self.T)
        self.rho_co2 = rho_co2
        return rho_co2
    
    def _heat_capacity_of_co2(self):
        Cp_co2 = -np.ones((self.t_step, self.NZ, self.NY, self.NX))
        P_mpa = self.P * 0.1 # conversion from bars to MPa
        Cp_co2_df = pd.read_excel(os.path.join(self.current_directory, self.heat_capacity_file), sheet_name='Heat capacity', index_col=0)
        for t in range(self.t_step):
            for k in range(self.NZ):
                for j in range(self.NY):
                    for i in range(self.NX):
                        P_index = np.where(Cp_co2_df.index<=P_mpa[t,k,j,i])[0][-1]
                        T_index = np.where(Cp_co2_df.columns<=self.T[t,k,j,i])[0][-1]
                        T_interval =Cp_co2_df.columns[T_index: T_index+2]
                        P_interval =Cp_co2_df.index[P_index: P_index+2]
                        Cp_co2_array = Cp_co2_df.iloc[P_index:P_index+2,T_index:T_index+2].to_numpy()
                        Cp_backward=np.interp(P_mpa[t,k,j,i], P_interval, Cp_co2_array[:,0])
                        Cp_forward=np.interp(P_mpa[t,k,j,i], P_interval, Cp_co2_array[:,1])
                        Cp_co2[t,k,j,i]=np.interp(self.T[t,k,j,i], T_interval, [Cp_backward, Cp_forward]) #Unit: kJ/(kg*K)
        return Cp_co2
    
    def _velocity_of_co2(self):
        coef=-(0.45+8*(0.56-1/self.Tpr_co2)**2)/self.Tpr_co2
        dE_dP=0.109*(3.85-self.Tpr_co2)**(2)*1.2/73.8*coef*np.exp(coef*self.Ppr_co2**1.2)*self.Ppr_co2**(0.2)
        dZ_dP=(0.03+0.00527*(3.5-self.Tpr_co2)**3)/73.8+dE_dP

        beta_co2=1/self.P-(1/self.Z_co2*dZ_dP)        #Unit: 1/bar
        V_co2 = np.sqrt(1/(beta_co2*self.rho_co2))/100   #Unit: km/s
        self.V_co2 = V_co2
        return V_co2
    
    def _bulk_moduli_of_co2(self):
        Cp_co2 = self._heat_capacity_of_co2()
        coef=-(0.45+8*(0.56-1/self.Tpr_co2)**2)/self.Tpr_co2
        dZ_dT=(-0.01581*self.Ppr_co2*(3.5-self.Tpr_co2)**2-0.028*self.Tpr_co2**3+0.642-\
            0.218*(3.85-self.Tpr_co2)*np.exp(-(0.45+8*(0.56-1/self.Tpr_co2)**2)*self.Ppr_co2**(1.2)/self.Tpr_co2)-\
            0.109*self.Ppr_co2**(1.2)*(3.85-self.Tpr_co2)**2*(-2.9588/self.Tpr_co2**2+17.92/self.Tpr_co2**3-24/self.Tpr_co2**4)*np.exp(-(0.45+8*(0.56-1/self.Tpr_co2)**2)*self.Ppr_co2**(1.2)/self.Tpr_co2))/304.25
        dE_dP=0.109*(3.85-self.Tpr_co2)**(2)*1.2/73.8*coef*np.exp(coef*self.Ppr_co2**1.2)*self.Ppr_co2**(0.2)
        dZ_dP=(0.03+0.00527*(3.5-self.Tpr_co2)**3)/73.8+dE_dP
        alpha_co2=dZ_dT/self.Z_co2+1/self.T
        beta_co2=1/self.P-(1/self.Z_co2*dZ_dP)        #Unit: 1/bar
        gamma = 1/(1-(self.T**2*self.R*self.Z_co2*alpha_co2**2)/(Cp_co2*beta_co2*self.P*44.01*10))
        K_co2 = (gamma*self.P)/(1-self.P*dZ_dP/self.Z_co2)/10_000       #Unit: GPa
        self.gamma = gamma
        self.dZ_dP = dZ_dP
        return K_co2
    
    def _differential_pressure(self):
        g_acc = 9.81 #m/s^2
        rho_sand = 2.6
        por_up = np.ones((self.NY, self.NX))*0.1
        P_over = -np.ones((self.t_step, self.NZ, self.NY, self.NX))
        rho_rock = (1-self.por)*rho_sand + self.por*(self.Sg*self.rho_co2+(1-self.Sg)*self.rho_b)
        rho_up = por_up*self.rho_b[0, 0] + (1-por_up)*rho_sand #g/cm^3
        P_up =  rho_up * self.TOP_D * g_acc / 100 #bar
        for k in range(self.NZ):
            for j in range(self.NY):
                for i in range(self.NX):
                    P_over[:, k, j, i] = P_up[j, i] + 0.5*(rho_rock[:, k, j, i]*g_acc*self.DZ[k, j, i]) + np.sum((rho_rock[:, 1:k, j, i]*g_acc*self.DZ[1:k, j, i])) / 100
        return P_over - self.P
    
    def velocity_of_saturated_rock(self, mu_min = 37, K_min = 44):
        self.Z_co2 = self._compressibility_of_co2()
        self.rho_b = self._density_of_brine()
        self.rho_co2 = self._density_of_co2()

        rho_sand = 2.6
        bar_to_GPa = 0.0001
        por_c = 0.4
        C = 2.8/por_c
        nu = (3*K_min-2*mu_min)/(6*K_min+2*mu_min) # Poisson's ratio of mineral
        
        V_b = self._velocity_of_brine()
        K_b = self.rho_b*V_b**2      #Unit: GPa

        Pd = self._differential_pressure()
        K_co2 = self._bulk_moduli_of_co2()
        K_mc = ((C**2*(1-por_c)**2*mu_min**2*Pd*bar_to_GPa)/(18*np.pi**2*(1-nu)**2))**(1/3) #GPa
        mu_mc = (5-4*nu)/(5*(1-nu))*((3*C**2*(1-por_c)**2*mu_min**2*Pd*bar_to_GPa)/(2*np.pi**2*(1-nu)**2))**(1/3) #GPa
        
        K_f = (self.por/por_c/(K_mc+4/3*mu_mc)+(1-self.por/por_c)/(K_min+4/3*mu_mc))**(-1)-4/3*mu_mc
        mu = (self.por/por_c/(mu_mc+mu_mc/6*((9*K_mc+8*mu_mc)/(K_mc+2*mu_mc)))+((1-self.por/por_c)/(mu_min+mu_mc/6*((9*K_mc+8*mu_mc)/(K_mc+2*mu_mc)))))**(-1)-mu_mc/6*((9*K_mc+8*mu_mc)/(K_mc+2*mu_mc))
        rho_fl = (1-self.Sg)*self.rho_b+self.Sg*self.rho_co2
        K_fl = 1/((1-self.Sg)/K_b+self.Sg/K_co2)
        
        K_sat = K_f + (1-K_f/K_min)**2/(self.por/K_fl+(1-self.por)/K_min-K_f/K_min**2)
        self.rho_sat = (1-self.por)*rho_sand + self.por*(self.Sg*self.rho_co2+(1-self.Sg)*self.rho_b)
        # self.rho_sat = rho_rock+self.por*(rho_fl-self.rho_b) #Unit: g/cm^3/
        self.K_b = K_b
        self.K_co2 = K_co2
        self.K_sat =K_sat
        self.mu = mu
        
        Vp = np.sqrt((K_sat+4*mu/3)/self.rho_sat)
        return Vp
        
    def ricker(self):
        t0 = np.arange(-self.length/2, (self.length-self.dt)/2, self.dt)
        y = (1.0 - 2.0*(np.pi**2)*(self.f**2)*(t0**2)) * np.exp(-(np.pi**2)*(self.f**2)*(t0**2))
        return t0, y

    def two_way_travel_time(self, DZ, TOP_D):
        TWT = -np.ones((self.t_step, self.NZ, self.NY, self.NX))
        for t in range(self.t_step):
            for i in range(self.NX):
                for j in range(self.NY):
                    for k in range(self.NZ):
                        twt = 0
                        for s in range(k+1):
                            twt += 2*DZ[s, j, i]/1000/self.Vp[t, s, j, i]
                        TWT[t, k, j, i] = 2*(TOP_D[j,i]/1000/self.v_above) + twt
        assert TWT.min() > 0, 'negative TWT is detected'
        return TWT
    
    def acoustic_impedance(self, DZ, TOP_D):
        TWT = self.two_way_travel_time(DZ, TOP_D)
        AI = self.Vp * self.rho_sat * 1000 # Pa*s/m3
        AI_tdom = -np.ones((self.t_step, self.t.shape[0], self.NY, self.NX))
        for t in range(self.t_step):
            for i in range(self.NX):
                for j in range(self.NY):
                    AI_tdom[t, :, j, i] = np.interp(x=self.t, xp=TWT[t, :,j,i], fp = AI[t, :,j,i])
        assert AI_tdom.min() > 0, 'negative AI is detected'
        
        return AI_tdom
        
    def reflectivity_coef(self, DZ, TOP_D):
        AI_tdom = self.acoustic_impedance(DZ, TOP_D)
        Rc_tdom = -np.ones_like(AI_tdom)
        for t in range(self.t_step):
            Rc_tdom[t, :-1] = (AI_tdom[t, 1:]-AI_tdom[t, :-1])/(AI_tdom[t, 1:]+AI_tdom[t, :-1]) 
            Rc_tdom[t,-1] = Rc_tdom[t, -2]
        return Rc_tdom
    
    def seismic(self, DZ, TOP_D):
        self.Vp = self.velocity_of_saturated_rock()
        TWT = self.two_way_travel_time(DZ, TOP_D)
        Rc_tdom = self.reflectivity_coef(DZ, TOP_D)
        seismic_tdom = -100*np.ones_like(Rc_tdom)
        t0, w = self.ricker()
        for t in range(self.t_step):
            for i in range(self.NX):
                for j in range(self.NY):
                    seismic_tdom[t, :, j, i] = np.convolve(w, Rc_tdom[t, :, j, i], mode='same')
        assert seismic_tdom.min() > -5 and seismic_tdom.max() < 5 , 'erroneous seismic is detected, seismic should be oscillated around zero'
        
        seismic = -np.ones((self.t_step, self.NZ, self.NY, self.NX))
        for t in range(self.t_step):
            for i in range(self.NX):
                for j in range(self.NY):
                    seismic[t, :, j, i] = np.interp(x=TWT[t, :, j, i], xp = self.t, fp = seismic_tdom[t, :, j, i])
        return seismic_tdom, seismic

def calc_noise_std(power_signal, snr):
    return np.sqrt(power_signal/snr)

def calc_noise_seismic(noise_std, nx=32, ny=32, nz=96):
    '''
    std1 = seismic.std()
    std2 = 0.001
    power_signal = np.mean(seismic**2)
    tol = 1000
    while tol > 0.0001:
        new_std = (std1+std2)/2
        new_power_noise = np.mean(np.random.normal(0, new_std, seismic.shape)**2)
        new_SNR_diff = power_signal/new_power_noise - SNR
        if new_SNR_diff <0:
            std1 = new_std
        else:
            std2 = new_std
        tol = abs(new_SNR_diff)
    print(f'SNR: {power_signal/new_power_noise}')
    return seismic + np.random.normal(0, new_std, seismic.shape)
    '''
    return np.random.normal(0, noise_std, size=(nz, ny, nx))

def add_noise_to_seismic(seismic, power_signal, snr):
    noise_std = calc_noise_std(power_signal, snr)
    (nz, ny, nx) = seismic.shape
    return seismic + calc_noise_seismic(noise_std, nx, ny, nz)