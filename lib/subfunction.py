import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import scipy.stats as stats
import torch
def density_of_brine(P = 100., S = 0.05, T = 70, ):
    """ to compute density of brine
    
    Args:
        P (float): Presssure in bar. Defaults to 100.
        S (float): salinity in mass frac. Defaults to 0.05.
        T (float): temperature in Kelvin. Defaults to 273.15.

    Returns:
        float: density of brine
    """
    P = P * 0.1 # conversion from bars to MPa
    T = T - 273.15 # conversion from kelvin to celcius
    rho_w = 1+10**(-6)*(-80*T-3.3*T**2+0.00175*T**3+489*P-2*T*P+0.016*T**2*P-1.3*10**(-5)*T**(3)*P-0.333*P**(2)-0.002*T*P**2) # Unit: g/cm3
    rho_b = rho_w + S*(0.668+0.44*S+10**(-6)*(300*P-2400*P*S+T*(80+3*T-3300*S-13*P+47*P*S)))
    return rho_b

def velocity_of_brine(P, S, T):
    """ to compute sonic velocity of brine based on (Wilson et al., 1959)
    This method can be used for Temperatures up to 100C and pressures up to 100MPa

    Args:
        P (float): Presssure in bar. Defaults to 100.
        S (float): salinity in mass frac. Defaults to 0.05.
        T (float): temperature in Kelvin. Defaults to 273.15.

    Returns:
        float: sonic velocity of brine in km/sec
    """
    T = T - 273.15 # conversion from kelvin to celcius
    P = P * 0.1 # conversion from bars to MPa
    coef = np.array([[1402.85, 1.524, 0.003437, -1.197*10**(-5)],
                    [4.871, -0.0111, 1.739*10**(-4), -1.628*10**(-6)],
                    [-0.04783, 2.747*10**(-4), -2.135*10**(-6), 1.237*10**(-8)],
                    [1.487*10**(-4), -6.503*10**(-7), -1.455*10**(-8), 1.327*10**(-10)],
                    [-2.197*10**(-7), 7.987*10**(-10), 5.23*10**(-11), -4.614*10**(-13)]])
    T_matrix = np.array([[1, T, T**2, T**3, T**4]]).T
    P_matrix = np.array([[1, P, P**2, P**3]])
    V_w = np.sum(coef*T_matrix*P_matrix)
    V_b = V_w + S*(1170-9.6*T+0.055*T**2-8.5*10**(-5)*T**3+ \
                2.6*P-0.0029*T*P-0.0476*P**2)+ \
                S**1.5*(780-10*P+0.16*P**2)-820*S**2
    return V_b/1000

def bulk_moduli_of_brine():
    return

def compressibility_of_co2(P, T):
    """calculate Z-factor of co2 based on (Thomas et al., 1970)

    Args:
        P (float): Presssure in bar. Defaults to 100.
        T (float): temperature in kelvin. Defaults to 273.15.
        
    Returns:
        float: compressibility(z) of co2
    """
    G = 1.8397145625/1.2041 # the ratio between densities of co2 and air at the surface condition
    Ppr_co2 = P/(48.92-4.048*G)
    Tpr_co2 = T/(94.72+170.75*G)
    E = 0.109*(3.85-Tpr_co2)**2*np.exp(-(0.45+8*(0.56-1/Tpr_co2)**2)*Ppr_co2**(1.2)/Tpr_co2)
    Z_co2 = (0.03+0.00527*(3.5-Tpr_co2)**3)*Ppr_co2+(0.642*Tpr_co2-0.007*Tpr_co2**4-0.52)+E
    return Z_co2

def density_of_co2(P,T):
    """calculate density of co2 based on (Thomas et al., 1970)

    Args:
        P (float): Presssure in bar. Defaults to 100.
        T (float): temperature in kelvin. Defaults to 273.15.
        
    Returns:
        float: density of co2 in g/cm3
    """
    R = 83.14462618 # gas constant (bar*cm^3/mol/K)
    Z_co2 = compressibility_of_co2(P, T)
    rho_co2 = 44.01*P/(R*Z_co2*T)
    
    return rho_co2

def heat_capacity_of_co2(directory:'Heat_capacity_of_co2.xlsx', P, T):
    """_summary_

    Args:
        P (float): presssure in bar. Defaults to 100.
        T (float): temperature in kelvin. Defaults to 273.15.
        directory (str): reference table directory of the heat capacity of co2 

    Returns:
        float: heat capacity of co2 in kJ/(kg*K)
    """
    Cp_co2 = -np.ones(P.shape)
    P = P * 0.1 # conversion from bars to MPa
    Cp_co2_df = pd.read_excel(directory, sheet_name='Heat capacity', index_col=0)
    for t in range(P.shape[0]):
        for k in range(P.shape[1]):
            for j in range(P.shape[2]):
                for i in range(P.shape[3]):
                    P_index = np.where(Cp_co2_df.index<=P[t,k,j,i])[0][-1]
                    T_index = np.where(Cp_co2_df.columns<=T)[0][-1]

                    T_interval =Cp_co2_df.columns[T_index: T_index+2]
                    P_interval =Cp_co2_df.index[P_index: P_index+2]
                    Cp_co2_array = Cp_co2_df.iloc[P_index:P_index+2,T_index:T_index+2].to_numpy()

                    Cp_backward=np.interp(P[t,k,j,i], P_interval, Cp_co2_array[:,0])
                    Cp_forward=np.interp(P[t,k,j,i], P_interval, Cp_co2_array[:,1])
                    Cp_co2[t,k,j,i]=np.interp(T, T_interval, [Cp_backward, Cp_forward]) #Unit: kJ/(kg*K)
                        
                    return Cp_co2

def velocity_of_co2(P, T):
    G = 1.8397145625/1.2041 # the ratio between densities of co2 and air at the surface condition
    Ppr_co2 = P/(48.92-4.048*G)
    Tpr_co2 = T/(94.72+170.75*G)
    
    coef=-(0.45+8*(0.56-1/Tpr_co2)**2)/Tpr_co2
    dE_dP=0.109*(3.85-Tpr_co2)**(2)*1.2/(48.92-4.048*G)*coef*np.exp(coef*Ppr_co2**1.2)*Ppr_co2**(0.2)
    dZ_dP=(0.03+0.00527*(3.5-Tpr_co2)**3)/(48.92-4.048*G)+dE_dP

    Z_co2 = compressibility_of_co2(P, T)
    rho_co2 = density_of_co2(P, T)
    
    beta_co2=1/P-(1/Z_co2*dZ_dP)        #Unit: 1/bar
    V_co2 = np.sqrt(1/(beta_co2*rho_co2))/100   #Unit: km/s
    return V_co2

def bulk_moduli_of_co2(P, T):
    R=0.08311446 # R: gas constant Unit: L*bar*K**(-1)*mol**(-1)
    G = 1.8397145625/1.2041 # the ratio between densities of co2 and air at the surface condition
    Ppr_co2 = P/(48.92-4.048*G)
    Tpr_co2 = T/(94.72+170.75*G)
    
    coef=-(0.45+8*(0.56-1/Tpr_co2)**2)/Tpr_co2
    dZ_dT=(-0.01581*Ppr_co2*(3.5-Tpr_co2)**2-0.028*Tpr_co2**3+0.642-\
        0.218*(3.85-Tpr_co2)*np.exp(-(0.45+8*(0.56-1/Tpr_co2)**2)*Ppr_co2**(1.2)/Tpr_co2)-\
        0.109*Ppr_co2**(1.2)*(3.85-Tpr_co2)**2*(-2.9588/Tpr_co2**2+17.92/Tpr_co2**3-24/Tpr_co2**4)*np.exp(-(0.45+8*(0.56-1/Tpr_co2)**2)*Ppr_co2**(1.2)/Tpr_co2))/(94.72+170.75*G)
    dE_dP=0.109*(3.85-Tpr_co2)**(2)*1.2/(48.92-4.048*G)*coef*np.exp(coef*Ppr_co2**1.2)*Ppr_co2**(0.2)
    dZ_dP=(0.03+0.00527*(3.5-Tpr_co2)**3)/(48.92-4.048*G)+dE_dP
    
    Z_co2 = compressibility_of_co2(P, T)
    Cp_co2=heat_capacity_of_co2('Heat_capacity_of_co2.xlsx', P, T)
    alpha_co2=dZ_dT/Z_co2+1/T
    beta_co2=1/P-(1/Z_co2*dZ_dP)        #Unit: 1/bar
    
    gamma = 1/(1-(T**2*R*Z_co2*alpha_co2**2*100)/(Cp_co2*beta_co2*P*44.01))
    
    K_co2 = (gamma*P)/(1-P*dZ_dP/Z_co2)/10000       #Unit: GPa
    
    return K_co2

def upscale_3d(target_array, upscale_factor = [6, 1, 1], method = 'mean'):
    assert method in ['mean'], "method should be 'mean', 'median', ..."
    assert target_array.shape[2]%upscale_factor[0] == 0, f'1st-dimension is not divided by {upscale_factor[0]}'
    assert target_array.shape[3]%upscale_factor[1] == 0, f'2nd-dimension is not divided by {upscale_factor[1]}'
    assert target_array.shape[4]%upscale_factor[2] == 0, f'3rd-dimension is not divided by {upscale_factor[2]}'
    
    n_sample, n_prop, nz, ny, nx = target_array.shape
    z_iter = nz//upscale_factor[0]
    y_iter = ny//upscale_factor[1]
    x_iter = nx//upscale_factor[2]
    upscaled_array = np.ones((n_sample, n_prop, z_iter, y_iter, x_iter)) * np.nan
    
    if method == 'mean':
        for sample in range(n_sample):
            for prop in range(n_prop):
                for k in range(z_iter):
                    for j in range(y_iter):
                        for i in range(x_iter):
                                upscaled_array[sample, prop, k, j, i] = np.mean(target_array[sample,
                                                                                             prop,
                                                                                             upscale_factor[0]*k:upscale_factor[0]*(k+1),
                                                                                             upscale_factor[1]*j:upscale_factor[1]*(j+1),
                                                                                             upscale_factor[2]*i:upscale_factor[2]*(i+1)])
                
    return upscaled_array

def _make_cdf(prop):
    prop_sorted = np.sort(prop.flatten())
    Cum_probability = np.linspace(0, 1, len(prop_sorted))
    
    return prop_sorted, Cum_probability

def cdf_mapping(pred,y):
    Cum_prob_pred = _make_cdf(pred)
    Cum_prob_y = _make_cdf(y)
    cdf_mapping = interp1d(pred, Cum_prob_pred.flatten(), bounds_error=False)
    inverse_cdf_mapping = interp1d(Cum_prob_y, y.flatten(), bounds_error=False)
    pred_mapped = inverse_cdf_mapping(cdf_mapping(pred))
    
    return pred_mapped

def cdf_mapping_by_facies(por, por_truth, facies, facies_truth):
    if  (por.shape == facies.shape):
        sand_por = por[facies != 0]
        shale_por= por[facies == 0]

        sand_por_truth= por_truth*facies_truth
        sand_por_truth=sand_por_truth[sand_por_truth!=0]
        
        shale_por_truth= por_truth*(1-facies_truth)
        shale_por_truth=shale_por_truth[shale_por_truth!=0]
        
        sand_por_sorted, Cum_prob_por_sand = _make_cdf(sand_por)
        shale_por_sorted, Cum_prob_por_shale = _make_cdf(shale_por)
        sand_por_truth_sorted, Cum_prob_por_sand_truth = _make_cdf(sand_por_truth)
        shale_por_truth_sorted, Cum_prob_por_shale_truth = _make_cdf(shale_por_truth)
        
        model = np.ones_like(por)*-999
        
        # for sand
        cdf_mapping = interp1d(sand_por_sorted, Cum_prob_por_sand, bounds_error=False)
        inverse_cdf_mapping = interp1d(Cum_prob_por_sand_truth, sand_por_truth_sorted, bounds_error=False)
        sand_por_mapped = inverse_cdf_mapping(cdf_mapping(sand_por))
        model[facies == 1] = sand_por_mapped
    
        # for shale
        cdf_mapping = interp1d(shale_por_sorted, Cum_prob_por_shale, bounds_error=False)
        inverse_cdf_mapping = interp1d(Cum_prob_por_shale_truth, shale_por_truth_sorted, bounds_error=False)
        shale_por_mapped = inverse_cdf_mapping(cdf_mapping(shale_por))
        model[facies != 1] = shale_por_mapped
        
        return model
    else:
        assert False, 'Something went wrong'
        

def normal_score_transform_by_facies(por, facies, por_truth, facies_truth, inverse = False):
    sand_por = por*facies
    sand_por = sand_por[sand_por!=0]
    
    shale_por= por*(1-facies)
    shale_por=shale_por[shale_por!=0]
    
    sand_por_truth = por_truth*facies_truth
    sand_por_truth=sand_por_truth[sand_por_truth!=0]
    
    shale_por_truth= por_truth*(1-facies_truth)
    shale_por_truth=shale_por_truth[shale_por_truth!=0]
    
    sand_por_truth_sorted, Cum_prob_por_sand_truth = _make_cdf(sand_por_truth)
    shale_por_truth_sorted, Cum_prob_por_shale_truth = _make_cdf(shale_por_truth)
    
    model = np.ones_like(por)*-999
    rv = stats.norm(0, 1)
    if inverse == True:
        sand_por_prob = rv.cdf(sand_por)
        shale_por_prob = rv.cdf(shale_por)
        inverse_cdf_map_sand = interp1d(Cum_prob_por_sand_truth, sand_por_truth_sorted, bounds_error=False)
        inverse_cdf_map_shale = interp1d(Cum_prob_por_shale_truth, shale_por_truth_sorted, bounds_error=False)
        sand_por_recovery = inverse_cdf_map_sand(sand_por_prob)
        shale_por_recovery = inverse_cdf_map_shale(shale_por_prob)
        model[facies == 1] = sand_por_recovery
        model[facies != 1] = shale_por_recovery
        
        return model
    else:
        cdf_mapping_sand = interp1d(sand_por_truth_sorted, Cum_prob_por_sand_truth, bounds_error=False)
        sand_por_prob = cdf_mapping_sand(sand_por)
        sand_por_prob[sand_por_prob == 1] = 0.999
        sand_por_prob[sand_por_prob == 0] = 0.001
        sand_por_nst = rv.ppf(sand_por_prob)
        cdf_mapping_shale = interp1d(shale_por_truth_sorted, Cum_prob_por_shale_truth, bounds_error=False)
        shale_por_prob = cdf_mapping_shale(shale_por)
        shale_por_prob[shale_por_prob == 1] = 0.999
        shale_por_prob[shale_por_prob == 0] = 0.001
        shale_por_nst = rv.ppf(shale_por_prob)
        model[facies == 1] = sand_por_nst
        model[facies != 1] = shale_por_nst
        
        return model
    
def normal_score_transform(seismic, seismic_truth, inverse = False):
    seismic_truth_sorted, Cum_prob_seismic_truth = _make_cdf(seismic_truth)
    rv = stats.norm(0, 1)
    if inverse == True:
        pass
        # sand_seismic_prob = rv.cdf(sand_seismic*3)
        # shale_seismic_prob = rv.cdf(shale_seismic*3)
        # inverse_cdf_map_sand = interp1d(Cum_prob_seismic_sand_truth, sand_seismic_truth_sorted)
        # inverse_cdf_map_shale = interp1d(Cum_prob_seismic_shale_truth, shale_seismic_truth_sorted)
        # sand_seismic_recovery = inverse_cdf_map_sand(sand_seismic_prob)
        # shale_seismic_recovery = inverse_cdf_map_shale(shale_seismic_prob)
        # model[facies == 1] = sand_seismic_recovery
        # model[facies != 1] = shale_seismic_recovery
        
        # return model
    else:
        cdf_mapping = interp1d(seismic_truth_sorted, Cum_prob_seismic_truth, bounds_error=False)
        seismic_prob = cdf_mapping(seismic)
        seismic_prob[seismic_prob == 1] = 0.999
        seismic_prob[seismic_prob == 0] = 0.001
        seismic_nst = rv.ppf(seismic_prob)

        return seismic_nst/3
    
def por_to_perm(por, facies):
    sand_por = por[facies == 1]
    shale_por = por[facies != 1]
    
    sand_perm = np.exp(3.379*np.log(sand_por)+8.53005018)
    shale_perm = np.exp(3.552*np.log(shale_por)+6.7300885)
    
    model = np.ones_like(por)*-9999
    model[facies == 1] = sand_perm
    model[facies != 1] = shale_perm
    
<<<<<<< HEAD
    return model
=======
    return model

def add_noise_to_seismic(seismic:np.ndarray, SNR:float) -> float:
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
>>>>>>> bc9ebff7d97c46052b17dd92ce71b11c0dd15fa6
