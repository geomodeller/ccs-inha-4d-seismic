######################################################################
## Description:                                                     ##
##   Fuctions are for calculating viscosity and density of CO2      ##
##   with the given pressure and temperature                        ##
## Authors: Jinwoo Eh, Eunsil Park, Honggeun at Inha University     ##
## Reference:  Heidaryan et al. (2011) and  Bahadori et al.(2009)   ##
## Last Update Date: 09/25/2023                                     ##
######################################################################

import numpy as np

def triple_point():
    P_t = 0.51795  # in MPa
    T_t = 216.592  # in Kelvin
    rho_t = [1178.53, 13.7614]  # in kg/m3 [Liquid and Gas] 
    return P_t, T_t, rho_t

def critical_point():
    P_c = 7.3773    # in MPa
    T_c = 304.1282  # in Kelvin
    rho_c = 467.6   # in kg/m3
    return P_c, T_c, rho_c

def melting_pressure(T):
    P_t, T_t, rho_t = triple_point()
    assert all(T > T_t), "the current T is lower than T_t, i.e., no melting curves exists"
    a_1 = 1955.5390
    a_2 = 2055.4593
    P_m = (1 + a_1*(T/T_t-1) + a_2*(T/T_t-1)**2)*P_t
    return P_m


def sublimation_pressure(T):
    P_t, T_t, rho_t = triple_point()
    assert all(T < T_t), "the current T is higher than T_t, i.e., no sublimation curves exists"
    a_1 = -14.740846
    a_2 = 2.4327015
    a_3 = -5.3061778
    P_sub = np.exp((a_1*(1-T/T_t) + a_2*(1-T/T_t)**(1.9) + a_3*(1-T/T_t)**(2.9))*T_t/T) * P_t
    return P_sub

def vapor_pressure(T):
    P_t, T_t, rho_t = triple_point()
    P_c, T_c, rho_c = critical_point()
    assert all(T > T_t), "the current T is lower than T_t, i.e., no vapor curves exists"
    a_1 = -7.0602087
    a_2 = 1.9391218
    a_3 = -1.6463597
    a_4 = -3.2995634
    t_1 = 1.0
    t_2 = 1.5
    t_3 = 2.0
    t_4 = 4.0
    P_s = np.exp((a_1*(1-T/T_c)**t_1 + a_2*(1-T/T_c)**t_2 + a_3*(1-T/T_c)**t_3 + a_4*(1-T/T_c)**t_4)*T_c/T) * P_c
    return P_s


def viscosity_of_co2(P, T, verbose = False):
    """this is to calculate viscosity of co2

    Args:
        p (float): Pressure [bar]
        T (float): Temperature [K]
        verbose (bool, optional): Whether to print out result. Defaults to False.

    Returns:
        viscosity (float): calculated viscosity of CO2 [cp]
    """
    A1 = -1.146067e-1
    A2 = 6.978380e-7
    A3 = 3.976765e-10
    A4 = 6.336120e-2
    A5 = -1.166119e-2
    A6 = 7.142596e-4
    A7 = 6.519333e-6
    A8 = -3.567559e-1
    A9 = 3.180473e-2
    
    viscosity = (A1 + A2*P + A3*P**2 + A4*np.log(T) + A5*np.log(T)**2 + A6*np.log(T)**3) / (1 + A7*P + A8*np.log(T) + A9* np.log(T)**2)

    if verbose == True: 
        print(f'viscosity: {viscosity}')

    return viscosity


def density_of_co2(P, T, verbose = False):
    """ this is to calculate density of co2

    Args:
        P (float): Pressure [bar]
        T (float): Temperature [K]
        verbose (bool, optional): Whether to print out result. Defaults to False.

    Returns:
        density (float): calculated density of CO2 [kg/m^3]
    """
    if P>25 and P<=100:
        A1 = 2.089800972761597e5;   B1 = -1.456286332143609e4; C1 = 2.885813588280259e2;   D1 = -1.597103845187521
        A2 = -1.675182353338921e3;  B2 = 1.16799554255704e2;   C2 = -2.31558333122805;     D2 = 1.284012022012305e-2
        A3 = 4.450600950630782;     B3 = -3.10430147581379e-1; C3 = 6.157718845508209e-3;  D3 = -3.420339567335051e-5
        A4 = -3.919844561756813e-3; B4 = 2.734973744483903e-4; C4 = -5.428007373890436e-6; D4 = 3.019572090945029e-8
        
        alpha = A1 + B1*P + C1*P**2 + D1*P**3
        beta  = A2 + B2*P + C2*P**2 + D2*P**3
        gamma = A3 + B3*P + C3*P**2 + D3*P**3
        theta = A4 + B4*P + C4*P**2 + D4*P**3

        density = alpha + beta*T + gamma*(T**2) + theta*(T**3)

    
    elif P>100 and P<700:
        A1 = 1.053293651041897e5;   B1 = -9.396448507019846e2;  C1 = 2.397414334181339;     D1 = -1.819046028481314e-3
        A2 = -8.253383504614545e2;  B2 = 7.618125848567747;     C2 = -1.963563757655062e-2; D2 = 1.497658394413360e-5
        A3 = 2.135712083402950;     B3 = -2.023128850373911e-2; C3 = 5.272125417813041e-5;  D3 = -4.043564072108339e-8
        A4 = -1.827956524285481e-3; B4 = 1.768297712855951e-5;  C4 = -4.653377143658811e-8; D4 = 3.586708189749551e-11

        alpha = A1 + B1*P + C1*P**2 + D1*P**3
        beta  = A2 + B2*P + C2*P**2 + D2*P**3
        gamma = A3 + B3*P + C3*P**2 + D3*P**3
        theta = A4 + B4*P + C4*P**2 + D4*P**3

        density = alpha + beta*T + gamma*(T**2) + theta*(T**3)

    else:
        assert False, "Warring: co2 pressure is not valid"
    
    if verbose == True:
        print(f'density: {density:.2f}')
    else: 
        return density

def density_of_co2_mesh(P, T, verbose = False):
    """ this is to calculate density of co2

    Args:
        P (float): Pressure [bar]
        T (float): Temperature [K]
        verbose (bool, optional): Whether to print out result. Defaults to False.

    Returns:
        density (float): calculated density of CO2 [kg/m^3]
    """
    A1 = 2.089800972761597e5;   B1 = -1.456286332143609e4; C1 = 2.885813588280259e2;   D1 = -1.597103845187521
    A2 = -1.675182353338921e3;  B2 = 1.16799554255704e2;   C2 = -2.31558333122805;     D2 = 1.284012022012305e-2
    A3 = 4.450600950630782;     B3 = -3.10430147581379e-1; C3 = 6.157718845508209e-3;  D3 = -3.420339567335051e-5
    A4 = -3.919844561756813e-3; B4 = 2.734973744483903e-4; C4 = -5.428007373890436e-6; D4 = 3.019572090945029e-8
    
    alpha = A1 + B1*P + C1*P**2 + D1*P**3
    beta  = A2 + B2*P + C2*P**2 + D2*P**3
    gamma = A3 + B3*P + C3*P**2 + D3*P**3
    theta = A4 + B4*P + C4*P**2 + D4*P**3

    density_1 = alpha + beta*T + gamma*(T**2) + theta*(T**3)

    
    A1 = 1.053293651041897e5;   B1 = -9.396448507019846e2;  C1 = 2.397414334181339;     D1 = -1.819046028481314e-3
    A2 = -8.253383504614545e2;  B2 = 7.618125848567747;     C2 = -1.963563757655062e-2; D2 = 1.497658394413360e-5
    A3 = 2.135712083402950;     B3 = -2.023128850373911e-2; C3 = 5.272125417813041e-5;  D3 = -4.043564072108339e-8
    A4 = -1.827956524285481e-3; B4 = 1.768297712855951e-5;  C4 = -4.653377143658811e-8; D4 = 3.586708189749551e-11

    alpha = A1 + B1*P + C1*P**2 + D1*P**3
    beta  = A2 + B2*P + C2*P**2 + D2*P**3
    gamma = A3 + B3*P + C3*P**2 + D3*P**3
    theta = A4 + B4*P + C4*P**2 + D4*P**3

    density_2 = alpha + beta*T + gamma*(T**2) + theta*(T**3)
    
    density = np.ones_like(P) * (-1)
    density[(P>25) & (P<=100)] = density_1[(P>25) & (P<=100)]
    density[(P>100) & (P<700)] = density_2[(P>100) & (P<700)] 

    if verbose == True: 
        print(f'viscosity: {density}')

    return density