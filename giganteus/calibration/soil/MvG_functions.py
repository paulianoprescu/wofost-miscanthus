import numpy as np
print("MvG_functions.py loaded, numpy version:", np.__version__)

def get_wc_MvG(H, WCR, WCS, ALPHA, NPAR):
    """
    Calculate the water content (wc) based on pressure head (H) using the Mualem van Genuchten model.

    Parameters:
    H: Pressure head [cm].
    WCR: Residual water content [cm3 cm-3].
    WCS: Saturated water content [cm3 cm-3].
    ALPHA: Curve shape parameter [-].
    NPAR: Curve shape parameter [-].

    Returns:
    wc: Water content [cm3 cm-3].
    """
    m = 1 - 1 / NPAR
    wc = WCR + (WCS - WCR) / ((1 + (ALPHA * H)**NPAR)**m)
    return wc

def get_cond_MvG(H, ALPHA, NPAR, LAMBDA, KSAT):
    """
    Calculate the hydraulic conductivity (cond) as a function of pressure head (H) using the Mualem van Genuchten model.
    
    Parameters:
    H: Pressure head [cm].
    ALPHA: Curve shape parameter [-].
    NPAR: Curve shape parameter [-].
    LAMBDA: Exponent in hydraulic conductivity function [-].
    KSAT: Hydraulic conductivity of saturated soil[cm day-1].

    Returns:
    cond: Hydraulic conductivity [cm day-1].
    """
    m = 1 - 1 / NPAR
    ah = ALPHA * H 
    h1 = (1 + ah**NPAR)**m
    h2 = ah**(NPAR - 1)
    denom = (1 + ah**NPAR)**(m*(LAMBDA + 2))
    cond = KSAT * (h1 - h2)**2 / denom
    return cond   

def generate_pcse_tables(WCR, WCS, ALPHA, NPAR, LAMBDA, KSAT):
    """
    Generate SMTAB and CONTAB lookup tables for PCSE.

    PCSE expects:
        SMTAB: volumetric soil moisture content as a function of pF [log (cm); cm3 cm-3].
        CONTAB: 10-log hydraulic conductivity as a function of pF [log (cm); log (cm day-1)].
    
    Returns:
    smtab
    contab
    smw: soil moisture content at wilting point [cm3/cm3].
    smfcf: soil moisture content at field capacity [cm3/cm3].
    sm0: soil moisture content at saturation [cm3/cm3].
    k0: hydraulic conductivity at saturation [cm day-1]. 
    """
    
    # Set sequence of pressure head values
    pF_values = np.array([-1.0, 1.0, 1.3, 1.491, 1.7, 2.0, 2.4, 2.7, 3.0, 3.4, 3.7, 4.0, 4.204, 6.0]) #standard used in WOFOST files
    H_values = 10 ** pF_values # convert pF to pressure head in cm

    smtab = []
    contab = []

    for pF, H in zip(pF_values, H_values):
        theta = get_wc_MvG(H, WCR, WCS, ALPHA, NPAR)
        k = get_cond_MvG(H, ALPHA, NPAR, LAMBDA, KSAT)
        smtab.append((pF, round(theta, 3)))
        contab.append((pF, round(k, 3)))

    smw = round(get_wc_MvG(10**4.2, WCR, WCS, ALPHA, NPAR), 3) # soil moisture content at wilting point (pF=4.2)
    smfcf = round(get_wc_MvG(10**2.0, WCR, WCS, ALPHA, NPAR), 3) # soil moisture content at field capacity (pF=2.0)
    sm0 = round(WCS, 3) # soil moisture content at saturation
    k0 = round(KSAT, 3) # hydraulic conductivity at saturation
    
    return smtab, contab, smw, smfcf, sm0, k0