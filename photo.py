# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 15:36:56 2020

@author: Kersti Leppä
"""

import numpy as np

# constants
# zero degrees celsius in Kelvin
DEG_TO_KELVIN = 273.15
# universal gas constant [J mol-1 K-1]
GAS_CONSTANT = 8.314
# O2 concentration in air [umol m-1]
O2_IN_AIR = 2.10e5
# H2O to CO2 diffusivity ratio [-]
H2O_CO2_RATIO = 1.6

TN = 25.0 + DEG_TO_KELVIN  # reference temperature [K]

def An_gs_model(photop, Qp, T, wa, ca, gb_c, gb_v, P=101300.0, fpheno=1.0, Rew=1.0):
    """
    Leaf gas-exchange by Farquhar-Medlyn model with description of mesophyll resistance.

    Args:
        photop (dict): parameter dict
        Qp (np.array): incident PAR at leaves (umol m-2 s-1).
        T (np.array): leaf temperature (degC).
        wa (np.array): mole fraction of water vapor in atmosphere (mol mol-1).
        ca (np.array): ambient CO2 (ppm).
        gb_c (np.array): boundary-layer conductance for CO2 (mol m-2 s-1).
        gb_v (np.array): boundary-layer conductance for H2O (mol m-2 s-1).
        P (np.array, optional): atm. pressure (Pa). Defaults to 101300.0.
        fpheno (np.array, optional): farctor accounting for phenological
        acclimation (-). Defaults to 1.0.
        Rew (np.array, optional): relative extractable water (-). Defaults to 1.0.

    Returns:
        An (np.array): net phosynthesis (umol m-2 s-1).
        Rd (np.array): mitochondrial respiration (umol m-2 s-1).
        Rp (np.array): photorespiration (umol m-2 s-1).
        E (np.array): transpiraiton (mol m-2 s-1).
        gs (np.array): stomatal conductance for CO2 (mol m-2 s-1).
        gm (np.array): mesophyll conductance for CO2 (mol m-2 s-1).
        cs (np.array): CO2 mole fraction at leaf surface (ppm).
        ci (np.array): CO2 mole fraction in the intercellular spaces (ppm).
        cc (np.array): CO2 mole fraction in the chloroplast (ppm).
        Tau_c (np.array): CO2 compensation point (ppm).

    """

    esat = 611.0 * np.exp((17.502 * T) / (T + 240.97))  # Pa
    VPD = 1e-3 * (esat - wa * P)  # kPa
    Tk = T + DEG_TO_KELVIN

    # drought response
    b = photop['drp']
    # g1 decrease with Rew
    fm = np.minimum(1.0, (Rew / b[0])**b[1])
    # apparent Vcmax, Jmax and Rd decrease with Rew
    fv = np.minimum(1.0, (Rew / b[2])**b[3])

    # parameters
    Vcmax = photop['Vcmax'] * fpheno * fv
    Jmax = photop['Jmax'] * fpheno * fv
    Rd = photop['Rd'] * fpheno * fv
    alpha = photop['alpha']
    theta = photop['theta']
    g0 = photop['g0']
    g1 = photop['g1'] * fm
    gmin= photop['gs_min']
    gm0 = photop['gm0']
    a2 = photop['a2']

    # CO2 compensation point (Bernacchi et al 2001)
    Tau_c = 42.75 * np.exp(37830*(Tk - TN) / (TN * GAS_CONSTANT * Tk))

    # Kc & Ko (umol/mol), Rubisco activity for CO2 & O2 (Bernacchi et al 2001)
    Kc = 404.9 * np.exp(79430.0*(Tk - TN) / (TN * GAS_CONSTANT * Tk))
    Ko = 2.784e5 * np.exp(36380.0*(Tk - TN) / (TN * GAS_CONSTANT * Tk))

    # adjust parameters for temperature
    Vcmax_T = photop['tresp']['Vcmax']
    Jmax_T = photop['tresp']['Jmax']
    Rd_T = photop['tresp']['Rd']
    Vcmax, Jmax, Rd = photo_temperature_response(Vcmax, Jmax, Rd, Vcmax_T, Jmax_T, Rd_T, Tk)

    #  model parameters k1_c, k2_c [umol/m2/s]
    Km = Kc*(1.0 + O2_IN_AIR / Ko)
    J = (Jmax + alpha*Qp -((Jmax + alpha*Qp)**2.0 - (4*theta*Jmax*alpha*Qp))**(0.5)) / (2*theta)

    # --- iterative solution for cs, ci and cc ---
    MaxIter = 50
    err = 9999.0
    cnt = 1
    cs = ca  # leaf surface CO2
    cc = 0.8*ca
    while err > 0.01 and cnt < MaxIter:

        # Assimilation rate / apparent photosynthesis (A = Vc - Rp)
        ## rubisco -limited rate
        Av = Vcmax * (cc - Tau_c) / (cc + Km)
        ## RuBP -regeneration limited rate
        Aj = J/4.0 * (cc - Tau_c) / (cc + 2.0*Tau_c)
        ## single limiting rate
        A = np.minimum(Av, Aj)

        # Carboxylation (true photosynthesis)
        Vc = A / (1 - Tau_c / cc)

        # Photorespiration
        Rp = Vc - A

        # Net photosynthesis
        An = Vc - Rp - Rd

        An1 = np.maximum(An, 0.0)

        # stomatal conductance (Medlyn et al. 2011)
        gs = (g0 + (1.0 + g1 / (VPD**0.5)) * An1 / cs)

        gs = np.maximum(gmin, gs)

        # mesophyll conductance (Dewar et al. 2018 Eq. 9)
        gm = np.maximum(gm0, a2 * An1 / (cc - Tau_c))

        # CO2 supply
        cc0 = cc
        cs = np.maximum(ca - An / gb_c, 0.5*ca)  # through boundary
        ci = np.maximum(cs - An / gs, 0.01*ca)  # through stomata
        cc = np.minimum(np.maximum(ci - An / gm, 0.95*cc0),1.05*cc0)  # through mesophyll

        if cnt > MaxIter - 10:
            cc = (cc0 + cc) / 2
        err = max(abs(cc0 - cc))
        cnt += 1

    gs_v = H2O_CO2_RATIO * gs

    geff = (gb_v*gs_v) / (gb_v + gs_v)  # molm-2s-1
    E = geff * VPD / (1e-3 * P)  # leaf transpiration rate

    return An, Rd, Rp, E, gs, gm, cs, ci, cc, Tau_c

def photo_temperature_response(Vcmax0, Jmax0, Rd0, Vcmax_T, Jmax_T, Rd_T, T):
    """
    Adjusts Farquhar model parameters for temperature.

    Args:
        Vcmax0 (float): maximum carboxylation capacity at 25 degC (umol m-2s-1).
        Jmax0 (float): maximum electron transport rate at 25 degC (umol m-2s-1).
        Rd0 (float): mitochondrial respiration rate at 25 degC (umol m-2s-1).
        Vcmax_T (list): temperature response parameters for Vcmax.
        Jmax_T (list): temperature response parameters for Jmax.
        Rd_T (list): temperature response parameters for Rd.
        T (np.array): leaf temperature (degC).

    Returns:
        Vcmax (np.array): Vcmax adjusted for T (umol m-2s-1).
        Jmax (np.array): Jmax  adjusted for T (umol m-2s-1).
        Rd (np.array): Rd adjusted for T (umol m-2s-1).

    """

    # --- Vcmax (Medlyn et al. 2002) ---
    Ha = 1e3 * Vcmax_T[0]  # J mol-1, activation energy Vcmax
    Hd = 1e3 * Vcmax_T[1]  # J mol-1, deactivation energy Vcmax
    Sd = Vcmax_T[2]  # entropy factor J mol-1 K-1

    NOM = np.exp(Ha * (T - TN) / (GAS_CONSTANT*TN*T)) * (1.0 + np.exp((TN*Sd - Hd) / (TN*GAS_CONSTANT)))
    DENOM = (1.0 + np.exp((T*Sd - Hd) / (T*GAS_CONSTANT)))
    Vcmax = Vcmax0 * NOM / DENOM

    # --- Jmax (Medlyn et al. 2002) ---
    Ha = 1e3 * Jmax_T[0]  # J mol-1, activation energy Vcmax
    Hd = 1e3 * Jmax_T[1]  # J mol-1, deactivation energy Vcmax
    Sd = Jmax_T[2]  # entropy factor J mol-1 K-1

    NOM = np.exp(Ha * (T - TN) / (GAS_CONSTANT*TN*T)) * (1.0 + np.exp((TN*Sd - Hd) / (TN*GAS_CONSTANT)))
    DENOM = (1.0 + np.exp((T*Sd - Hd) / (T*GAS_CONSTANT)))
    Jmax = Jmax0*NOM / DENOM

    # --- Rd (Bernacchi et al 2001) ---
    Ha = 1e3 * Rd_T[0]  # J mol-1, activation energy Rd
    Rd = Rd0 * np.exp(Ha*(T - TN) / (TN * GAS_CONSTANT * T))

    return Vcmax, Jmax, Rd

def pheno_cycle(phenop, Tdaily, doy):
    """
    Computes phenology modifier based on temperature acclimation (Mäkelä et al. 2008).

    Args:
        phenop (dict): parameters of acclimation model.
            'tau': time constant (days)
            'Tbase': base temperature (degC)
            'smax': threshold for full acclimation (degC)
            'fmin': minimum photocapacity (-)
        Tdaily (np.array): mean daily air temperature (degC).
        doy (np.array): day of year.

    Returns:
        fpheno (np.array): phenology modifier [fmin...1].

    """

    X = np.zeros(len(Tdaily))
    for k in range(1,len(Tdaily)):
        if doy[k] != doy[k-1]:
            X[k] = X[k - 1] + 1.0 / phenop['tau'] * (Tdaily[k-1] - X[k - 1])
        else:
            X[k] = X[k - 1]
    S = np.maximum(X - phenop['Tbase'], 0.0)

    fpheno = np.maximum(phenop['fmin'],
                        np.minimum(S / (phenop['smax'] - phenop['Tbase']), 1.0))
    return fpheno