# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 15:47:16 2020

@author: Kersti Leppä
"""

import numpy as np

# constants
EPS = np.finfo(float).eps
# zero degrees celsius in Kelvin
DEG_TO_KELVIN = 273.15
# VPDB standard for 13C/12C ratio
D13C_STD = 0.01123720
# VSMOW standard for 18O/16O ratio
D18O_STD = 2.0052e-3
# H2O molar density [mol m-3]
H2O_MOLARDENSITY = 55.5e3
# H2 18O diffusivity [m2 s-1]
H2_18O_DIFFYSIVITY = 2.66e-9

def carbon_fractination(ca, cs, ci, cc, An, Rd, tau_c, d13ca, dt, ratio_pinitol,
                        d13c_pinitol, a_b=2.9e-3, a_s=4.4e-3, a_m=1.8e-3,
                        b=29.0e-3, e=-6.0e-3, f=8.0e-3,
                        sugar_sto=1e5, init_d13c=-32.0):
    """
    Solves net photosynthetic discrimination following Wingate et al. (2007)
    and the accumulation of the d13C signal in the needle sugar and WSC pool.

    Args:
        ca (np.array): CO2 mole fraction in ambient air (ppm).
        cs (np.array): CO2 mole fraction at leaf surface (ppm).
        ci (np.array): CO2 mole fraction in the intercellular spaces (ppm).
        cc (np.array): CO2 mole fraction in the chloroplast (ppm).
        An (np.array): net phosynthesis (umol m-2 s-1).
        Rd (np.array): mitochondrial respiration (umol m-2 s-1).
        tau_c (np.array): CO2 compensation point (ppm).
        d13ca (np.array): Carbon isotope composition in CO2 of ambient air.
        dt (np.array): timestep (s).
        ratio_pinitol (float): ratio of pinitol in WSC (sugar+pinitol) (-).
        d13c_pinitol (float): d13C of pinitiol (permil).
        a_b (float, optional): fractionation during diffusion through the
                boundary layer. Defaults to 2.8e-3.
        a_s (float, optional): fractionation during diffusion through the
                stomata. Defaults to 4.4e-3.
        a_m (float, optional): fractionation during transfer through mesophyll.
                Defaults to 1.8e-3.
        b (float, optional): fractionation by Rubisco. Defaults to 29.0e-3.
        e (float, optional): fractionation during mitochondrial respiration.
                Defaults to -6.0e-3.
        f (float, optional): fractionation during photorespiration. Defaults
                to 8.0e-3.
        sugar_sto (float, optional): size of sugar store (umol C m-2).
                Defaults to 1e5.
        init_d13c (float, optional): initial sugar store d13c (permil).
                Defaults to -32.0.

    Returns:
        Delta (np.array): discrimination of net photosynthesis (permil).
        d13c_flux (np.array): d13C of net CO2 flux (permil).
        d13c_sugar (np.array): d13C of needle sugar pool (permil).
        d13c_bulk (np.array): d13C of needle bulk WSC (sugar+pinitol) (permil).
    """

    # carboxylation efﬁciency
    k = (An + Rd)/(cc - tau_c)

    R_ca = (d13ca/1000 + 1) * D13C_STD

    # initialize arrays
    R_sugar = An * 0.0
    Delta = An * 0.0

    # initial state
    R_sugar[0] = (init_d13c/1000 + 1) * D13C_STD

    for i in range(1, len(An)):

        if k[i] * ca[i] - Rd[i] == 0.0:
            Delta[i] = 0.0
        else:
            Delta[i] = 1 / (k[i] * ca[i] - Rd[i]) * (
                (a_b * (ca[i] - cs[i]) / ca[i]
                + a_s * (cs[i] - ci[i]) / ca[i]
                + a_m * (ci[i] - cc[i]) / ca[i]
                + b * cc[i] / ca[i] - f * tau_c[i] / ca[i]) * (k[i] * ca[i])
                + (1 - R_ca[i] / R_sugar[i-1] * (1 + e)) * Rd[i])

        # weird values when An + Rd > 0 and An < 0
        Delta[i] = np.minimum(np.maximum(Delta[i], -0.5), 0.5)

        # sugar discharge
        discharge = An[i]

        # sugar isotope ratio (implicit)
        R_sugar[i] = ((sugar_sto * R_sugar[i-1] +
                       An[i] * dt[i-1] * R_ca[i] / (1 + Delta[i])) / (
                       sugar_sto + discharge * dt[i-1]))

    # in permil
    d13c_sugar = (R_sugar / D13C_STD - 1) * 1000
    d13c_flux = (R_ca / (1 + Delta)/ D13C_STD - 1) * 1000
    d13c_bulk = ((1-ratio_pinitol) * d13c_sugar
                 + ratio_pinitol * d13c_pinitol)

    return 1000 * Delta, d13c_flux, d13c_sugar, d13c_bulk

def oxygen_fractination(E, An, Rd, w_a, gb, gs, d18o_xylem, d18o_vapor, T, dt, P,
                        ratio_pinitol, d18o_pinitol, e_kb=19e-3, e_k=28e-3,
                        peclet=True, L_eff=15e-3, e_wc=27e-3,
                        nonsteadystate=False, Vlw=15, f1=1.0,
                        sugar_sto=1e5, init_d18O=24.0):
    """
    Solves oxygen fractionation in leaf water, new assimilates and the
    accumulation of the d18O signal in the needle sugar and WSC pool.

    Args:
        E (np.array): transpiraiton (mol m-2 s-1).
        An (np.array): net phosynthesis (umol m-2 s-1).
        Rd (np.array): mitochondrial respiration (umol m-2 s-1).
        w_a (np.array): mole fraction of water vapor in atmosphere (mol mol-1).
        gb (np.array): boundary-layer conductance for water vapor (mol m-2 s-1).
        gs (np.array): stomatal conductance for water vapor (mol m-2 s-1).
        d18o_xylem (np.array): d18O of water vapor (permil).
        d18o_vapor (np.array): d18O of source water (permil).
        T (np.array): leaf temperature (degC).
        dt (np.array): timestep (s).
        P (np.array): atmospheric pressure (Pa).
        ratio_pinitol (float): ratio of pinitol in WSC (sugar+pinitol) (-).
        d18o_pinitol (float): d18O of pinitiol (permil).
        e_kb (float, optional): fractionation during diffusion of water vapor
                through boundary layer. Defaults to 19e-3.
        e_k (float, optional): fractionation during diffusion of water vapor
                through stomata. Defaults to 28e-3.
        peclet (boolean, optional): if true applies peclet model, else
                two-pool model. Defaults to True.
        L_eff (float, optional): leaf mesophyll effective mixing length (m).
                Defaults to 15e-3.
        e_wc (float or np.array, optional): biochemical fractionation factor.
                Defaults to 27e-3.
        nonsteadystate (boolean, optional): if true applies non-steady state
                for leaf water modeling, else steady-state approach. Defaults
                to False.
        Vlw (float, optional): leaf mesophyll water volume (mol m-2). Defaults to 15.
        f1 (float, optional): parameter of two-pool model (-). Defaults to 1.0.
        sugar_sto (float, optional): size of sugar store (umol C m-2).
                Defaults to 1e5.
        init_d18O (float, optional): initial sugar store d18O (permil).
                Defaults to  24.0.

    Returns:
        d18o_e (np.array): d18o at evaporative sites (permil).
        d18o_lw (np.array): d18o of leaf water (permil).
        d18o_sugar (np.array): d18o of needle sugar pool (permil).
        d18o_bulk (np.array): d18O of needle bulk WSC (sugar+pinitol) (permil).

    """

    esat = 611.0 * np.exp((17.502 * T) / (T + 240.97))  # Pa
    w_i = esat / P

    e_star = np.exp(1137/(T + DEG_TO_KELVIN)**2 - 0.4156 /(T + DEG_TO_KELVIN) - 0.0020667) - 1

    e_kkb = (e_k * gb + e_kb * gs) / (gb + gs)

    # isotopic ratios in xylem and ambient vapor
    R_x = (d18o_xylem/1000 + 1) * D18O_STD
    R_v = (d18o_vapor/1000 + 1) * D18O_STD

    # isotopic ratio at evaporative sites
    R_e = (1 + e_star)*((1+e_kkb)*R_x*(1 - w_a / w_i) + R_v *(w_a / w_i))  # w_a/w_i = RH / 100

    # peclet effect
    if peclet:
        p = np.maximum(EPS, L_eff * E / (H2O_MOLARDENSITY * H2_18O_DIFFYSIVITY))
        f1 = (1 - np.exp(-p)) / p
    else:
        # no peclet effec, two pool model
        f1 = An * 0.0 + f1

    R_lw_ss = (R_e - R_x) * f1 + R_x

    # initialize array
    R_lw = R_lw_ss + 0

    if nonsteadystate:
        a = dt / Vlw * E[1:] * w_i[1:] / (w_i[1:] - w_a[1:]) / ((1 + e_star[1:])*(1+e_kkb[1:])) * 1 / f1[1:]
        for i in range(len(R_lw)-1):
            R_lw[i+1] = (a[i] * R_lw_ss[i+1] + R_lw[i]) / (1 + a[i])
        R_e = (R_lw - R_x) / f1 + R_x

    R_sugar = An * 0.0

    # initial state
    R_sugar[0] = (init_d18O/1000 + 1) * D18O_STD

    R_CO2 = (1 + e_wc) * R_lw

    for i in range(1,len(An)):

        # sugar discharge
        discharge = An[i]

        # sugar isotope ratio (implicit)
        R_sugar[i] = ((sugar_sto* R_sugar[i-1] +
                        (An[i] + Rd[i]) * dt[i-1] * R_CO2[i]) / (
                            sugar_sto + (Rd[i] + discharge) * dt[i-1]))

    d18o_e = (R_e / D18O_STD - 1) * 1000
    d18o_lw = (R_lw / D18O_STD - 1) * 1000
    d18o_sugar = (R_sugar / D18O_STD - 1) * 1000
    d18o_bulk = ((1-ratio_pinitol) * d18o_sugar
                 + ratio_pinitol * d18o_pinitol)

    return d18o_e, d18o_lw, d18o_sugar, d18o_bulk

def sourcewater_d18O(df, layerdepth=0.2, d18o_init=-12.):
    """
    Models oxygen isotope composition of source water based on
    a mass balance approach for the soil rooting zone using daily timestep.

    Args:
        df (pd.dataframe): 'ET': evapotranspiration (mm d-1)
                           'Prec': precipitation (mm d-1)
                           'd18O_prec': d18O of precipition (permil)
                           'Ws': volumetric top soi moisture (m3 m-3)
                           'Tair': air temperature (degC)
        layerdepth (float, optional): layer depth (m). Defaults to 0.2.
        d18o_init (float, optional): initial d18o of source water (permil).
                Defaults to -12.

    Returns:
        d18O_sw (np.array): source water d18O (permil).

    """
    # only precipitation during unfrozen conditions considered
    df['Prec'] = np.where(df['Tair'] > 1, df['Prec'], 0.0)

    # resample to daily
    WB_daily = df[['ET','Prec','d18O_prec','Ws','Tair']].resample('D').mean()

    # oxygen isotope ratio of rain
    WB_daily['R_rain'] = (WB_daily['d18O_prec']/1000 + 1) * D18O_STD

    # water storage and its change based on measured soil moisture
    WB_daily['Wsto'] = WB_daily['Ws'] * layerdepth * 1e3  # [mm]
    WB_daily['dWsto'] = 0.0
    WB_daily['dWsto'][:-1] = WB_daily['Wsto'].values[1:] - WB_daily['Wsto'].values[:-1]

    # drainage solved from rooting zone water budget
    WB_daily['Drainage'] = WB_daily['Prec'] - WB_daily['ET'] - WB_daily['dWsto']

    # oxygen isotope ratio in rooting zone solved from the budget of
    # the rooting zone water 18O
    WB_daily['R_sw'] = 0.0
    WB_daily['R_sw'][0] = (d18o_init/1000 + 1) * D18O_STD #
    for i in range(len(WB_daily)-1):
        WB_daily['R_sw'][i+1] = ((WB_daily['Prec'][i] * WB_daily['R_rain'][i] + WB_daily['R_sw'][i] * WB_daily['Wsto'][i]) /
                                 (WB_daily['ET'][i] + WB_daily['Drainage'][i] + WB_daily['Wsto'][i+1]))

    # delta value and to original timescale
    WB_daily['d18O_sw'] = (WB_daily['R_sw']/D18O_STD - 1)*1000
    df.loc[:,'d18O_sw'] = WB_daily['d18O_sw']
    df.loc[:,'d18O_sw'] = df['d18O_sw'].ffill()

    return df['d18O_sw'].values