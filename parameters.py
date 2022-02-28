# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 16:09:15 2020

@author: Kersti Leppä
"""

photop = {
        'Vcmax': 32.,  # (umol m-2 s-1)
        'Jmax': 64.,  # (umol m-2 s-1)
        'Rd':0.45,  # (umol m-2 s-1)
        'tresp': { # temperature response parameters
               'Vcmax': [53., 200.0, 640.],
               'Jmax': [38., 200.0, 656.],
              'Rd': [36.]
              },
        # quantum yield parameter
        'alpha': 0.05,  # (mol/mol) -- NOTE, value also accounts for converting PAR (umol m-2(ground) s-1) to incident PAR (umol m-2(leaf) s-1)
        # curvature parameter
        'theta': 0.7,  # (-)
        # boundary layer conductance
        'gb_c': 1.5,  # (mol m-2 s-1)
        'gb_v': 1.6*1.5,  # (mol m-2 s-1)
        # gs-model parameters (Medlyn et al. 2011)
        'g0': 1.0e-3,  # (mol m-2 s-1)
        'g1': 2.0,  # (kPa^0.5)
        'gs_min': 3.0e-3,  # (mol m-2 s-1)
        # gm-model parameters
        'gm0': 0.01,  # (mol m-2 s-1)
        'a2': 4.5,  # (-)
        # Rew-based drought response
        'drp': [0.39, 0.83, 0.31, 3.0]
        }

phenop = {
    'fmin': 0.0,  # (-)
    'Tbase': -3.1,  # (degC)
    'tau': 18.,  # (degC)
    'smax': 17.3-3.1,  # (degC)
    }

