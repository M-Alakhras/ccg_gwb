# ATNF_utilities.py
"""
Utility functions to query ATNF catalouge.
"""

import os
import pickle
from psrqpy import QueryATNF
from tqdm.auto import tqdm
from ccg_gwb import CCG_CACHEDIR, CCG_PATH, CCG_ENV
from ccg_gwb.simulation.timing_model_parameters import validate_parameters, Parameter

def query_ATNF(condition=None):
    if condition is None:
        condition = ''
    query = QueryATNF(condition=condition)
    return query

def ATNF_ephemeris(query):
    CACHE_file = CCG_CACHEDIR + '/ephems.pkl'
    if os.path.exists(CACHE_file):
        with open(CACHE_file, 'rb') as file:
            ephems = pickle.load(file)
        return ephems
    psrs = query.get_pulsars()
    ephems = []
    for psr in tqdm(psrs):
        try:
            ephems.append(query.get_ephemeris(psr))
        except:
            pass
    with open(CACHE_file, 'wb') as file:
        pickle.dump(ephems, file)
    return ephems

def ATNF2PINT(param):
    pint_param = param
    if param == 'NAME':
        pint_param = 'PSR'
    if param == 'ECCDOT':
        pint_param = 'EDOT'
    if param == 'EDOT':
        pint_param = '_EDOT'
    return pint_param

def parse_ephem(ephem):
    lines = ephem.split('\n')
    all_params = []
    for line in lines:
        if len(line) < 1:
            continue
        param = line.split()
        name = param[0]
        value = param[1]
        error = None
        if len(param) > 2:
            error = param[2]
        param = Parameter(ATNF2PINT(name), value=value, error=error)
        all_params.append(param.auto_detect())
    pint_params, extra_params = validate_parameters(all_params)
    return pint_params