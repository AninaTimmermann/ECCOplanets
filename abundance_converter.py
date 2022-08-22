#%% Package and Data Loading
import sys
import os
from pathlib import Path

path = os.path.dirname(os.path.realpath(__file__))

sys.path.append(path)



from pltStyle import *
import data
import classes 
import simulation



import pandas as pd
import numpy as np
from copy import deepcopy

#%% 

def dex2cosmo(sys_name,abundance_dict,sol_ref='Solar'):
    abundance_cosmo = deepcopy(abundance_dict)
    if type(sol_ref) == str:
        sol = data.abundance.loc[sol_ref]
        SiH_sol = sol['Si']/sol['H']
    
    H_star = 1e6 / 10**abundance_dict['Si'] / SiH_sol
    for e in abundance_dict:
        if e in sol:
            abundance_cosmo[e] = int(sol[e]/sol['H'] * 10**abundance_dict[e] * H_star)
        else: del abundance_cosmo[e]
    abundance_cosmo['H'] = int(H_star)
    abundance_cosmo['He'] = int(H_star * 0.08)
    abundance_cosmo['Si'] = 1e6
    abundance_cosmoDF = pd.DataFrame(abundance_cosmo, index=[sys_name])
    return abundance_cosmoDF



def dex2cosmo_nonH(sys_name,abundance_dict,norm_el,sol_ref='Solar'):
    abundance_cosmo = deepcopy(abundance_dict)
    if type(sol_ref) == str:
        sol = data.abundance.loc[sol_ref]
        SiNorm_sol = sol['Si']/sol[norm_el]
    
    Norm_star = 1e6 / 10**abundance_dict['Si'] / SiNorm_sol
    for e in abundance_dict:
        if e in sol:
            abundance_cosmo[e] = int(sol[e]/sol[norm_el] * 10**abundance_dict[e] * Norm_star)
        else: del abundance_cosmo[e]
    abundance_cosmo[norm_el] = int(Norm_star)
    abundance_cosmo['He'] = int(abundance_cosmo['H'] * 0.08)
    abundance_cosmo['Si'] = 1e6
    abundance_cosmoDF = pd.DataFrame(abundance_cosmo, index=[sys_name])
    return abundance_cosmoDF

        
def astro2cosmo(sys_name,abundance_dict):
    abundance_cosmo = deepcopy(abundance_dict)
    Si = abundance_dict['Si']
    x = 1e6/10**Si
    sol = data.abundance.loc['Solar']
    for e in abundance_dict:
        if e in sol:
            p = abundance_dict[e]
            abundance_cosmo[e] = int(x*10**p)
        else: del abundance_cosmo[e]
    
    if not 'He' in abundance_dict:
        abundance_cosmo['He'] = int(abundance_cosmo['H'] * 0.08)
    
    abundance_cosmoDF = pd.DataFrame(abundance_cosmo, index=[sys_name])
    
    return abundance_cosmoDF
    




