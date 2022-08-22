import sys
import os
from pathlib import Path

full_path = os.path.dirname(os.path.realpath(__file__))
path = str(Path(full_path).parents[0])

sys.path.append(path)

from pltStyle import *
import data
import classes 
import simulation
from abundance_converter import dex2cosmo, astro2cosmo, dex2cosmo_nonH

import numpy as np
import pandas as pd



#%%

molecules = data.molecules
abundance = data.abundance
shomate = data.shomate
nist = data.nist
atom_wt = data.atom_wt


#%% adding new stellar abundance data
#%%% new abundance data (dict)
'''
Kepler 100 data (Schuler+2015, table 3)
in dex ([X/H] = log(X/H)_star - log(X/H)_sun)
'''
Kepler100_dex = {'C':	0.11,
                'O':	0.19,
                'Na':	0.19,
                'Mg':	0.13,
                'Al':	0.16,
                'Si':	0.13,
                'S':	0.11,
                'K':	0.21,
                'Ca':	0.10,
                'Sc':	-0.01,
                'Ti':	0.08,
                'V':	0.07,
                'Cr':	0.09,
                'Mn':	0.06,
                'Fe':	0.06,
                'Co':	0.07,
                'Ni':	0.10,
                'Zn':	0.18}

'''
Solar data (Brewer+2016, table 4)
in astro-log-scaling (log(H) := 12, log(X) = log(n_X/n_H) + 12)
'''
Solar_Brewer_astro = {'H':   12,
                    'C':     8.39,
                    'N':     7.78,
                    'O':     8.66,
                    'Na':    6.17,
                    'Mg':    7.53,
                    'Al':    6.37,
                    'Si':    7.51,
                    'Ca':    6.31,
                    'Ti':    4.90,
                    'V':     4.0,
                    'Cr':    5.64,
                    'Mn':    5.39,
                    'Fe':    7.45,
                    'Ni':    6.23,
                    'Y':     2.21}


#%%% convert to cosmochemical scale (N(Si)=10^6)

'''
for dex2cosmo specify:
    sys_name (str): name of star
    new data (dict): dictionary of new data in dex
    sol_ref (str): sys_name from abundance data frame to used as reference for 
                    conversion; default is Solar (Lodders2003)
'''
Kepler100_conv = dex2cosmo('Kepler 100',Kepler100_dex,sol_ref='Solar')


'''
for dex2cosmo specify:
    sys_name (str): name of star
    new data (dict): dictionary of new data in astro log-scaling
'''
Solar_Brewer_conv = astro2cosmo('Solar_Brewer',Solar_Brewer_astro)


#%%% append to abundance data frame and save

data.save_new_data('abundance',abundance,Kepler100_conv)

data.save_new_data('abundance',abundance,Solar_Brewer_conv)


#%% extend thermochemical database
#%%% nist

# new tabulated data, e.g. NIST-JANAF data
thermdata_path = r'D:\PhD\Condensation\NIST_JANAF\rock'
thermdata_file= r'\Al2S3_s.txt'

# importing T, S, H-H0, DfH, and adding molecule name
new_tab = pd.read_csv(thermdata_path+thermdata_file,sep='\t',header=0,skiprows=1,usecols=[0,2,4,5],names=['Temperature','Entropy','Enthalpy','DfH'])
new_tab = new_tab.assign(molecule_name=len(new_tab)*['Al2S3_s'])

# add to nist DataFrame
nist = pd.concat([nist,new_tab])

# save updated DataFrame
nist.to_pickle(path+r'\data\nist.pkl')

#%%% molecules

# specifying information of molecule
new_mol = {'molecule_name': 'Al2S3_s',
           'molecule_common_name': 'Aluminum Sulfide',
           'data_source': 'NIST-JANAF',
           'molecule_phase': 's',
           'Al': 2,
           'S': 3}

# adding data to molecule DataFrame
new_mol_DF = pd.DataFrame(new_mol, index=[0]).set_index('molecule_name')
molecules = molecules.append(new_mol_DF)

# save updated DataFrame
molecules.to_pickle(path+r'\data\molecules.pkl')

#%%% Shomate
#%%%% new Shomate fit
from scipy.optimize import curve_fit

mol = 'Al2S3_s'

H298 = nist.loc[(nist['Temperature']==298.15)&(nist['molecule_name']==mol),'DfH'].values[0]
data = nist.loc[(nist['molecule_name']==mol),['Temperature','Entropy','Enthalpy']].iloc[1:]

# primarily fit enthalpy
def fit_enthalpy(T,A,B,C,D,E,F):
    t = T/1000
    return A*t + B*t**2 /2 + C*t**3 /3 + D*t**4 /4 - E/t + F - H298

def fit_G(T,G):
    t = T/1000
    return A*np.log(t) + B*t + C*t**2 /2 + D*t**3 /3 - E/(2*t**2) + G

params,pcov = curve_fit(fit_enthalpy,data['Temperature'].values,data['Enthalpy'].values)
A,B,C,D,E,F = params
popt,pcov = curve_fit(fit_G,data['Temperature'].values,data['Entropy'].values)
params = np.append(params,np.array([popt[0],H298]))


# alt: primarily fit entropy
def fit_entropy(T,A,B,C,D,E,G):
    t = T/1000
    return A*np.log(t) + B*t + C*t**2 /2 + D*t**3 /3 - E/(2*t**2) + G

def fit_F(T,F):
    t = T/1000
    return A*t + B*t**2 /2 + C*t**3 /3 + D*t**4 /4 - E/t + F - H298

params2,pcov = curve_fit(fit_entropy,data['Temperature'].values,data['Entropy'].values)
A,B,C,D,E,G = params2
popt,pcov = curve_fit(fit_F,data['Temperature'].values,data['Enthalpy'].values)
params2 = np.insert(params2,-1,popt[0])
params2 = np.append(params2,H298)


# linear extension to cover T range to 6000K
def lin_ext(T,A,F):
    t = T/1000
    return A*t + F - H298

popt,pcov = curve_fit(lin_ext,data['Temperature'].values[-3:],data['Enthalpy'].values[-3:])
A,B,C,D,E,F = [popt[0],0,0,0,0,popt[1]]
popt,pcov = curve_fit(fit_G,data['Temperature'].values[-3:],data['Entropy'].values[-3:])
params_ext = np.array([A,0,0,0,0,F,popt[0],H298])

T_ext = np.arange(data['Temperature'].values[-1],6100,100)

#%%%% quality control: inspect fit
import matplotlib.pyplot as pl

pl.figure()
pl.plot(data['Temperature'].values,data['Entropy'].values,'d',mfc='b',mec='k',label='S (tab)')
pl.plot(data['Temperature'].values,fit_entropy(data['Temperature'].values, params[0], params[1], params[2], params[3], params[4], params[6]),c='b',label='S (Shomate)')
pl.plot(T_ext,fit_entropy(T_ext, params_ext[0], params_ext[1], params_ext[2], params_ext[3], params_ext[4], params_ext[6]),c='b',alpha=0.6,label='S (lin ext.)')

pl.plot(data['Temperature'].values,data['Enthalpy'].values,'d',mfc='r',mec='k',label='H/1000 (tab)')
pl.plot(data['Temperature'].values,fit_enthalpy(data['Temperature'].values, params[0], params[1], params[2], params[3], params[4], params[5]),c='r',label='H/1000 (Shomate)' )
pl.plot(T_ext,fit_enthalpy(T_ext, params_ext[0], params_ext[1], params_ext[2], params_ext[3], params_ext[4], params_ext[5]),c='r',alpha=0.6,label='H/1000 (lin ext.)' )

pl.legend()
pl.xlabel('Temperature [K]')
pl.ylabel(r'$\frac{J}{K \cdot mol}$')

#%%%% add to DataFrame
shomate_new = {'T_lb': [100,1500],
               'T_ub': [1500,6000],
               'A': [params[0],params_ext[0]],
               'B': [params[1],params_ext[1]],
               'C': [params[2],params_ext[2]],
               'D': [params[3],params_ext[3]],
               'E': [params[4],params_ext[4]],
               'F': [params[5],params_ext[5]],
               'G': [params[6],params_ext[6]],
               'H': [params[7],params_ext[7]],
               'molecule_name': mol}

shomate_new_DF = pd.DataFrame(shomate_new)

shomate = shomate.append(shomate_new_DF)

# save updated DataFrame
shomate.to_pickle(path+r'\data\shomate.pkl')




#%% save tables

#nist.to_pickle(path+r'\data\nist.pkl')

#shomate.to_pickle(path+r'\data\shomate.pkl')

#abundance.to_pickle(path+r'\data\abundance.pkl')