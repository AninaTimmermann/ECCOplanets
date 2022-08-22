import sys
import os
from pathlib import Path

full_path = os.path.dirname(os.path.realpath(__file__))
path = str(Path(full_path).parents[0])

sys.path.append(path)

from pltStyle import *
import data
import classes 

#%% Sample 0: The Molecule DataBase
'''
Sample Description:
-------------------    
In this sample, we analyse the molecule data.

1.) First we create a new instance of the molecule class and retrieve its basic attributes.
    
2.) Second, we analyse its thermochemical data
    2.a) tabulated data
    2.b) continuous data (based on Shomate eq.)
    2.c) plotting data
''' 

#%% 1.) Creating a Molecule Instance and Retrieving Basic Info

# creating a new molecule class instance for the MgSiO3
mol = classes.molecule('MgSiO3_s')

# retrieving the basic attributes of the molecule
print(mol.name, mol.phase, mol.source, mol.composition)

#%% 2.) Thermochemical Data
#%%% 2.a) tabulated data

# investgating the tabulated data
tab = mol.nist

print(tab)

print(tab[['Temperature','Entropy']])

# investigating the parameters of the Shomate function
shom = mol.shomate

print(shom)

#%%% 2.b) continuous data (func of T)

# retrieving the (continous) values of entropy, enthalpy and gibbs free energy at T=1505.5 K
print(mol.S(1505.5))

print(mol.H_H0(1505.5))

print(mol.G0(1505.5))

#%%% 2.c) plotting data

"""
plot(therm_prop,T_l,T_u,new_fig=True,tab=False)

Parameters
------------
therm_prop : termochemical property 
    -> 'H_H0' : Enthalpy
    -> 'S' : Entropy
    -> 'G0' : Gibbs Energy
    -> 'all' : plot all three in subfigures
T_l : lower Temperature bound
T_u : upper Temperature bound
new_fig (opt) : set False to plot in previous figure
tab (opt) : set True to compare polynomial to tabulated data (only for H_H0 and S)
"""   

mol.plot('all',300,2000,tab=True)



