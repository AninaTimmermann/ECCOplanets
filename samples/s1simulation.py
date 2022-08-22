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


#%% Sample 1: Running a Simulation
'''
Sample Description:
-------------------    
In this sample, a new condensation simulation is run.

1.) First, we need to specify the parameters of the simulation:
    --> T_start
    --> T_end
    --> T_int
    --> p
    --> system
    --> SoI (species of interest)

2.) Second, we run the simulation
    Running the simulation will take quite a bit of time - the code will automatically
    estimate the time it will take and print the result (very rough estimate!).
    
    The simulation returns the variable "sim_new".
    This is an instance of the condensation_simulation-class.
    It can either be analysed straight away 
        --> e.g. sim_new.plot_s(), sim_new.s_composition_el(1400,True),
    or the saved .pkl-file can later be loaded for analysis.
        --> it is saved in the folder "simulations"
        --> filename: sim_[system]_p[p]_T[T_start]-[T_end]-[T_int]_specs[len(SoI)]_y-m-d.pkl
''' 


# 1.) specify parameters of simulation:
T_start = 2000 # start temperature in K
T_end = 300 # end temperature in K
T_int = 1 # temperature increment in K
p = 1e-4 # Disk pressure in bar
system = 'Solar' # element abundance of which stellar system

SoI = ['He_g','Al_g','Al2O_g','Ca_g','CaO_g','Mg_g','MgO_g','O2_g','H2O_g','H2_g','Na_g','NaOH_g','Si_g','SiO_g','Ti_g','TiO_g','C_g','CH4_g','CO_g','CO2_g','Fe_g','Cr_g','CrO_g','P_g','PO_g','Ni_g','S_g','SO_g','SO2_g','S2_g','HS_g','H2S_g',
       'Al2O3_s','CaAl12O19_s','CaAl4O7_s','Ca2Al2SiO7_s','CaTiO3_s','MgAl2O4_s','CaAl2Si2O8_s','NaAlSi3O8_s','Mg2SiO4_s','MgSiO3_s','CaMgSi2O6_s','C_s','TiC_s','Ti2O3_s','TiO2_s','SiC_s','Fe_s','Fe3O4_s','Cr_s','Cr2FeO4_s','P_s','Ni_s','FeS_s','MgS_s','CaS_s']


# 2.) run simulation:
sim_new = simulation.Temperature_Progression(T_start,T_end,T_int,p,system,SoI,save=True)

