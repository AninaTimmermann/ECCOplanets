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

#%% Sample 3: Getting from Simulation to Planet Composition Analysis
'''
Sample Description:
-------------------    
In this sample, we show how to get all the way from running a simulation to 
analysing the composition of a planet.

1.) First, we run the simulation (alternative: load an existing simulation)

2.) Second, we improve the simulation results:
    2.a) We crop the gas phase
        --> after looking at the T-progression of the solid phases, we cut of
            everything at higher Ts
    2.b) We recompute small sections if necessary
        --> if we see any numerical glitches, we recompute the simulation with
            a higher resolution in that T-range

3.) Third, we compute the elemental composition of a rocky planet as a func of 
    T, the condensation Ts of all elements, and save the info (alternative: 
    load an existing planet)
    
4.) Finally, we analyse the planet data
    4.a) we retrieve the data of the simulation that produced the planet and 
         we look at the condTs of the elements
    4.b) we plot the composition of a planet as a func of T and look at the 
         "devolatilization" pattern
''' 

#%% 1.) Running the simulation

T_start = 2000 # start temperature in K
T_end = 500 # end temperature in K
T_int = 1 # temperature increment in K
p = 1e-4 # Disk pressure in bar
SoI = ['He_g','Al_g','Al2O_g','Ca_g','CaO_g','Mg_g','MgO_g','O2_g','H2O_g','H2_g','Na_g','NaOH_g','Si_g','SiO_g','Ti_g','TiO_g','C_g','CH4_g','CO_g','CO2_g','Fe_g','Cr_g','CrO_g','Ni_g',
       'Al2O3_s','CaAl12O19_s','CaAl4O7_s','Ca2Al2SiO7_s','CaTiO3_s','MgAl2O4_s','CaAl2Si2O8_s','NaAlSi3O8_s','Mg2SiO4_s','MgSiO3_s','CaMgSi2O6_s','C_s','TiC_s','Ti2O3_s','TiO2_s','SiC_s','Fe_s','Fe3O4_s','Cr_s','Cr2FeO4_s','Ni_s']
 
system = 'Solar Brewer'

sim = simulation.Temperature_Progression(T_start,T_end,T_int,p,system,SoI,save=True)

'''
# Alternative: Loading the simulation
sim_name = r'sim_Solar Brewer_p0.0001_T2000-500-1_specs45_2022-08-15.pkl'
sim_path = path+r'\simulations'

sim = simulation.load_sim(sim_name, sim_path)
'''

star = sim.info['system']

sim.plot_s() # overview of solid phases

#%% 2.) Improve simulation result
#%%% 2.a) crop gas phase

sim.crop(1800,update=True,save=True)
sim.plot_s()

#%%% 2.b) recompute small section
T_u = 664
T_l = 659

sim_rec = sim.recompute(T_u,T_l,0.1,save=True)
sim_rec.plot_s()

sim = sim_rec

#%% 3.) Computing elemental make-up of planet, condTs, saving (+ inspection)

planet_comp = sim.rop_el(normalise='Al',wt=False,save=True)

my_planet = classes.planet(planet_comp)

'''
# Alternative: Loading planet
plan_name = r'\Solar Brewer_solids-comp_norm-Al_wt-False.pkl'
plan_path = path+r'\simulations\el_amounts'

my_planet = simulation.load_planet(bc_files_sol,path_bc)
'''

#%% 4.) Analysing Planet
#%%% 4.a) Basic Info

print(my_planet.info)

print(my_planet.condTs)

#%%% 4.b) Plotting

# plot bulk composition as a func of T
my_planet.plot_bc(wt=True)

# plot devolatilization as a func of T
my_planet.plot_devol(condTs=True)




