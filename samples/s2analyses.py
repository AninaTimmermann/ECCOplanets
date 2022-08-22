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


#%% Sample 2: Analysing a Simulation
'''
Sample Description:
-------------------    
In this sample, we analyse the output of a simulation.

1.) First, we need load a simulation:

2.) Second, we analyse the simulation in various regards:
    2.a) Basic info
        --> simulation parameters
        --> condensation sequence (list)
    2.b) T-progression line plots
        --> of all solids (overview)
        --> of a particular molecule (e.g. MgSiO3)
        --> of all species containing a particular element
    2.c) Condensation T
        --> basic query of cond T of a species
        --> 50% cond T or an element
    2.d) Rocky Planet Composition
        --> elemental composition
        --> molecular composition

3.) Finally, we will improve the simulation result. We will:
    3.a) recompute a section for all species
    3.b) smooth a section of a particular curve
    3.c) remove isolated values
''' 

#%% 1.) Loading the simulation
sim_path = path+r'\simulations'
sim_name = r'\sim_Solar_p0.0001_T2000-300-1_specs55_2022-03-30.pkl'

sim = simulation.load_sim(sim_name, sim_path)


#%% 2.) Analysis
#%%% 2.a) Basic info

# simulation parameters:
print(sim.info)
#or
for i in sim.info:
    print(i,' : ',sim.info[i])
    
# condensation sequence:
print(sim.condensation_sequence())

#%%% 2.b) T-progression line plots

# overview over solid species:
sim.plot_s()

# curve of specific molecule (e.g. MgSiO3):
sim.plot_m('MgSiO3_s')

# curves of all species containing specified element (e.g. Al):
sim.plot_el('Al')
#adjust plot parameters as usual in matplotlib (NB: use abbr. pl instead of plt)
pl.ylim(0,0.0006)

#%%% 2.c) Condensation T

# (see already sim.plot_m() for visual representation)
# print condensation T of species
print(sim.cond_T('MgSiO3_s')[0])

# condensation T of element
sim.condensation_element('Fe')


#%%% 2.d) Rocky Planet Composition

# elemental composition, e.g. at T=1350K, given in atomic-weight-%
sim.s_composition_el(1350,wt=True)

# molecular composition, e.g. at T=1350K, shown in a pie-chart
sim.s_composition_mol(1350,out='bar')

# => for elemental composition as a func of T, move on to s3planet_composition.py



#%% 3.) Result Improvement
#%%% 3.a) recompute a section
'''
Imagine there is numerical glitch between T=428K and T=425K.
We select the last T value that looks ok and the first value after the glitch,
and set them as 'upper_boundary' and 'lower_boundary'.
Usually, a resolution of dt=0.1 will remove any glitches.
'''

# recompute section T = [429,424]K with dT = 0.1K
upper_boundary = 429
lower_boundary = 424
resolution = 0.1

sim_rec = sim.recompute(upper_boundary,lower_boundary,resolution,save=False)

# plot the solid species with the recomputed section to see the improvement
sim_rec.plot_s()


#%%% 3.b) smooth curve

# we are looking at a very minor gas species: CaO
sim.plot_m('CaO_g')

# we want to smooth the section between 1750 and 1500 K
# the noise is all below the curve, thus we use quantiles 0.9 and 1
sim.smooth('CaO_g',1800,1400,0.9,1,update=False)

# we investigate the smoothed curve and decide the quality of the new curve is good enough
# so we update the DataFrame
sim.smooth('CaO_g',1800,1400,0.9,1,update=True)


#%%% 3.b) remove values

# we are still looking at CaO gas, now we remove the values at [2000, 1999, 426] K
sim.remove([2000, 1999, 426],['CaO_g'],save=False)
