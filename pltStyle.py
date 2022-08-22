import matplotlib.pyplot as pl
from cycler import cycler

cc = (cycler(linestyle=['-', '--', '-.',':'])*
      cycler(color=['k','royalblue','orange','forestgreen','firebrick','hotpink']))

pl.rc('lines',lw=2)
pl.rc('axes', prop_cycle=cc)
pl.rc('font',size=18)
pl.rc('figure',figsize=(12,8))
pl.rc('figure',autolayout=True) 

c_dic = {'Al':'darkgreen',
        'C':'dimgrey',
        'Ca':'darkorange',
        'Fe':'darkred',
        'Mg':'white',
        'O':'lightskyblue',
        'Si':'yellowgreen',
        'Ti':'black',
        'Ni':'hotpink',
        'N':'purple',
        'Na':'peachpuff',
        'Cr':'gold',
        'P':'red',
        'H':'navy',
        'S':'yellow'}