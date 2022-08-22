# ECCOplanets
We provide a python code to simulate the formation of rocky planets in chemical equilibrium (based on a Gibbs free energy minimisation). 
We also provide tools for the analysis of the simulated planet. 
The code includes a database of thermochemical data and a database of stellar abundance patterns.

## Prerequisites
The code is built in python and uses the common packages pandas, numpy and matplotlib.
The samples are written to be read, used and adapted in an IDE like spyder.

required:
```
python
pandas
numpy
matplotlib
```

recommended IDE, like:
```
spyder
```


## Repository Contents
The repository includes the core python files of our ECCOplanets code:
```
abundance_converter.py
classes.py
condensation.py
data.py
pltStyle.py
simulation.py
```

If you just want to use the code to run condensation simulations and analyse the results, you don't have to look at any of these files.
Quick overview:

- abundance_converter.py:
  includes functions to convert between different stellar elemental abundance formats (used in s4adding_data.py)  
  
- classes.py
  contains all information about the three classes used in the code: 
  -- the molecule class: functions to explore the thermochemical databases (mostly shown in sample s0moleculeDB.py)
  -- the simulation class: functions to analyse the condensation simulation (mostly shown in sample s2analyses.py)
  -- the planet class: functions to analyse the simulated rocky planets (mostly shown in sample s3planet_composition.py
  This document is most helpful if you want to check the available functions and how to use them.
  
- condensation.py
  this is a dummy file, just containing the header necessary to integrate the components of the code. You can use it as a starting point for projects with the
  ECCOplanets code. But we recommend starting with the sample files instead.
  
- data.py
  this code only reads the .pkl files containing all relevant data into python
  
- pltStyle.py
  here, we define properties of the pyplot figures. If you don't like our colour scheme, change it here!

- simulation.py
  this is the backbone of our condensation simulation. It contains the mathematical and thermochemical formulation of our code and the parameters of the
  optimisation procedure. If you want to change the scientific assumption of the condensation simulation, do it here!
  
We also have some folders:

```
data
samples
simulations
```

- data: contains the databases for the thermochemical data and the stellar abundance data
- samples: contains commented code snippets to show how to use the code (look below at "Getting started")
- simulations: these include some already run simulations you can check out
  -- el_amounts: this subfolder contains some already simulated rocky planets


## Getting Started
The best way to get to know the functions of the code is by using the samples provided in the "samples" folder.

### s0moleculeDB.py
This is probably not the most interesting sample to explore the code. Only look at this sample if you are interested in thermochemical data we use in the code.
In the sample, you create an instance of the molecule class and explore the functionalities of the class.
You will also look at our different types of data, i.e. tabulated thermochemical data and fitted Shomate data.
Finally, you will plot some thermochemical information of a molecule.

### s1simulation.py
This might be a good place to start! In this sample we introduce the parameters of the condensation simulation, i.e. the temperature range and resolution, 
the molecule selection, and the elemental abundance pattern.

### s2analyses.py
This is the natural second step, after running a simulation. You load a simulation, retrieve its basic parameters, and analyse the simulation results.
This includes
- producing different lines plots, for the development of the amount of molecules as a function of temperature
- analysing the condensation of molecules and elements
- and improving the simulation result by recomputing a section, removing numerical glitches, and smoothing the curve

### s3planet_composition.py
You can also skip s1 and s2 and start at s3. This sample takes you from running a simulation, over improving the result (cropping the gas phase and recomputing a
section), to computing the composition of rocky planets as a function of temperature and analysing these planets, with regard to the condensation temperature of
elements, bulk composition and devolatilisation relative to the star.

### s4adding_data.py
It is, what it says on the tin. This sample shows you how to expand the thermochemical databases and the stellar abundance database.


## Prospects
The idea of this code is to provide a very simplified starting point to get an approximate idea of the variety of planetary compositions based on the variety of
stellar compositions.
This first release only includes abundance data of F, G, and K stars from the Brewer+2016 database. There is absolutely no reason, though, not to include data of
other stellar types, just be sure to consider uncertainty. Also you might need to think about the sensibility of running simulations on very incomplete abundance
data sets (e.g. if Al, Ca, Mg, Si, or O data is not reported).
We only consider pure minerals at the moment, but there are plans to add solid solutions at a later point. 

## Authors
* **Anina Timmermann** - *Initial work* - [AninaTimmermann](https://github.com/AninaTimmermann)
Please report any problems to me!

## License
This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details

## Acknowledgments
The code was part of my PhD project at the Georg-August-University in Göttingen, Germany.
Many thanks to my supervisors Ansgar Reiners, Andreas Pack, and Yutong Shan!
