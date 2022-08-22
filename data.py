import pandas as pd
import numpy as np
import sys
import os

path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(path)

data_path = path
#%%


molecules = pd.read_pickle(path+r'\data\molecules.pkl')
abundance = pd.read_pickle(path+r'\data\abundance.pkl')
shomate = pd.read_pickle(path+r'\data\shomate.pkl')
nist = pd.read_pickle(path+r'\data\nist.pkl')
atom_wt = pd.read_pickle(path+r'\data\atom_wt.pkl')




#%% save tables

def save_new_data(name,old_DataFrame,new_DataFrame):
    DataFrame = old_DataFrame.append(new_DataFrame)
    DataFrame.to_pickle(data_path+r'\data\\'+name+'.pkl')



#molecules.to_pickle(path+r'\data\molecules.pkl')

#nist.to_pickle(path+r'\data\nist.pkl')

#shomate.to_pickle(path+r'\data\shomate.pkl')

#abundance.to_pickle(path+r'\data\abundance.pkl')
