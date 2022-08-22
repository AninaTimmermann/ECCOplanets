import pandas as pd
import numpy as np

from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import linprog

from scipy import constants as const
from numpy import log as ln

from datetime import datetime as dt

import sys
import os

path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(path)

from data import molecules, abundance
from classes import molecule, condensation_simulation, planet


path += r'\simulations'

def Temperature_Progression(T_start,T_end,T_int,p,system,SoI,save=True):
    start_time = dt.now()
    '''data pre-processing'''
    SoI_composition = molecules.loc[SoI].drop(['molecule_common_name','data_source'],axis=1).dropna(axis=1,how='all').fillna(0).sort_values(by=['molecule_phase'])
    EoI = list(SoI_composition.columns)[1:]
    x_g_ind = np.where(SoI_composition['molecule_phase'].values=='g')[0]
    SoI = SoI_composition.index.values
    
    col = ['T']+list(SoI)
    dat = np.zeros((len(np.arange(T_start,T_end,-T_int)),len(col)))
    dat[:,0] = np.round(np.arange(T_start,T_end,-T_int),int(abs(np.floor(np.log10(T_int))))+1)
    comp_prog = pd.DataFrame(dat, columns = col)
    comp_prog = comp_prog.set_index('T')

    
    '''functions and derivatives'''
    def const_T(Temperature,spec_list):
        c = []
        r = const.R*Temperature
        rp = r * ln(p)
        for s in spec_list:
            m = molecule(s)
            mu0 = m.G0(Temperature)
            if m.phase == 'g': 
                c.append(mu0 + rp)
            else:
                c.append(mu0)
        return c, r

    def G_sys(x,Temperature):
        lin = np.array(c).dot(x)
        trans = []
        X = np.sum(x[x_g_ind])
        for k in x_g_ind:
            if x[k] > 0:
                trans.append(r*x[k]*ln(x[k]/X))
        return lin + np.sum(trans)

    def J_G_sys(x,Temperature):
        jac = []
        X = np.sum(x[x_g_ind])
        for k in range(len(x)):
            if k in x_g_ind and x[k] > 0:
                jac.append(r*ln(x[k]/X))
            else:
                jac.append(0)
        return np.array(c) + np.array(jac)
    
    def H_G_sys(x,Temperature):
        hess = np.zeros((len(x),len(x)))
        X = np.sum(x[x_g_ind])
        for k in x_g_ind:
            for l in x_g_ind:
                if k == l and  x[k] != 0:
                    hess[k][l] = 1/x[k] - 1/X
                else:
                    hess[k][l] = 1/X
        return r*hess

    '''constraints and bounds'''
    A = SoI_composition.drop('molecule_phase',axis=1).to_numpy().T
    b =  abundance[EoI].loc[system].values
    
    non_neg = Bounds(0,np.inf)
    num_bal = LinearConstraint(A,b,b)
    
    '''initial guess (linear problem)'''
    c,r = const_T(T_start,SoI)
    x0 = linprog(c, A_eq=A, b_eq=list(b), method='simplex')
    if x0['success'] == True:
        x0 = x0['x']
    else:
        if len(x_g_ind)>0:
            c_g = np.array(c)[x_g_ind]
            SoI_g_composition = SoI_composition[SoI_composition['molecule_phase']=='g']
            A_g = SoI_g_composition.drop('molecule_phase',axis=1).to_numpy().T
            x_g = linprog(c_g, A_eq=A_g, b_eq=list(b), method='simplex')
            if x_g['success'] == True:
                x0 = np.zeros(len(SoI))
                x0[x_g_ind] = x_g['x']
            else:
                return 'Error:Computation of initial composition failed'
        else:
            return 'Error:Computation of initial composition failed'
    
    '''Temperature Progression'''
    for T in comp_prog.index:
        
        c,r = const_T(T,SoI)

        x_opt = minimize(G_sys, x0, args=(T), method='trust-constr', jac=J_G_sys, hess=H_G_sys, bounds=non_neg, constraints=num_bal)
        
        comp_prog.loc[T] = list(x_opt['x'])
        
        x0 = x_opt['x']
        
        if T == comp_prog.index[10]:
            time_10 = dt.now()
            dur = (time_10-start_time)/10 * len(comp_prog.index)
            print('expected end of simulation: ', start_time+dur)
        
    comp_prog.attrs = {'system':system,'pressure':p}
    
    if save:
        name = r'\sim_'+system+'_p'+str(p)+'_T'+str(T_start)+'-'+str(T_end)+'-'+str(T_int)+'_specs'+str(len(SoI))+'_'+dt.now().strftime('%Y-%m-%d')
        comp_prog.to_pickle(path+name+'.pkl') 
        
    return condensation_simulation(comp_prog)
            

def load_sim(file,path):
    sim_DF = pd.read_pickle(path+file)      
    sim = condensation_simulation(sim_DF)
    return sim

#%%

def load_planet(file,path):
    planet_DF = pd.read_pickle(path+file)
    plan = planet(planet_DF)
    return plan
