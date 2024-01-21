#%% Package and Data Loading
import pandas as pd
import numpy as np

from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import linprog
from scipy.stats import norm
from scipy import interpolate


from scipy import constants as const
from numpy import log as ln

from datetime import datetime as dt

import matplotlib.pyplot as pl

import sys
import os

path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(path)

from data import molecules, abundance, nist, shomate, atom_wt
from pltStyle import cc, c_dic


path += r'\simulations'

#%% Molecule Class

class molecule:
    def __init__(self,name):
        self.name = name+', '+molecules['molecule_common_name'].loc[name]
        self.source = molecules['data_source'].loc[name]
        self.phase = molecules['molecule_phase'].loc[name]
        self.shomate = shomate.loc[shomate['molecule_name'] == name].drop(['molecule_name'],axis=1)
        self.nist = nist.loc[nist['molecule_name'] == name].drop(['molecule_name'],axis=1)
        self.composition = molecules.loc[name].drop(['molecule_common_name','data_source','molecule_phase']).dropna(axis=0) 
    
    def H_H0(self,T):
        if T > self.shomate['T_ub'].max() or T < self.shomate['T_lb'].min():
            return 'Error: Select T between '+str(self.shomate['T_lb'].min())+' and '+str(self.shomate['T_ub'].max())
        for i in range(len(self.shomate)):
            if T >= self.shomate.iloc[i]['T_lb'] and T <= self.shomate.iloc[i]['T_ub']:
                T_l,T_u,A,B,C,D,E,F,G,H = self.shomate.iloc[i].values
        t = T/1000
        return A*t + B*t**2 /2 + C*t**3 /3 + D*t**4 /4 - E/t + F - H 
    
    def S(self,T):
        if T > self.shomate['T_ub'].max() or T < self.shomate['T_lb'].min():
            return 'Error: Select T between '+str(self.shomate['T_lb'].min())+' and '+str(self.shomate['T_ub'].max())
        for i in range(len(self.shomate)):
            if T >= self.shomate.iloc[i]['T_lb'] and T <= self.shomate.iloc[i]['T_ub']:
                T_l,T_u,A,B,C,D,E,F,G,H = self.shomate.iloc[i].values
        t = T/1000
        return A*ln(t) + B*t + C*t**2 /2 + D*t**3 /3 - E/(2*t**2) + G
    
    def G0(self,T):
        if T > self.shomate['T_ub'].max() or T < self.shomate['T_lb'].min():
            return 'Error: Select T between '+str(self.shomate['T_lb'].min())+' and '+str(self.shomate['T_ub'].max())
        for i in range(len(self.shomate)):
            if T >= self.shomate.iloc[i]['T_lb'] and T <= self.shomate.iloc[i]['T_ub']:
                T_l,T_u,A,B,C,D,E,F,G,H = self.shomate.iloc[i].values
        t = T/1000
        return 1000*(A*t*(1-np.log(t))-1/2*B*t**2-1/6*C*t**3-1/12*D*t**4 - 1/2* E/t + F - G*t)
    
    def plot(self,therm_prop,T_l,T_u,new_fig=True,tab=False):
        """
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
        
        if T_l < self.shomate['T_lb'].min():
            T_l = self.shomate['T_lb'].min()
        if T_u > self.shomate['T_ub'].max():
            T_u = self.shomate['T_ub'].max()
        T = np.arange(T_l,T_u,10)
        
        if therm_prop == 'all':
            fig, axs = pl.subplots(3, sharex=True)
            #fig.suptitle('Thermochemical properties of '+self.name)
            axs[0].plot(T,[self.S(k) for k in T],label=r'$S^o$ (Shomate)')
            axs[0].plot(self.nist['Temperature'],self.nist['Entropy'],'d',label=r'$S^o$ (tabulated)')
            axs[0].set_ylabel(r'$S^o$ $\left[ \frac{\rm{J}}{\rm{K} \cdot \rm{mol}} \right]$')
            axs[0].legend()
            axs[1].plot(T,[self.H_H0(k) for k in T],label=r'$H-H^o$ (Shomate)')
            axs[1].plot(self.nist['Temperature'],self.nist['Enthalpy'],'d',label=r'$H-H^o$ (tabulated)')
            axs[1].set_ylabel(r'$H-H^o$ $\left[ \frac{\rm{kJ}}{\rm{mol}} \right]$')
            axs[1].legend()
            axs[2].plot(T,[self.G0(k) for k in T],label=r'$G^o$ (Shomate)')
            axs[2].legend()
            axs[2].set_ylabel(r'$G^o$ $\left[ \frac{\rm{J}}{\rm{mol}} \right]$')
            axs[2].set_xlabel('$T$ [K]')
            axs[2].set_xlim(T[0],T[-1])
        else:
            if new_fig:
                pl.figure()
            pl.plot(T,[getattr(self, therm_prop)(k) for k in T],label=getattr(self, 'name'))
            pl.xlabel('$T$ [K]')
            pl.ylabel(therm_prop)
            pl.legend()
            if tab:
                if therm_prop == 'H_H0': p = 'Enthalpy'
                elif therm_prop == 'S': p = 'Entropy'
                else: 
                    return 'no tabulated data for G0 available'
                pl.plot(self.nist['Temperature'],self.nist[p],'d',label='tabulated data')
                pl.legend()
        


#%% simulation class

class condensation_simulation:
    def __init__(self,simulation_result):
        self.result = simulation_result
        self.sim_mols = simulation_result.columns.to_list()
        self.info = {'T_start':self.result.index.max(),
                     'T_end':self.result.index.min(),
                     'T_int':self.result.index[0]-self.result.index[1],
                     'molecules':self.sim_mols,
                     'system':self.result.attrs['system'],
                     'disk pressure':self.result.attrs['pressure']}    

    def create_filename(self):
        name0 = r'\sim_'+self.info['system']
        name_p = r'_p'+str(self.info['disk pressure'])
        name_T = r'_T'+str(int(self.info['T_start']))+'-'+str(int(self.info['T_end']-1))+'-'+str(self.info['T_int'])
        name_s = r'_specs'+str(len(self.info['molecules']))
        name_d = r'_'+dt.now().strftime('%Y-%m-%d')
        name = name0+name_p+name_T+name_s+name_d
        return name
    
    def crop(self,temp,update=False,save=False):
        simDF = self.result.copy()
        newDF = simDF.loc[temp:]
        if save:
            self.result = newDF
            self.info['T_start'] = newDF.index.max()
            simDF.attrs = {'system':self.info['system'],'pressure':self.info['disk pressure']}
            name = self.create_filename()
            simDF.to_pickle(path+name+'.pkl')
        
        if update:
            self.result = newDF
            self.info['T_start'] = newDF.index.max()
        
        else: 
            solids = []
            for m in sorted(newDF.columns):
                if molecule(m).phase == 's':
                    solids.append(m)
            pl.figure(figsize=(12,12))
            for s in solids:
                pl.plot(newDF.index.values,newDF[s].values/newDF.sum(axis=1) * 100,label=s)
            pl.xlim(newDF.index[0],newDF.index[-1])
            pl.xlabel('$T$ [K]')
            pl.ylabel('mol-%')
            pl.legend(ncol= 4,bbox_to_anchor=[0.5,-0.5],loc='lower center')
        
    
    def remove(self,temps,spec,update=False,save=False):
        """
        Paramenters
        ------------
        temps : (ls) temperatures at which the value should be removed
        specs : (ls) species for which the values should be removed
        save : (bool) if True, updated DataFrame will be saved under new name
        """   
        data = self.result.copy()
        for t in temps:
            data.at[t,spec] = np.nan
        data[spec] = data[spec].interpolate()
        pl.figure()
        pl.plot(data.index,data[spec]/self.result.sum(axis=1) * 100,ls='-', label=spec+' new')
        pl.plot(self.result.index,self.result[spec]/self.result.sum(axis=1) * 100,ls=':',alpha=0.3, label=spec+' old')
        pl.xlim(np.max(temps)+100,np.min(temps)-100)
        pl.xlabel('$T$ [K]')
        pl.ylabel('mol-%')
        pl.yscale('log')
        pl.legend()
        if update:
            self.result.update(data)
        if save:
            data.attrs = {'system':self.info['system'],'pressure':self.info['disk pressure']}
            name = self.create_filename()
            data.to_pickle(path+name+'.pkl') 
    
    def recompute(self,T_start,T_end,T_int,save=False):
        """
        Paramenters
        ------------
        T_start : (float) highest temperature of the recomputation range
        T_end : (float) lowest temperature of the recomputation range
        T_int : (float) temperature increment for recomputation - should be ~order of magnitude lower than original resolution
        save : (bool) if True, updated DataFrame will be saved under new name
        """  
        p = self.info['disk pressure'] # Disk pressure in bar
        system = self.info['system'] # element abundance of which stellar system
        SoI = self.info['molecules']
        
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
        
        '''initial guess'''
        x0 = self.result[SoI].loc[T_start].values
                
        '''Temperature Progression'''
        for T in comp_prog.index:
            
            c,r = const_T(T,SoI)
            x_opt = minimize(G_sys, x0, args=(T), method='trust-constr', jac=J_G_sys, hess=H_G_sys, bounds=non_neg, constraints=num_bal)
            comp_prog.loc[T] = list(x_opt['x'])
            x0 = x_opt['x']
        
        sim_rec = self.result.copy()
        sim_rec.update(comp_prog)
        sim_rec.attrs = {'system':self.info['system'],'pressure':self.info['disk pressure']}
        
        if save:
            name = self.create_filename()
            sim_rec.to_pickle(path+name+'.pkl') 
        
        return condensation_simulation(sim_rec)
    
    def smooth(self,molecule,T_start,T_end,quant1,quant2,win=20,update=False,save=False):
        """
        Paramenters
        ------------
        molecule : (str) species of which the curve is to be smoothed
        T_start : (float) highest temperature of the smoothing range
        T_end : (float) lowest temperature of the smoothing range
        quant1 : (float) value between 0 and 1, lower quantile of values to be removed
        quant2 : (float) value between 0 and 1, upper quantile of values to be removed
        update : (bool) if True, original DataFrame will be updated (not saved) for analysis
        save : (bool) if True, updated DataFrame will be saved under new name
        """  
        simDF = self.result.copy()
        dat = simDF[molecule].loc[T_start:T_end]
        T = dat.index.values
        newdat = []
        for i in range(int(len(dat)/win)):
            ds = dat.values[i*win:i*win+win]
            q1,q2 = np.quantile(ds,[quant1,quant2])
            for d in ds:
                if d >= q1 and d <= q2:
                    newdat.append(d)
                else:
                    newdat.append(np.nan)
       
        newDF = pd.DataFrame(data={molecule:newdat,'T':T[:len(newdat)]})
        newDF = newDF.set_index('T')
        newDF = newDF.interpolate(axis=0)
        newDF = newDF.dropna()
        pl.figure()
        pl.plot(newDF.index,newDF[molecule]/self.result.loc[newDF.index].sum(axis=1) * 100,ls='-',c='k', label=molecule+' new')
        pl.plot(self.result.index,self.result[molecule]/self.result.sum(axis=1) * 100,ls='-',c='k',alpha=0.3,label=molecule+' old')
        pl.xlim(6000,300)
        pl.xlabel('$T$ [K]')
        pl.ylabel('mol-%')
        pl.legend()
        
        if save:
            simDF.update(newDF)
            simDF.attrs = {'system':self.info['system'],'pressure':self.info['disk pressure']}
            name = self.create_filename()
            simDF.to_pickle(path+name+'.pkl')
        
        if update:
            self.result.update(newDF)
        
        else: 
            return newDF
    
                   
    def sta_lta(self, x, fx): #auxiliary-function for cond_T
        if abs(np.diff(x)[0]) > 1:
            int1D_fun = interpolate.interp1d(x, fx)
            T_1K = np.arange(x[0],x[-1],-1)
            b = int1D_fun(T_1K)
            nsta = 5
            nlta = 30*nsta
        else:
            DT = abs(np.diff(x)[0])
            b = fx.copy()
            nsta = int(5*1/DT)
            nlta = 30*nsta
        e = abs(np.log10(min(b[b>0])))
        b[b < 0] = 0
        sta = np.cumsum((b*10**e) ** 2)
        sta = np.require(sta, dtype=float)
        lta = sta.copy()
        
        
        sta[nsta:] = sta[nsta:] - sta[:-nsta]
        sta /= nsta
        lta[nlta:] = lta[nlta:] - lta[:-nlta]
        lta /= nlta

        sta[:nlta - 1] = 0

        dtiny = np.finfo(0.0).tiny
        idx = lta < dtiny
        lta[idx] = dtiny
        
        cf = sta / lta
        if abs(np.diff(x)[0]) > 1: cf_series = pd.Series(data=cf,index=T_1K)
        else: cf_series = pd.Series(data=cf,index=x)
        
        return cf_series
        
    def trigger(self, x, fx): #auxiliary-function for cond_T
        cf = self.sta_lta(x, fx)
        DT = cf.index.values[0]-cf.index.values[1]
        th_h = np.mean([30,cf.max()])*0.85
        th_l = th_h*0.2
        
        h_idx = cf.index[cf >= th_h]
        if len(h_idx) == 0:
            return pd.Index([]),pd.Index([])
        else:
            on = [h_idx[0]]
            
            for i,dT in enumerate(h_idx[:-1]-h_idx[1:]):
                if not np.allclose(dT ,DT): on.append(h_idx[i+1])
            off = []
            for o in on:
                if any(cf.loc[o:]<=th_l): off.append(cf.loc[o:].index[cf.loc[o:]<=th_l][0])
                else: off.append(cf.index[-1])
            on_c = []
            off_c = []
            for i in range(len(on)):
                if on[i]-off[i] > 6:
                    on_c.append(on[i])
                    off_c.append(off[i])
            on_c = pd.Index(on_c).drop_duplicates()
            off_c = pd.Index(off_c).drop_duplicates()        
            return on_c, off_c

    def cond_T(self, m): #auxiliary-function for condensation_sequence and plotting
        """
        Paramenters
        ------------
        m : (str) species for which the condensation temperature is to be found
        """  
        fx = self.result[m].values/self.result.sum(axis=1) * 100
        T = self.result.index.values
        dT = T[0]-T[1]
        
        a,b = self.trigger(T,fx)
        if dT > 1:
            int1D_fun = interpolate.interp1d(T, fx)
            x = np.arange(T[0],T[-1],-1)
            fx = int1D_fun(x)
            T = x
        
        dat = pd.Series(data=fx,index=T)
        
        cond_Ts = []
        cond_vals = []
        
        max_Ts = []
        max_vals = []
        
        if len(a) > len(b):
            for i in range(len(a)-1):
                if not (a[i] > b[i] and b[i] > a[i+1]): a = np.delete(a,i)

        if len(a) == 0:
            return [], []
        
        else:
            for i in range(len(a)):
                T_max_amount = dat.loc[a[i]+10*dT:2*b[i]-a[i]].max()
                T_max = dat.loc[a[i]+10*dT:2*b[i]-a[i]].idxmax()
                T_hm_amount = T_max_amount/2
                T_hm = (abs(dat.loc[a[i]:T_max] - T_hm_amount)).idxmin()
                cond_Ts.append(T_hm)
                cond_vals.append(dat.loc[T_hm])
                max_Ts.append(T_max)
                max_vals.append(dat.loc[T_max])
                
            return cond_Ts, cond_vals
    
    
    def condensation_sequence(self):
        seq = {}
        for m in self.sim_mols:
            if molecule(m).phase == 's':
                cTm = self.cond_T(m)[0]
                seq[m] = [cTm]
        condDf = pd.DataFrame.from_dict(seq, orient='index', columns=['condensation_T'])
        return condDf.sort_values(by=['condensation_T'],ascending=False)
    
    def condensation_element(self,element,plot=True):
        """
        Paramenters
        ------------
        element : (str) element for which condensation temperature is to be found
        """  
        el_mols = []
        for m in self.sim_mols:
            if element in molecule(m).composition.index.tolist():
                el_mols.append(m)
        solid_form = np.zeros(len(self.result.index))
        gas_form = np.zeros(len(self.result.index))
        for m in el_mols:
            struc = molecule(m).composition
            amount = self.result[m].to_numpy()
            #amount[amount<0] = 0
            if molecule(m).phase == 's':
                solid_form += amount*struc[element]
            else:
                gas_form += amount*struc[element]
           
        total = gas_form + solid_form
        gas_pc = gas_form/total * 100
        solid_pc = solid_form/total * 100
        
        
        if len(np.where(solid_pc>=50)[0]) > 0:
            cond_T_ind = np.where(solid_pc>=50)[0][0]
            cond_T = self.result.index.values[cond_T_ind]
        
        if plot:
            fig, ax1 = pl.subplots()
            #color = 'black'
            color='navy'
            ax1.set_xlabel('$T$ [K]')
            #ax1.set_ylabel('gas form [%]', color=color)
            ax1.plot(self.result.index.values,gas_pc, color=color,lw=2)
            
            ax1.tick_params(axis='y', labelcolor=color)
            ax1.set_xlim(self.result.index[0],self.result.index[-1])
            ax1.set_ylim(-1,101)
            ax2 = ax1.twinx()  
            color = 'royalblue'
            #ax2.set_ylabel('solid form [%]', color=color)  
            ax2.plot(self.result.index.values,solid_pc, color=color,lw=2)
            ax2.tick_params(axis='y', labelcolor=color)
            ax2.set_ylim(-1,101)
            ax1.set_ylabel(r'${ n_{\rm{Ca \,(g)}} }\;/\;{ n_{\rm{Ca \,(tot)}} }$ [%]', color='k')
            ax2.set_ylabel(r'${ n_{\rm{Ca \,(s)}} }\;/\;{ n_{\rm{Ca \,(tot)}} }$ [%]', color='k')
            if len(np.where(solid_pc>=50)[0]) > 0:
                ax2.plot(cond_T,50,'d',mfc='r',mec='k',ms=10,label=r'$T_{c}$ ('+element+') = '+str(int(cond_T))+' K')
                ax2.legend(loc='center left')   
            else:
                ax2.legend([element+' stays mostly in gas phase'],loc='center left')
        
            ind10 = np.where(solid_pc>=10)[0][0]
            ind90 = np.where(solid_pc>=90)[0][0]
            DTC = self.result.index.values[ind10] - self.result.index.values[ind90]
            print(DTC)
        
        
        else: 
            if len(np.where(solid_pc>=50)[0]) > 0: return cond_T
            else: pass 
    
    def rop_el(self,normalise=False,wt=False,save=True):
        solids = []
        elements = {'T':[]}
        for m in self.sim_mols:
            if molecule(m).phase == 's':
                solids.append(m)

        Tu = min(self.info['T_start'],2500)
        Tl = max(self.info['T_end'],300)
        Ti = self.info['T_int']
        
        Temps = np.round(np.arange(Tu,Tl,-Ti),int(abs(np.floor(np.log10(Ti))))+1)

        for i, T in enumerate(Temps):
            elements['T'].append(T)
            for s in solids:
                struc = molecule(s).composition
                amount = self.result[s].loc[T]
                for e in struc.index:
                    if i == 0 and not e in elements: elements[e] = [0]
                        
                    elif len(elements[e]) <= i: elements[e].append(0)
                                       
                    if wt: elements[e][i] += amount*struc[e]*atom_wt.loc[e].values[0]
                    else: elements[e][i] += amount*struc[e]
        composition = pd.DataFrame(elements)
        composition = composition.set_index('T')
        
        total_solids_ppm = composition.sum(axis=1) / abundance.loc[self.info['system']].sum() * 1e6
        
        if normalise == 'percent': composition = composition.divide(composition.sum(axis=1),axis='index')*100
        elif normalise == False: pass
        else: composition = composition.divide(composition[normalise],axis='index')
        
        

        el_condTs = {}
        for e in composition.columns:
            T = self.condensation_element(e,plot=False)
            el_condTs[e] = T
        
        condT_DF = pd.DataFrame(el_condTs, index=['cond T'])
        composition = pd.concat([condT_DF,composition])
        composition['atoms in solids (total) [ppm]'] = np.insert(total_solids_ppm.values,0,np.nan)        

        composition.attrs = {'system':self.info['system'],'pressure':self.info['disk pressure'],'normalisation':normalise,'wt':wt}

        if save:
            name = self.info['system']+'_p'+str(self.info['disk pressure'])+'_solids-comp_'+'norm-'+str(normalise)+'_wt-'+str(wt)
            composition.to_pickle(path+r'\el_amounts\\'+name+'.pkl')
        
        return composition
      
    
    def s_composition_el(self,T,wt=False):
        """
        Paramenters
        ------------
        T : (float) temperature at which solids-composition shall be shown 
        wt : (bool) set True to show in atomic-weight-% instead of mol-%
        """  
        solids = []
        elements = {}
        for m in self.sim_mols:
            if molecule(m).phase == 's':
                solids.append(m)
        for s in solids:
            struc = molecule(s).composition
            amount = self.result[s].loc[T]
            for e in struc.index:
                if not e in elements:
                    elements[e] = 0
                elements[e] += amount*struc[e]
        elements_df = pd.Series(elements)
        elem_val = elements_df.to_list()
        elem_name = elements_df.index.tolist()
        if wt:
            for i,e in enumerate(elem_name):
                weight = atom_wt.loc[e].values[0]
                elem_val[i] *=  weight
        
        labels_leg = []
        labels_pl = []
        def make_autopct(values):
            def my_autopct(pct):
                if pct >= 1:
                    return '{p:.2f}% '.format(p=pct)
                else:
                    return 
            return my_autopct

        for i,v in enumerate(elem_val):
            labels_leg.append(elem_name[i]+':   '+str(round(v/sum(elem_val)*100,2))+' %')
            if v >= 1:
                labels_pl.append(elem_name[i])
            else:
                labels_pl.append('')
        pl.figure(figsize=(15,9))
        pl.pie(elem_val,colors=[c_dic[key] for key in elem_name],wedgeprops={'linewidth': 2.0, 'edgecolor': 'k'})

        pl.title(self.info['system']+': Constituents of Solids at T = '+str(T)+' K [mol-%]')
        if wt:
            pl.title(self.info['system']+': Constituents of Solids at T = '+str(T)+' K [wt-%]')
        pl.legend(labels_leg,loc="best",bbox_to_anchor=(1.2, 0.8),title='Included Elements')
        pl.tight_layout()
        return (elements_df/elements_df.sum()*100).round(decimals=2)
    
    def s_composition_mol(self,T,out='bar'):
        """
        Paramenters
        ------------
        T : (float) temperature at which solids-composition shall be shown 
        out : (bool) representation of output
            -> 'bar' : bar chart (Default)
            -> 'pie' : pie chart
            -> 'ls' : list
        """ 
        solids = []
        for m in self.sim_mols:
            if molecule(m).phase == 's':
                solids.append(m)
        sol_comp = self.result[solids].loc[T]
        sol_tot = sol_comp.sum()
        sol_comp_pc = sol_comp.div(sol_tot) * 100
        sol_comp_pc_v = np.squeeze(sol_comp_pc.values)
        tot = self.result.loc[T].sum()
        if out == 'pie':
            cs = pl.cm.get_cmap('jet')(np.linspace(0,1,len(sol_comp_pc_v)))
            labels_leg = []
            labels_pl = []
            def make_autopct(values):
                def my_autopct(pct):
                    if pct >= 1:
                        return '{p:.2f}% '.format(p=pct)
                    else:
                        return 
                return my_autopct

            for i,v in enumerate(sol_comp_pc_v):
                labels_leg.append(solids[i]+':   '+str(round(v,2))+' %')
                if v >= 1:
                    labels_pl.append(solids[i])
                else:
                    labels_pl.append('')
            pl.figure(figsize=(15,9))
            pl.pie(sol_comp_pc_v,autopct=make_autopct(sol_comp_pc_v),pctdistance=0.8,colors=cs,labels=labels_pl,rotatelabels=False,wedgeprops={'linewidth': 2.0, 'edgecolor': 'k'})
            pl.title(self.info['system']+': Composition of Solids at T = '+str(T)+' K; total solids: '+str(round(sol_tot/tot *1000,2))+' mol-‰')
            pl.legend(labels_leg,loc="best",bbox_to_anchor=(1.2, 0.8),title='Included Species')
            pl.tight_layout()
        elif out == 'bar':
            pl.figure()
            pl.bar(np.arange(0,len(sol_comp_pc_v),1),sol_comp_pc_v,ec='k',fc='skyblue')
            pl.xticks(np.arange(0,len(sol_comp_pc_v),1), solids, rotation='vertical')
            pl.ylabel('mol-% of total solids')
            pl.title('Composition of Solids at T = '+str(T)+' K')
        else:
            return sol_comp_pc, sol_tot/tot *100
    
    def plot_all(self):
        pl.figure()
        for s in self.sim_mols:
            pl.plot(self.result.index.values,self.result[s].values/self.result.sum(axis=1) * 100,label=s)
        pl.yscale('log')
        pl.xlabel('$T$ [K]')
        pl.ylabel('mol-%')
        #pl.legend()
        pl.xlim(self.result.index[0],self.result.index[-1])
    
    def plot_s(self):
        solids = []
        '''
        for m in sorted(self.sim_mols):
            if molecule(m).phase == 's':
                solids.append(m)
        '''
        pl.figure(figsize=(16,10))
        
        seq = self.condensation_sequence()
        for i in seq.index:
            if not seq['condensation_T'][i] == []:
                solids.append(i)
        for s in solids:
            pl.plot(self.result.index.values,self.result[s].values/self.result.sum(axis=1) * 100,label=s)
        pl.xlabel('$T$ [K]')
        pl.ylabel('mol-%')
        pl.legend(ncol= 4,bbox_to_anchor=[0.5,-0.5],loc='lower center',fontsize=22)
        pl.xlim(1800,self.info['T_end']) 

    
    def plot_m(self,mol,nf=True):
        if nf==True:
            pl.figure()
        pl.plot(self.result.index.values,self.result[mol].values/self.result.sum(axis=1) * 100,label=mol)
        pl.xlim(self.result.index.values[0],self.result.index[-1])
        if molecule(mol).phase == 's':
            c_T, c_mol = self.cond_T(mol)
            pl.plot(c_T,c_mol,'d',mfc='r',mec='k',ms=10,label=r'$T_{cond}$ = '+str(c_T)+' K')
            pl.xlim(2000,self.result.index[-1])
        pl.xlabel('$T$ [K]')
        pl.ylabel('mol-%')
        pl.legend()

    def plot_ls(self,mol_ls,nf=True):
        if nf == True: pl.figure()
        for s in mol_ls:
            pl.plot(self.result.index.values,self.result[s].values/self.result.sum(axis=1) * 100,label=s)
        pl.xlabel('$T$ [K]')
        pl.ylabel('mol-%')
        pl.legend()
        pl.xlim(self.info['T_start'],self.info['T_end'])    
    
    def plot2(self,mol1,mol2):
        fig, ax1 = pl.subplots()
        color = 'black'
        ax1.set_xlabel('$T$ [K]')
        ax1.set_ylabel(mol1+' mol-%', color=color)
        ax1.plot(self.result.index.values, self.result[mol1].values/self.result.sum(axis=1) * 100, color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.set_xlim(self.result.index[0],self.result.index[-1])
        ax2 = ax1.twinx()  
        color = 'royalblue'
        ax2.set_ylabel(mol2+' mol-%', color=color)  
        ax2.plot(self.result.index.values, self.result[mol2].values/self.result.sum(axis=1) * 100, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
    
    def plot_el(self,element,nf=True):
        el_mols = []
        for m in self.sim_mols:
            if element in molecule(m).composition.index.values:
                el_mols.append(m)
        pl.figure()
        for s in el_mols:
            pl.plot(self.result.index.values,self.result[s].values/self.result.sum(axis=1) * 100,label=s)
        pl.xlabel('$T$ [K]')
        pl.ylabel('mol-%')
        pl.legend(ncol= 4,bbox_to_anchor=[0.5,-0.5],loc='lower center', borderaxespad=0.)
        pl.xlim(self.info['T_start'],self.info['T_end'])


#%% planet class
        
class planet:
    def __init__(self,composition_result):
        self.result = composition_result
        self.condTs = self.result.iloc[0]
        self.composition = self.result.iloc[1:].drop(columns=['atoms in solids (total) [ppm]'])
        self.atoms_s_ppm = self.result['atoms in solids (total) [ppm]'].iloc[1:]
        self.info = {'T_start':self.composition.index.max(),
                     'T_end':self.composition.index.min(),
                     'T_int':self.composition.index[0]-self.composition.index[1],
                     'elements':self.composition.columns,
                     'system':self.result.attrs['system'],
                     'disk pressure':self.result.attrs['pressure'],
                     'normalisation':self.result.attrs['normalisation'],
                     'wt':self.result.attrs['wt']}
        
    def feedingzone(self,dT,wt=True):
        els = self.info['elements'].sort_values(ascending=False)
        bc_pc = self.composition.copy(deep=True)
        bc_pc = bc_pc.loc[self.condTs.max():]
        
        if wt and not self.info['wt']:
            for e in els:
                bc_pc[e] *= atom_wt.loc[e].values
        
        bc_pc_fz = bc_pc.rolling(dT,center=True).mean()
        
        tot = bc_pc_fz.sum(axis=1)
        bc_pc_fz = bc_pc_fz.divide(tot,axis=0) * 100
        return bc_pc_fz
    
    
    def plot_bc(self,condTs=False,wt=True):
        els = self.info['elements'].sort_values(ascending=False)
        bc_pc = self.composition.copy(deep=True)
        
        if wt and not self.info['wt']:
            for e in els:
                bc_pc[e] *= atom_wt.loc[e].values
        
        tot = bc_pc.sum(axis=1)
        bc_pc = bc_pc.divide(tot,axis=0) * 100
        
        
        pl.figure(figsize=(12,6))
        for i in range(1,len(els)+1):
            pl.fill_between(np.array(bc_pc.index, dtype=float), np.array(bc_pc[els[:i-1]].sum(axis=1).values, dtype=float),np.array(bc_pc[els[:i]].sum(axis=1), dtype=float), color= c_dic[els[i-1]],alpha=0.8,ec='k',label=els[i-1])
            if condTs:
                if type(self.condTs[els[i-1]]) == np.float64:
                    pl.axvline(self.condTs[els[i-1]], color= c_dic[els[i-1]])
        
        pl.xlim(self.condTs.max()+5,self.info['T_end'])
        pl.ylim(0,100)
        pl.legend(ncols=int(len(els)/2),bbox_to_anchor =(0.5,-0.45), loc='lower center')
        pl.xlabel('$T$ [K]')
        pl.ylabel('bulk composition [wt-%]\n ') 
        
    def plot_bc_fz(self,dT,method='Gauss',condTs=False,wt=True):
        els = self.info['elements'].sort_values(ascending=False)
        bc_pc = self.composition.copy(deep=True)
        bc_pc = bc_pc.loc[self.condTs.max():]
        ylab = 'bulk composition [wt-%]'
                
        if wt and not self.info['wt']:
            for e in els:
                bc_pc[e] *= atom_wt.loc[e].values
        
        if not wt: 'bulk composition [%]'
        
        
        if method == 'mean':
            bc_pc_fz = bc_pc.multiply(self.atoms_s_ppm,axis=0).rolling(dT,center=True).mean()
            tot = bc_pc_fz.sum(axis=1)
            bc_pc_fz = bc_pc_fz.divide(tot,axis=0) * 100
            ylab += '\nblock feeding zone'
            
        if method == 'Gauss':
            bc_pc_fz = bc_pc.copy()
            ylab += '\nGaussian feeding zone'
            for T in bc_pc.index.tolist():
                Gauss = norm.pdf(bc_pc.index.tolist(),T,dT)

                bc_pc_norm = bc_pc.multiply(Gauss,axis=0).multiply(self.atoms_s_ppm,axis=0).sum(axis=0)
                
                bc_pc_fz.loc[T] =  bc_pc_norm / bc_pc_norm.sum() * 100
            dT *= 2
        
        
        pl.figure(figsize=(12,5))
        for i in range(1,len(els)+1):
            pl.fill_between(np.array(bc_pc_fz.index, dtype=float), np.array(bc_pc_fz[els[:i-1]].sum(axis=1).values, dtype=float),np.array(bc_pc_fz[els[:i]].sum(axis=1), dtype=float), color= c_dic[els[i-1]],alpha=0.8,ec='k',label=els[i-1])
            if condTs:
                if type(self.condTs[els[i-1]]) == np.float64:
                    pl.axvline(self.condTs[els[i-1]], color= c_dic[els[i-1]])
        
        pl.xlim(self.condTs.max()-dT/2,self.info['T_end']+dT/2)
        pl.ylim(0,100)
        #pl.legend()
        pl.xlabel(r'T$_{central}$ [K] (feeding zone width: '+str(dT)+'K)')
        pl.ylabel(ylab)
     
                
        
    
    def plot_devol(self,norm='Al',feedingzone=False,condTs=False):
        els = self.info['elements']
        devol = self.composition.copy(deep=True)
        
        if feedingzone:
            devol = devol.rolling(feedingzone,center=True).mean()
        
        stell_original = abundance.loc[self.info['system']]
        stell = stell_original.copy(deep=True)
        if self.info['wt']:
            for e in els:
                stell.at[e] *= atom_wt.loc[e].values
        
        if norm != self.info['normalisation']:
            for e in els:
                devol[e] /= self.composition[norm]
        
        for e in els:
            stell[e] /= stell_original[norm]                
        
        pl.figure()
        for e in els:
            pl.plot(devol.index,devol[e]/stell[e],ls='-',c=c_dic[e],label=e)
            if condTs:
                if type(self.condTs[e]) == np.float64:
                    pl.axvline(self.condTs[e], ls=':', color= c_dic[e])
        
        pl.xlim(self.condTs.max()+5,self.info['T_end'])
        pl.ylim(0,1)
        pl.legend()
        pl.xlabel('$T$ [K]')
        pl.ylabel(r'$\left(\frac{X_i}{X_{Al}}\right)_p / \left(\frac{X_i}{X_{Al}}\right)_\star$') 
        return devol/stell
    
    def plot_el(self,element,condT=True,nf=True):
        if nf: pl.figure()
        pl.plot(self.composition.index,self.composition[element],'-',c='k')
        if condT: pl.axvline(self.condTs[element],ls=':',c='k')
        pl.xlabel('$T$ [K]')
        if self.info['wt']:
            if self.info['normalisation'] == 'percent': ytit1 = '[wt-%]'
            elif not self.info['normalisation']: ytit1 = '[u]'
            else: ytit1 = '[m/m('+self.info['normalisation']+')]'
        elif not self.info['wt']:
            if self.info['normalisation'] == 'percent': ytit1 = '[mol-%]'
            elif not self.info['normalisation']: ytit1 = '[mol]'
            else: ytit1 = '[X/X('+self.info['normalisation']+')]'
        pl.ylabel(element+' '+ytit1)
        pl.xlim(self.condTs.max()+5,self.info['T_end'])
        pl.ylim(0,self.composition[element].max())
        
        
        
        
        

            
        
        