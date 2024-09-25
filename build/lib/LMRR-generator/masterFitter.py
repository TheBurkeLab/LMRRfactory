"""
Class that allows for fitting of rate constants at various temperatures and pressures (k(T,P))
"""
import sys, os
# sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
print(sys.path)
from ext.cantera import cantera as ct
from ext.numpy import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from numpy.polynomial.chebyshev import chebval
# import scipy.optimize
from ext.scipy.scipy.optimize import curve_fit
# from io import StringIO
# import matplotlib as mpl
from ext.pyyaml import yaml
# import re

class masterFitter:
    def __init__(self, T_ls, P_ls, inputFile,n_P=7, n_T=7, M_only=False):
        self.T_ls = T_ls
        self.P_ls = P_ls
        self.n_P=n_P
        self.n_T=n_T
        self.P_min = P_ls[0]
        self.P_max = P_ls[-1]
        self.T_min = T_ls[0]
        self.T_max = T_ls[-1]
        self.M_only=M_only
    
        keyReactions = {}
        collider_input = self.openyaml(inputFile)
        for rxn in collider_input['reactions']:
            newdict={}
            for collider in rxn['colliders']:
                newdict[list(collider.keys())[0]]=collider[list(collider.keys())[0]]
            keyReactions[rxn['equation']]=newdict
        chemical_input = self.openyaml(collider_input['chemical-mechanism'])
        if len(chemical_input['phases'])>1:
            print("Error: multiphase kinetics are not supported.")
        shortMechanism={
            'units': chemical_input['units'],
            'phases': chemical_input['phases'],
            'species': chemical_input['species'],
            'reactions': []
            }
        for i, rxn in enumerate(chemical_input['reactions']):
            if rxn['equation'] in keyReactions.keys():
                colliderList = []
                if rxn['type'] == 'pressure-dependent-Arrhenius': 
                    colliderList.append({
                        'collider': 'M',
                        'eps': {'A': 1, 'b': 0, 'Ea': 0},
                        'rate-constants': rxn['rate-constants'],
                    })
                elif rxn['type'] == 'falloff' and 'Troe' in rxn:
                    colliderList.append({
                        'collider': 'M',
                        'eps': {'A': 1, 'b': 0, 'Ea': 0},
                        'low-P-rate-constant': rxn['low-P-rate-constant'],
                        'high-P-rate-constant': rxn['high-P-rate-constant'],
                        'Troe': rxn['Troe'],
                    })
                for j, col in enumerate(keyReactions[rxn['equation']].keys()):
                    colliderList.append({
                        'collider': col,
                        'eps': keyReactions[rxn['equation']][col]
                    })
                shortMechanism['reactions'].append({
                    'equation': rxn['equation'],
                    'type': 'linear-burke',
                    'collider-list': colliderList
                })
        self.shortMech = yaml.safe_dump(shortMechanism,default_flow_style=None,sort_keys=False, allow_unicode=True)
        self.chemical_input = chemical_input
        self.keyReactions = keyReactions
    
    def openyaml(self,fname):
        with open(fname) as f:
            data = yaml.safe_load(f)
        def fix_no(data):
            if isinstance(data, dict):
                return {k: fix_no(v) for k, v in data.items()}
            elif isinstance(data, list):
                return [fix_no(item) for item in data]
            elif data is False and isinstance(data, bool):
                return 'NO'
            return data
        return fix_no(data)
        
    def get_Xvec(self,reaction):
        Prange = self.P_ls
        Xvec=[]
        for i,P in enumerate(Prange):
            for j,T in enumerate(self.T_ls):
                Xvec.append([P,T])
        Xvec=np.array(Xvec)
        return Xvec.T

    def get_Xdict(self,reaction):
        Prange = self.P_ls
        Xdict={}
        for i,P in enumerate(Prange):
            Xdict[P]=self.T_ls
        return Xdict

    def get_YAML_kTP(self,reaction,collider):
        # gas = ct.Solution("shortMech.yaml")
        gas = ct.Solution(yaml=self.shortMech)
        k_TP = []
        for T in self.T_ls:
            temp = []
            for P in self.P_ls:
                gas.TPX = T, P*ct.one_atm, {collider:1.0}
                val=gas.forward_rate_constants[gas.reaction_equations().index(reaction['equation'])]
                temp.append(val)
            k_TP.append(temp)
        return np.array(k_TP)
    
    def first_cheby_poly(self, x, n):
        '''Generate n-th order Chebyshev ploynominals of first kind.'''
        if n == 0: return 1
        elif n == 1: return x
        result = 2. * x * self.first_cheby_poly(x, 1) - self.first_cheby_poly(x, 0)
        m = 0
        while n - m > 2:
            result = 2. * x * result - self.first_cheby_poly(x, m+1)
            m += 1
        # print(result)
        return result
    def reduced_P(self,P):
        '''Calculate the reduced pressure.'''
        P_tilde = 2. * np.log10(P) - np.log10(self.P_min) - np.log10(self.P_max)
        P_tilde /= (np.log10(self.P_max) - np.log10(self.P_min))
        return P_tilde
    def reduced_T(self,T):
        '''Calculate the reduced temperature.'''
        T_tilde = 2. * T ** (-1) - self.T_min ** (-1) - self.T_max ** (-1)
        T_tilde /= (self.T_max ** (-1) - self.T_min ** (-1))
        return T_tilde
    def cheby_poly(self,reaction,collider):
        '''Fit the Chebyshev polynominals to rate constants.
            Input rate constants vector k should be arranged based on pressure.'''
        k_TP = self.get_YAML_kTP(reaction,collider)
        cheb_mat = np.zeros((len(k_TP.flatten()), self.n_T * self.n_P))
        for m, T in enumerate(self.T_ls):
            for n, P in enumerate(self.P_ls):
                for i in range(self.n_T):
                    for j in range(self.n_P):
                        T_tilde = self.reduced_T(T)
                        P_tilde = self.reduced_P(P)
                        T_cheb = self.first_cheby_poly(T_tilde, i)
                        P_cheb = self.first_cheby_poly(P_tilde, j)
                        cheb_mat[m * len(self.P_ls) + n, i * self.n_P + j] = T_cheb * P_cheb
        coef = np.linalg.lstsq(cheb_mat, np.log10(k_TP.flatten()),rcond=None)[0].reshape((self.n_T, self.n_P))
        return coef
    def get_cheb_table(self,reaction,collider,label,epsilon,kTP='off'):
        coef = self.cheby_poly(reaction,collider)
        if kTP=='on':
            colDict = {
                'collider': label,
                'eps': epsilon,
                'temperature-range': [float(self.T_min), float(self.T_max)],
                'pressure-range': [f'{self.P_min:.3e} atm', f'{self.P_max:.3e} atm'],
                'data': []
            }
            for i in range(len(coef)):
                row=[]
                for j in range(len(coef[0])):
                    # row.append(f'{coef[i,j]:.4e}')
                    row.append(float(coef[i,j]))
                colDict['data'].append(row)
        else:
            colDict = {
                'collider': label,
                'eps': epsilon,
            }
        
        return colDict
    
    def get_PLOG_table(self,reaction,collider,label,epsilon,kTP='off'):
        if kTP=='on':
            colDict = {
                'collider': label,
                'eps': epsilon,
                'rate-constants': []
            }
            def arrhenius(T, A, n, Ea):
                return np.log(A) + n*np.log(T)+ (-Ea/(1.987*T))
            gas = ct.Solution(yaml=self.shortMech)
            Xdict = self.get_Xdict(reaction)
            for i,P in enumerate(Xdict.keys()):
                k_list = []
                for j,T in enumerate(Xdict[P]):
                    gas.TPX = T, P*ct.one_atm, {collider:1}
                    k_T = gas.forward_rate_constants[gas.reaction_equations().index(reaction['equation'])]
                    k_list.append(k_T)
                k_list=np.array(k_list)
                popt, pcov = curve_fit(arrhenius, self.T_ls, np.log(k_list),maxfev = 2000)
                newDict = {
                    'P': f'{P:.3e} atm',
                    'A': float(popt[0]),
                    'b': float(popt[1]),
                    'Ea': float(popt[2]),
                }
                colDict['rate-constants'].append(newDict)
        else:
            colDict = {
                'collider': label,
                'eps': epsilon,
            }
        
        return colDict

    def get_Troe_table(self,reaction,collider,label,epsilon,kTP='off'):
        def f(X,a0,n0,ea0,ai,ni,eai,Fcent):
            N= 0.75 - 1.27 * np.log10(Fcent) 
            c= -0.4 - 0.67 * np.log10(Fcent)
            d=0.14
            Rcal=1.987
            Rjoule=8.3145
            M = X[0]*ct.one_atm/Rjoule/X[1]/1000000.0
            k0 = a0 * (X[1] ** n0) * np.exp(-ea0 / (Rcal * X[1]))
            ki = ai * (X[1] ** ni) * np.exp(-eai / (Rcal * X[1]))
            logps = np.log10(k0) + np.log10(M) - np.log10(ki)
            den = logps + c
            den = den / (N - d * den)
            den = np.power(den, 2) + 1.0
            logF = np.log10(Fcent) / den
            logk_fit = np.log10(k0) + np.log10(M) + np.log10(ki) + logF - np.log10(ki + k0 * M)
            return logk_fit
        Xdict=self.get_Xdict(reaction)
        gas = ct.Solution(yaml=self.shortMech)
        logk_list=[]
        for i,P in enumerate(Xdict.keys()):
            for j,T in enumerate(Xdict[P]):
                gas.TPX=T,P*ct.one_atm,{collider:1.0}
                k_TP=gas.forward_rate_constants[gas.reaction_equations().index(reaction['equation'])]
                logk_list.append(np.log10(k_TP))
        # NEED TO GENERALIZE THE FOLLOWING LINES
        # if "H + OH (+M)" in reaction:
        k0_g = [4.5300E+21, -1.8100E+00, 4.9870E+02]
        ki_g = [2.5100E+13, 0.234, -114.2]
        # # elif "H + O2 (+M)" in reaction:
        #     k0_g = [6.366e+20, -1.72, 524.8]
        #     ki_g = [4.7e+12,0.44,0.0]
        # # elif "H2O2 (+M)" in reaction:
        #     k0_g = [2.5e+24,-2.3, 4.8749e+04]
        #     ki_g = [2.0e+12,0.9,4.8749e+04]
        # # elif "NH2 (+M)" in reaction:
        #     k0_g = [1.6e+34,-5.49,1987.0]
        #     ki_g = [5.6e+14,-0.414,66.0]
        # # elif "NH3 <=>" in reaction:
        #     k0_g = [2.0e+16, 0.0, 9.315e+04]
        #     ki_g = [9.0e+16, -0.39, 1.103e+05]
        # # elif "HNO" in reaction:
        #     k0_g = [2.4e+14, 0.206, -1550.0]
        #     ki_g = [1.5e+15, -0.41, 0.0]
        guess = k0_g+ki_g+[1]
        bounds = (
                [1e-100, -np.inf, -np.inf, 1e-100, -np.inf, -np.inf, 1e-100],  # Lower bounds
                [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 1]         # Upper bounds
            )
        Xvec=self.get_Xvec(reaction)
        popt, pcov = curve_fit(f,Xvec,logk_list,p0=guess,maxfev=1000000,bounds=bounds)
        a0,n0,ea0,ai,ni,eai=popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]
        def numFmt(val):
            return round(float(val),3)
        if kTP=='on':
            colDict = {
                'collider': label,
                'eps': epsilon,
                'low-P-rate-constant': {'A':numFmt(a0), 'b': numFmt(n0), 'Ea': numFmt(ea0)},
                'high-P-rate-constant': {'A': numFmt(ai), 'b': numFmt(ni), 'Ea': numFmt(eai)},
                'Troe': {'A': round(float(popt[6]),3), 'T3': 1.0e-30, 'T1': 1.0e+30}
            }
        else:
            colDict = {
                'collider': label,
                'eps': epsilon,
            }
        return colDict
    
    
    
    def final_yaml(self,foutName,fit_fxn): # returns PLOG in LMRR YAML format
        newMechanism={
                'units': self.chemical_input['units'],
                'phases': self.chemical_input['phases'],
                'species': self.chemical_input['species'],
                'reactions': []
                }
        sM = yaml.safe_load(self.shortMech)
        for rxn in self.chemical_input['reactions']:
            if rxn['equation'] in self.keyReactions.keys():
                colliderList=[]
                for rxn2 in sM['reactions']:
                    if rxn2['equation']==rxn['equation']:
                        for j, col in enumerate(self.keyReactions[rxn['equation']].keys()):
                            if j == 0:
                                eps = self.keyReactions[rxn['equation']][col]
                                colliderList.append(fit_fxn(rxn2,col,"M",eps,kTP='on'))
                            elif len(list(rxn2['collider-list'][j].keys()))>2:
                                eps = self.keyReactions[rxn['equation']][col]
                                colliderList.append(fit_fxn(rxn2,col,col,eps,kTP='on'))
                            else:
                                eps = self.keyReactions[rxn['equation']][col]
                                colliderList.append(fit_fxn(rxn2,col,col,eps,kTP='off'))
                newMechanism['reactions'].append({
                    'equation': rxn['equation'],
                    'type': 'linear-burke',
                    'collider-list': colliderList
                })
            else:
                newMechanism['reactions'].append(rxn)
        with open(foutName, 'w') as outfile:
            yaml.dump(newMechanism, outfile, default_flow_style=None,sort_keys=False)
    def Troe(self,foutName): # returns PLOG in LMRR YAML format
        self.final_yaml(foutName,self.get_Troe_table)
    def PLOG(self,foutName): # returns PLOG in LMRR YAML format
        self.final_yaml(foutName,self.get_PLOG_table)
    def cheb2D(self,foutName): # returns Chebyshev in LMRR YAML format
        self.final_yaml(foutName,self.get_cheb_table)