import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
import cantera as ct
import yaml

def get_Xdict(self,reaction):
    Prange = self.P_ls
    Xdict={}
    for i,P in enumerate(Prange):
        Xdict[P]=self.T_ls
    return Xdict

def get_PLOG_table(self,reaction,collider,label,epsilon,kTP='off'):
    if kTP=='on':
        colDict = {
            'name': label,
            'eps': epsilon,
            'rate-constants': []
        }
        def arrhenius(T, A, n, Ea):
            return np.log(A) + n*np.log(T)+ (-Ea/(1.987*T))
        # gas = ct.Solution("shortMech.yaml")
        gas = ct.Solution(yaml=yaml.safe_dump(self.outMech))
        Xdict = self.get_Xdict(reaction)
        for i,P in enumerate(Xdict.keys()):
            k_list = []
            for j,T in enumerate(Xdict[P]):
                gas.TPX = T, P*ct.one_atm, {collider:1}
                k_T = gas.forward_rate_constants[gas.reaction_equations().index(reaction['equation'])]
                k_list.append(k_T)
            k_list=np.array(k_list)
            popt, pcov = curve_fit(arrhenius, self.T_ls, np.log(k_list),maxfev = 30000)
            newDict = {
                'P': f'{P:.3e} atm',
                'A': float(popt[0]),
                'b': float(popt[1]),
                'Ea': float(popt[2]),
            }
            colDict['rate-constants'].append(newDict)
    else:
        colDict = {
            'name': label,
            'eps': epsilon,
        }
    return colDict

def plog():
    return get_PLOG_table