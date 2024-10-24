import numpy as np
# import sys, os
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../cantera/build/python')))
import cantera as ct
import yaml

def get_YAML_kTP(self,reaction,collider):
    gas = ct.Solution(yaml=yaml.safe_dump(self.data))
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
    result = 2. * x * first_cheby_poly(x, 1) - first_cheby_poly(x, 0)
    m = 0
    while n - m > 2:
        result = 2. * x * result - first_cheby_poly(x, m+1)
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
    k_TP = get_YAML_kTP(reaction,collider)
    cheb_mat = np.zeros((len(k_TP.flatten()), self.n_T * self.n_P))
    for m, T in enumerate(self.T_ls):
        for n, P in enumerate(self.P_ls):
            for i in range(self.n_T):
                for j in range(self.n_P):
                    T_tilde = reduced_T(T)
                    P_tilde = reduced_P(P)
                    T_cheb = first_cheby_poly(T_tilde, i)
                    P_cheb = first_cheby_poly(P_tilde, j)
                    cheb_mat[m * len(self.P_ls) + n, i * self.n_P + j] = T_cheb * P_cheb
    coef = np.linalg.lstsq(cheb_mat, np.log10(k_TP.flatten()),rcond=None)[0].reshape((self.n_T, self.n_P))
    return coef
def get_cheb_table(self,reaction,collider,label,epsilon,kTP='off'):
    coef = cheby_poly(reaction,collider)
    if kTP=='on':
        colDict = {
            'name': label,
            'efficiency': epsilon,
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
            'name': label,
            'efficiency': epsilon,
        }
    return colDict

def chebyshev(self,foutName):
    return get_cheb_table