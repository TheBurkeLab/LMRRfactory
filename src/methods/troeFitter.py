import numpy as np
from scipy.optimize import curve_fit
# import sys, os
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../cantera/build/python')))
import cantera as ct
import yaml

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

def troe(self,reaction,collider,label,epsilon=None,kTP='off'):
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
    Xdict=get_Xdict(self,reaction)
    # gas = ct.Solution("shortMech.yaml")
    # with open(self.foutName+".yaml") as f:
    #     test = yaml.safe_load(f)
    #     # print(test)
    # gas = ct.Solution(yaml=yaml.safe_dump(self.data))

    gas = ct.Solution(self.foutName+".yaml")
    logk_list=[]
    for i,P in enumerate(Xdict.keys()):
        for j,T in enumerate(Xdict[P]):
            gas.TPX=T,P*ct.one_atm,{collider:1.0}
            k_TP=gas.forward_rate_constants[self.rxnIdx]
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
    Xvec=get_Xvec(self,reaction)
    popt, pcov = curve_fit(f,Xvec,logk_list,p0=guess,maxfev=1000000,bounds=bounds)
    a0,n0,ea0,ai,ni,eai=popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]
    def numFmt(val):
        return round(float(val),3)
    if kTP=='on' and not epsilon:
        colDict = {
            'name': label,
            'low-P-rate-constant': {'A':numFmt(a0), 'b': numFmt(n0), 'Ea': numFmt(ea0)},
            'high-P-rate-constant': {'A': numFmt(ai), 'b': numFmt(ni), 'Ea': numFmt(eai)},
            'Troe': {'A': round(float(popt[6]),3), 'T3': 1.0e-30, 'T1': 1.0e+30}
        }
    elif kTP=='on' and epsilon:
        colDict = {
            'name': label,
            'efficiency': epsilon,
            'low-P-rate-constant': {'A':numFmt(a0), 'b': numFmt(n0), 'Ea': numFmt(ea0)},
            'high-P-rate-constant': {'A': numFmt(ai), 'b': numFmt(ni), 'Ea': numFmt(eai)},
            'Troe': {'A': round(float(popt[6]),3), 'T3': 1.0e-30, 'T1': 1.0e+30}
        }
    else:
        colDict = {
            'name': label,
            'efficiency': epsilon,
        }
    return colDict