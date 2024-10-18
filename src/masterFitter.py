"""
Class that allows for fitting of rate constants at various temperatures and pressures (k(T,P))
"""

# import pkg_resources
# sys.path.append(pkg_resources.resource_filename('LMRRfactory', 'ext/cantera/build/python'))

from methods.generateYAML import generateYAML
from methods.chebyshevFitter import chebyshev
from methods.troeFitter import troe
from methods.plogFitter import plog
import yaml
import os
import numpy as np

import warnings

warnings.filterwarnings("ignore")


class masterFitter:
    def __init__(self,baseInput=None,lmrrInput=None,allPdep=False,date=""):
        self.T_ls = None
        self.P_ls = None
        self.n_P= None
        self.n_T= None
        self.P_min = None
        self.P_max = None
        self.T_min = None
        self.T_max = None
        # self.allPdep is option to apply generic 3b-effs to all p-dep reactions in
        # mechanism that haven't already been explicitly specified in either
        # "thirdbodydefaults.yaml" or self.colliderInput
        self.allPdep = False

        path='test\\outputs'
        if date!="":
            path+=f'\\{date}'
        os.makedirs(path,exist_ok=True)
        path+='\\'

        if baseInput:
            try:
                self.colliderInput = baseInput['colliders']
                self.mechInput = baseInput['mechanism']
                self.foutName = os.path.basename(self.mechInput).replace(".yaml","")
                self.foutName = path + self.foutName + "_LMRR"
                try:
                    # create a YAML in the LMRR format
                    if allPdep:
                        self.allPdep = True
                        self.foutName = self.foutName + "_allP"
                    self.data = generateYAML(self)
                except ValueError:
                    print(f"An LMR-R mechanism could not be generated using the baseInput files.")
            except ValueError:
                print("Base input must take the form: "
                        "{'colliders': 'filename1', 'mechanism': 'filename2'}")
        if lmrrInput:
            try:
                with open(lmrrInput) as f:
                    self.data = yaml.safe_load(f)
                self.foutName = path + lmrrInput.replace(".yaml","")
            except FileNotFoundError:
                print(f"Error: The file '{lmrrInput}' was not found.")

    def convertToTroe(self,T_ls, P_ls): # returns Troe in LMRR YAML format
        self.T_ls = T_ls
        self.P_ls = P_ls
        foutName2 = self.foutName+"_Troe.yaml"
        self._fittedYAML(foutName2,troe)

    def convertToPLOG(self,T_ls, P_ls): # returns PLOG in LMRR YAML format
        try:
            self.T_ls = T_ls
            self.P_ls = P_ls
            foutName = self.foutName+"_PLOG"
            self._fittedYAML(foutName,plog)
        except ValueError:
            print(f"Error: no LMR-R mechanism detected. If one already exists, it can be imported using LMRRfactory.load() -- otherwise, a new one can be generated using LMRRfactory.generate()")

    def convertToChebyshev(self,T_ls, P_ls,n_P=7, n_T=7): # returns Chebyshev in LMRR YAML format
        try:
            self.T_ls = T_ls
            self.P_ls = P_ls
            self.n_P=n_P
            self.n_T=n_T
            self.P_min = P_ls[0]
            self.P_max = P_ls[-1]
            self.T_min = T_ls[0]
            self.T_max = T_ls[-1]
            foutName = self.foutName+"_Chebyshev"
            self._fittedYAML(foutName,chebyshev)
        except ValueError:
            print(f"Error: no LMR-R mechanism detected. If one already exists, it can be imported using LMRRfactory.load() -- otherwise, a new one can be generated using LMRRfactory.generate()")

    def _fittedYAML(self,foutName,fit_fxn): # KEEP
        newMechanism={
                'units': self.data['units'],
                'phases': self.data['phases'],
                'species': self.data['species'],
                'reactions': []
                }
        for reaction in self.data['reactions']:
            if reaction.get('type')=='linear-Burke':
                colliderList=[]
                for i, col in enumerate(reaction['colliders']):
                    # print(col)
                    if i == 0:
                        colliderList.append(fit_fxn(self,reaction,reaction['reference-collider'],"M",col['eps'],kTP='on'))
                    elif len(list(reaction['colliders'][i].keys()))>3:
                        colliderList.append(fit_fxn(self,reaction,col['name'],col['name'],col['eps'],kTP='off'))
                    else:
                        colliderList.append(fit_fxn(self,reaction,col['name'],col['name'],col['eps'],kTP='off'))
                newMechanism['reactions'].append({
                    'equation': reaction['equation'],
                    'type': 'linear-Burke',
                    'reference-collider': reaction['reference-collider'],
                    'colliders': colliderList
                })
            else:
                newMechanism['reactions'].append(reaction)
        with open(foutName, 'w') as outfile:
            yaml.dump(newMechanism, outfile, default_flow_style=None,sort_keys=False)

########################################################################################


models = {
    'Alzueta': 'test\\data\\alzuetamechanism.yaml',
    # 'Mei': 'G:\\Mon disque\\Columbia\\Burke Lab\\07Mechanisms\\Mei-2019\\mei-2019.yaml',
    # 'Glarborg': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Glarborg-2018\\glarborg-2018.yaml",
    # 'Zhang': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Zhang-2017\\zhang-2017.yaml",
    # 'Otomo': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Otomo-2018\\otomo-2018.yaml",
    # 'Stagni': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Stagni-2020\\stagni-2020.yaml",
    # 'Shrestha': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Shrestha-2021\\shrestha-2021.yaml",
    # 'Han': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Han-2021\\han-2021.yaml"
    }

T_list=np.linspace(200,2000,100)
P_list=np.logspace(-1,2,num=10)

for m in models.keys():
    base = {'mechanism': models[m]}
    base['colliders'] = 'test\\testinput.yaml'
    mF = masterFitter(baseInput=base,allPdep=True,date='Oct17')
    mF.Troe(T_list,P_list)