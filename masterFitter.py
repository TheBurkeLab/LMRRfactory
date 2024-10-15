"""
Class that allows for fitting of rate constants at various temperatures and pressures (k(T,P))
"""
# import sys, os
# import pkg_resources
# sys.path.append(pkg_resources.resource_filename('LMRRfactory', 'ext/cantera/build/python'))
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np

from generateYAML import generateYAML
from chebyshevFitter import chebyshev
from troeFitter import troe
from plogFitter import plog
import yaml


class masterFitter:
    def __init__(self, T_ls, P_ls, colliderInput,mechInput,yamlName,n_P=7, n_T=7, M_only=False):
        self.M_only=M_only
        self.colliderInput = colliderInput
        self.mechInput = mechInput
        self.yamlName = yamlName
        self.data = generateYAML() # create a YAML in the LMRR format

    def Troe(self,T_ls, P_ls): # returns PLOG in LMRR YAML format
        self.T_ls = T_ls
        self.P_ls = P_ls
        foutName = self.yamlName+"_Troe"
        self.fittedYAML(foutName,troe)
    def PLOG(self,T_ls, P_ls): # returns PLOG in LMRR YAML format
        self.T_ls = T_ls
        self.P_ls = P_ls
        foutName = self.yamlName+"_PLOG"
        self.fittedYAML(foutName,plog)
    def Chebyshev(self,T_ls, P_ls,n_P=7, n_T=7): # returns Chebyshev in LMRR YAML format
        self.T_ls = T_ls
        self.P_ls = P_ls
        self.n_P=n_P
        self.n_T=n_T
        self.P_min = P_ls[0]
        self.P_max = P_ls[-1]
        self.T_min = T_ls[0]
        self.T_max = T_ls[-1]
        foutName = self.yamlName+"_Chebyshev"
        self.fittedYAML(foutName,chebyshev)

    def fittedYAML(self,foutName,fit_fxn): # KEEP
        newMechanism={
                'units': self.mech['units'],
                'phases': self.mech['phases'],
                'species': self.mech['species'],
                'reactions': []
                }
        for reaction in self.outMech['reactions']:
            if reaction.get('type')=='linear-Burke':
                colliderList=[]
                for i, col in enumerate(reaction['colliders']):
                    if i == 0:
                        colliderList.append(fit_fxn(reaction,reaction['reference-collider'],"M",col['eps'],kTP='on'))
                    elif len(list(reaction['colliders'][i].keys()))>3:
                        colliderList.append(fit_fxn(reaction,col['name'],col['name'],col['eps'],kTP='on'))
                    else:
                        colliderList.append(fit_fxn(reaction,col['name'],col['name'],col['eps'],kTP='off'))
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


models = [
    {'name': 'Alzueta', 'path': 'alzuetamechanism.yaml'},
    # {'name': 'Mei', 'path': 'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Mei-2019\\mei-2019.yaml'},
    # # {'name': 'Glarborg', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Glarborg-2018\\glarborg-2018.yaml"},
    # {'name': 'Zhang', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Zhang-2017\\zhang-2017.yaml"},
    # {'name': 'Otomo', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Otomo-2018\\otomo-2018.yaml"},
    # {'name': 'Stagni', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Stagni-2020\\stagni-2020.yaml"},
    # # {'name': 'Shrestha', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Shrestha-2021\\shrestha-2021.yaml"},
    # {'name': 'Han', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Han-2021\\han-2021.yaml"},
    ]
# colours = ["xkcd:grey","xkcd:purple", "xkcd:teal", "orange", "r", "b", "xkcd:lime green", "xkcd:magenta", "xkcd:navy blue"]

date='Oct15'

for i, model in enumerate (models):
    # # INPUTS
    T_list=np.linspace(200,2000,100)
    # P_list=np.logspace(-12,12,num=120)
    P_list=np.logspace(-1,2,num=10)
    mF = masterFitter(T_list,P_list,colliderInput="testinput.yaml",mechInput=model['path'],n_P=7,n_T=7,M_only=True)
    path=f'outputs\\{date}'
    os.makedirs(path,exist_ok=True)
    mF.colliders(path+f"\\{model['name']}_LMRR.yaml")
    mF.Troe("troeTest.yaml")
    # mF.Troe(path+"\\LMRtest_Troe_M")
    # mF.PLOG(path+"\\LMRtest_PLOG_M")
    # mF.cheb2D(path+"\\LMRtest_cheb_M")



    # def generalizedEquations(self,mech):
    #     newMech={'reactions': []}
    #     for rxn in mech['reactions']:
    #         eqn = rxn['equation']
    #         reactants, products = eqn.split('<=>')
    #         reactants = reactants.strip()
    #         products = products.strip()
    #         newMech['reactions'].append(rxn) #append the fwd direction
    #         rxn = copy.deepcopy(rxn)
    #         rxn['equation']=f"{products} <=> {reactants}"
    #         newMech['reactions'].append(rxn) #append the rev direction
    #         rxn = copy.deepcopy(rxn)
    #         products = products.replace("(+M)","").strip()
    #         reactants = reactants.replace("(+M)","").strip()
    #         rxn['equation']=f"{reactants} <=> {products}"
    #         newMech['reactions'].append(rxn)
    #         rxn = copy.deepcopy(rxn)
    #         rxn['equation']=f"{products} <=> {reactants}"
    #         newMech['reactions'].append(rxn)
    #         if "+" in reactants:
    #             spec1, spec2 = reactants.split("+")
    #             spec1=spec1.strip()
    #             spec2=spec2.strip()
    #             if spec1==spec2:
    #                 reactants = f"2 {spec1}"
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} <=> {products}"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} <=> {reactants}"
    #                 newMech['reactions'].append(rxn)
    #             else:
    #                 reactants = f"{spec2} + {spec1}"
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} <=> {products}"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} <=> {reactants}"
    #                 newMech['reactions'].append(rxn)
                    
    #         elif "+" in products:
    #             spec1, spec2 = products.split("+")
    #             spec1=spec1.strip()
    #             spec2=spec2.strip()
    #             if spec1==spec2:
    #                 products = f"2 {spec1}"
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} <=> {products}"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} <=> {reactants}"
    #                 newMech['reactions'].append(rxn)
                    
    #             else:
    #                 products = f"{spec2} + {spec1}"
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} <=> {products}"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} <=> {reactants}"
    #                 newMech['reactions'].append(rxn)
                    
    #     #     elif re.search(r'2\b',reactants) == True:
    #     #         spec = reactants.replace("2","")
    #     #         spec = spec.strip()
    #     #         reactants = f"{spec} + {spec}"
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
    #     #         newMech['reactions'].append(rxn)
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
    #     #         newMech['reactions'].append(rxn)
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{reactants} <=> {products}"
    #     #         newMech['reactions'].append(rxn)
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{products} <=> {reactants}"
    #     #         newMech['reactions'].append(rxn)
            
    #     #     elif re.search(r'2\b',products) == True:
    #     #         spec = products.replace("2","")
    #     #         spec = spec.strip()
    #     #         products = f"{spec} + {spec}"
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
    #     #         newMech['reactions'].append(rxn)
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
    #     #         newMech['reactions'].append(rxn)
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{reactants} <=> {products}"
    #     #         newMech['reactions'].append(rxn)
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{products} <=> {reactants}"
    #     #         newMech['reactions'].append(rxn)
    #     # # with open('blend_double.yaml', 'w') as outfile:
    #     # #         yaml.dump(newMech, outfile, default_flow_style=None,sort_keys=False)
    #     return newMech