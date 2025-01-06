# datapath = pkg_resources.resource_filename('LMRRfactory', 'data') + "/"
import yaml
import numpy as np
from scipy.optimize import least_squares
import copy
from collections import Counter
import re
import os
import pkg_resources
import yaml
import os
import numpy as np

import warnings

warnings.filterwarnings("ignore")

class makeYAML:
    def __init__(self,mechInput=None, colliderInput=None, lmrrInput=None,outputPath=None,allPdep=False,allPLOG=False):
        self.T_ls = None
        self.P_ls = None
        self.n_P= None
        self.n_T= None
        self.P_min = None
        self.P_max = None
        self.T_min = None
        self.T_max = None
        self.rxnIdx = None
        self.colliderInput=None
        # self.allPdep is option to apply generic 3b-effs to all p-dep reactions in
        # mechanism that haven't already been explicitly specified in either
        # "thirdbodydefaults.yaml" or self.colliderInput
        self.allPdep = False
        self.allPLOG = False

        # path='USSCI/factory_mechanisms'
        # path = outputPath
        # if date!="":
        #     path+=f'/{date}'
        os.makedirs(outputPath,exist_ok=True)
        path=outputPath+'/'

        if mechInput:
            if colliderInput:
                self.colliderInput = colliderInput
            self.mechInput = mechInput
            self.foutName = os.path.basename(self.mechInput).replace(".yaml","")
            self.foutName = path + self.foutName + "_LMRR"
            # try:
            # create a YAML in the LMRR format
            if allPdep:
                self.allPdep = True
                self.foutName = self.foutName + "_allP"
            if allPLOG:
                self.allPLOG = True
                self.foutName = self.foutName + "_allPLOG"
            self.data = self.generateYAML()
            # except ValueError:
            #     print(f"An LMR-R mechanism could not be generated using the "
            #             "baseInput files.")
        if lmrrInput:
            try:
                with open(lmrrInput) as f:
                    self.data = yaml.safe_load(f)
                self.foutName = path + lmrrInput.replace(".yaml","")
            except FileNotFoundError:
                print(f"Error: The file '{lmrrInput}' was not found.")

    def generateYAML(self):
        data_path = pkg_resources.resource_filename('LMRRfactory', '/')
        data = {
            'mech': self.loadYAML(self.mechInput), # load input mechanism}
            'defaults': self.loadYAML(data_path+"thirdbodydefaults.yaml"), # load default colliders
            'allPdep': self.allPdep, # True or False
            'allPLOG': self.allPLOG, # True or False
        }
        if self.colliderInput is not None:
            data['input']=self.loadYAML(self.colliderInput)
        else:
            data['input']=None
        
        self.cleanMechInput(data) # clean up 'NO' parsing errors in 'mech'
        self.lookForPdep(data) # Verify that 'mech' has >=1 relevant p-dep reaction
        if data.get('input') is not None:
            if data['input'].get('reactions') is not None:
                for reaction in data['input']['reactions']:
                    reaction['equation'] = self.normalize(reaction['equation'])
        if data.get('defaults') is not None:
            if data['defaults'].get('reactions') is not None:
                for reaction in data['defaults']['reactions']:
                    reaction['equation'] = self.normalize(reaction['equation'])
        # Remove defaults colliders and reactions that were explictly provided by user
        self.deleteDuplicates(data)
        # Blend the user inputs and remaining collider defaults into a single YAML
        self.blendedInput(data)
        # Sub the colliders into their corresponding reactions in the input mechanism
        self.zippedMech(data)
        self.saveYAML(data['output'], self.foutName+".yaml")
        print(f"LMR-R mechanism successfully generated and stored at "
            f"{self.foutName}.yaml")
        return data['output']

    def cleanMechInput(self, data):
        # Prevent 'NO' from being misinterpreted as bool in species list
        data['mech']['phases'][0]['species'] = [
            "NO" if str(molec).lower() == "false" else molec
            for molec in data['mech']['phases'][0]['species']
        ]
        for species in data['mech']['species']:
            if str(species['name']).lower() == "false":
                species['name']="NO"
        # Prevent 'NO' from being misinterpreted as bool in efficiencies list found in
        # Troe falloff reactions
        for reaction in data['mech']['reactions']:
            effs = reaction.get('efficiencies')
            if effs:
                reaction['efficiencies'] = {
                    "NO" if str(key).lower() == "false" else key: effs[key]
                    for key in effs
                }

    def lookForPdep(self, data):
        # Raise an error if the input mech has no Troe, PLOG, or Chebyshev reactions
        if not any(
            reaction.get('type') in ['pressure-dependent-Arrhenius', 'Chebyshev'] or
            (reaction.get('type') == 'falloff' and 'Troe' in reaction)
            for reaction in data['mech']['reactions']
        ):
            raise ValueError("No pressure-dependent reactions found in mechanism."
                            " Please choose another mechanism.")

    def normalize(self, equation):
        # Split the equation into reactants and products
        reactants, products = equation.split('=')
        reactants = reactants.strip().replace('(+M)', '').replace(' ', '').replace('<','').replace('>','')
        products = products.strip().replace('(+M)', '').replace(' ', '').replace('<','').replace('>','')
        def normalize_side(side):
            # Split into species and their coefficients
            species_list = re.split(r'\s*\+\s*', side)
            species_counter = Counter()
            for species in species_list:
                # Handle cases with coefficients like '2H' and '2 H'
                match = re.match(r'(\d*)\s*([^\d\s]\w*)', species)
                if not match:
                    raise ValueError(f"Incorrect formula for {equation} in input YAML.")
                coeff, name = match.groups()
                coeff = int(coeff) if coeff else 1  # Default to 1 if no coefficient
                species_counter[name] += coeff
            normalized_species = []
            for name in sorted(species_counter.keys()):  # Sort species alphabetically
                count = species_counter[name]
                normalized_species.extend([name]*count)
            normalized_side = ' + '.join(normalized_species)
            return normalized_side
        norm_reactants = normalize_side(reactants)
        norm_products = normalize_side(products)
        # Make it so that equations inputted in reverse directions are still deemed the same
        if norm_reactants > norm_products:
            norm_equation = f"{norm_reactants} <=> {norm_products}"
        else:
            norm_equation = f"{norm_products} <=> {norm_reactants}"
        return norm_equation

    def deleteDuplicates(self, data): # delete duplicates from thirdBodyDefaults
        newData = {'generic-colliders': data['defaults']['generic-colliders'],
                'reactions': []}
        inputRxnNames = None
        if data.get('input') is not None:
            if data['input'].get('reactions') is not None:
                inputRxnNames = [rxn['equation'] for rxn in data['input']['reactions']]
                inputColliderNames = [[col['name'] for col in rxn['colliders']]
                                    for rxn in data['input']['reactions']]
        for defaultRxn in data['defaults']['reactions']:
            if inputRxnNames is not None and defaultRxn['equation'] in inputRxnNames:
                idx = inputRxnNames.index(defaultRxn['equation'])
                inputColliders = inputColliderNames[idx]
                newColliderList = [col for col in defaultRxn['colliders']
                                if col['name'] not in inputColliders]
                if len(newColliderList)>0:
                    newData['reactions'].append({
                        'equation': defaultRxn['equation'],
                        'reference-collider': defaultRxn['reference-collider'],
                        'colliders': newColliderList
                    })
            else: # reaction isn't in input, so keep the entire default rxn
                newData['reactions'].append(defaultRxn)
        data['defaults']=newData

    def blendedInput(self, data):
        blendData = {'reactions': []}
        speciesList = data['mech']['phases'][0]['species']
        # first fill it with all of the default reactions and colliders (which have valid
        # species)
        for defaultRxn in data['defaults']['reactions']:
            newCollList = []
            for col in defaultRxn['colliders']:
                if col["name"] in speciesList:
                    newCollList.append(col)
            defaultRxn['colliders'] = newCollList
            blendData['reactions'].append(defaultRxn)
        defaultRxnNames = [rxn['equation'] for rxn in blendData['reactions']]
        if data.get('input') is not None:
            if data['input'].get('reactions') is not None:
                for inputRxn in data['input']['reactions']:
                    # Check if input reaction also exists in defaults file, otherwise add the entire
                    # input reaction to the blend as-is
                    if inputRxn['equation'] in defaultRxnNames:
                        idx = defaultRxnNames.index(inputRxn['equation'])
                        blendRxn = blendData['reactions'][idx]
                        # If reference colliders match, append new colliders, otherwise override
                        # with the user inputs
                        if inputRxn['reference-collider'] == blendRxn['reference-collider']:
                            newColliders = [col for col in inputRxn['colliders']
                                            if col['name'] in speciesList]
                            blendRxn['colliders'].extend(newColliders)
                        else:
                            print(f"User-provided reference collider for {inputRxn['equation']}, "
                                f"({inputRxn['reference-collider']}) does not match the program "
                                f"default ({blendData['reactions'][idx]['reference-collider']})."
                                f"\nThe default colliders have thus been deleted and the reaction"
                                f" has been completely overrided by (rather than blended with) "
                                f"the user's custom input values.")
                            blendRxn['reference-collider'] = inputRxn['reference-collider']
                            newColliders = [col for col in inputRxn['colliders']
                                            if col['name'] in speciesList]
                            blendRxn['colliders'] = newColliders
                            # blendRxn['colliders'] = inputRxn['colliders']
                    else:
                        if all(col['name'] in speciesList for col in inputRxn['colliders']):
                            blendData['reactions'].append(inputRxn)
        # # Convert collision efficiencies to arrhenius format and save to blended YAML
        # for reaction in blendData['reactions']:
        #     reaction['colliders']=[self.arrheniusFit(col) for col in reaction['colliders']]
        data['blend']=blendData

    def arrheniusFit(self, col):
        newCol = copy.deepcopy(col)
        temps=np.array(newCol['temperatures'])
        eps = np.array(newCol['efficiency'])
        def arrhenius_rate(T, A, beta, Ea):
            # R = 8.314  # Gas constant in J/(mol K)
            R = 1.987 # cal/molK
            return A * T**beta * np.exp(-Ea / (R * T))
        def fit_function(params, T, ln_eps):
            A, beta, Ea = params
            return np.log(arrhenius_rate(T, A, beta, Ea)) - ln_eps
        initial_guess = [3, 0.5, 50.0]
        result = least_squares(fit_function, initial_guess, args=(temps, np.log(eps)))
        A_fit, beta_fit, Ea_fit = result.x
        newCol['efficiency'] = {'A': round(A_fit.item(),8),'b': round(beta_fit.item(),8),'Ea': round(Ea_fit.item(),8)}
        newCol.pop('temperatures', None)
        return newCol

    def colliders(self,data,mech_rxn,blend_rxn=None,generic=False):
        speciesList = data['mech']['phases'][0]['species']
        divisor = 1
        colliders=[]
        colliderNames=[]
        is_M_N2 = False
        troe_efficiencies = mech_rxn.get('efficiencies', {})
        for name, val in troe_efficiencies.items():
            # Check if N2 is the reference collider instead of Ar
            if name.lower() =='ar' and val!=0 and val !=1:
                is_M_N2 = True
                divisor = 1/val #ratio of N2:Ar
            # Give warning if both Ar and N2 are non-unity colliders
            if is_M_N2 and name.lower() == 'n2' and val!=0 and val !=1:
                print(f"Warning: {mech_rxn['equation']} has both Ar and N2 as non-unity colliders!")

        if is_M_N2:
            if blend_rxn:
                divisors=[]
                # Extract T-dependent values for N2 if blend_rxn is provided
                for col in blend_rxn['colliders']:
                    if col['name'].lower()=='n2':
                        # divisor_Tdep=np.ones(len(col['efficiency']))
                        divisors.append(col['efficiency']) #T-dep divisor of length 2 or 3
                # Make reaction-specific colliders wrt N2 and append to collider list 
                for col in blend_rxn['colliders']:
                    #Convert N2:Ar database entry to Ar:N2
                    if col['name'].lower() == 'n2': 
                        col['name']='AR'
                        col['efficiency']=np.divide(1,col['efficiency'])
                        colliders.append(self.arrheniusFit(col))
                        colliderNames.append(col['name'].lower())
                    elif col['name'] in speciesList:
                        print(col['efficiency'])
                        for i in range(len(divisors)):
                            try:
                                col['efficiency'] = np.divide(col['efficiency'], divisors[i])
                                break
                            except:
                                pass
                        print(col['efficiency'])
                        colliders.append(self.arrheniusFit(col))
                        colliderNames.append(col['name'].lower())
            # Add troe efficiencies that haven't already been given a value
            for name, val in mech_rxn.get('efficiencies', {}).items():
                # already_given = any(col['name'] == name for col in colliders)
                already_given = name.lower() in colliderNames
                if not already_given and not name.lower()=='n2': #ignores the redundant n2=1 entry
                    colliders.append({
                        'name': name,
                        'efficiency': {'A':val,'b':0,'Ea':0 },
                        'note': 'present work',
                    })
                    colliderNames.append(name.lower())
            if generic:
                for col in data['defaults']['generic-colliders']:
                    already_given = str(col['name']).lower() in colliderNames
                    if col['name'] in speciesList and not already_given and not col['name'].lower()=='n2':
                        if col.get('temperatures') is not None:
                            col['efficiency'] = np.divide(col['efficiency'],divisor)
                            colliders.append(self.arrheniusFit(col))
                        else:
                            colliders.append({
                                'name': col['name'],
                                'efficiency': {'A': col['efficiency']/divisor,'b':0,'Ea':0},
                                'note': col['note']
                            })
        
        else:
            if blend_rxn:
                # Make reaction-specific colliders wrt Ar and append to collider list 
                for col in blend_rxn['colliders']:
                    if col['name'] in speciesList:
                        colliders.append(self.arrheniusFit(col))
                        colliderNames.append(col['name'].lower())
            # Add troe efficiencies that haven't already been given a value
            for name, val in troe_efficiencies.items():
                # already_given = any(col['name'] == name for col in colliders)
                already_given = name.lower() in colliderNames
                if not already_given and not name.lower()=='ar':
                    colliders.append({
                        'name': name,
                        'efficiency': {'A':val,'b':0,'Ea':0 },
                        'note': 'present work',
                    })
                    colliderNames.append(name.lower())
            if generic:
                for col in data['defaults']['generic-colliders']:
                    already_given = str(col['name']).lower() in colliderNames
                    if col['name'] in speciesList and not already_given and not col['name'].lower()=='ar':
                        if col.get('temperatures') is not None:
                            colliders.append(self.arrheniusFit(col))
                        else:
                            colliders.append({
                                'name': col['name'],
                                'efficiency': {'A': col['efficiency'],'b':0,'Ea':0},
                                'note': col['note']
                            })
        return colliders

    def zippedMech(self, data):
        newData={
            'units': data['mech']['units'],
            'phases': data['mech']['phases'],
            'species': data['mech']['species'],
            'reactions': []
            }
        blendRxnNames = [rxn['equation'] for rxn in data['blend']['reactions']]
        
        for mech_rxn in data['mech']['reactions']:
            # print(mech_rxn['equation'])
            pDep = False
            PLOG = False
            # troeEfficiencies={}
            # Create the M-collider entry for the pressure-dependent reactions
            if mech_rxn.get('type') == 'falloff' and 'Troe' in mech_rxn:
                pDep = True
                colliderM = {
                    'name': 'M',
                    'type': 'falloff',
                    'low-P-rate-constant': mech_rxn['low-P-rate-constant'],
                    'high-P-rate-constant': mech_rxn['high-P-rate-constant'],
                    'Troe': mech_rxn['Troe'],
                }
                # if mech_rxn.get('efficiencies') is not None:
                #     for key in mech_rxn['efficiencies'].keys():
                #         troeEfficiencies[key]=mech_rxn['efficiencies'][key]
            elif mech_rxn.get('type') == 'pressure-dependent-Arrhenius':
                pDep = True
                PLOG = True
                colliderM = {
                    'name': 'M',
                    'type': 'pressure-dependent-Arrhenius',
                    'rate-constants': mech_rxn['rate-constants']
                }
            elif mech_rxn.get('type') == 'Chebyshev':
                pDep = True
                colliderM = {
                    'name': 'M',
                    'type': 'Chebyshev',
                    'temperature-range': mech_rxn['temperature-range'],
                    'pressure-range': mech_rxn['pressure-range'],
                    'data': mech_rxn['data'],
                }
            # print(troeEfficiencies)
            # print(mech_rxn['equation'])
            print(self.normalize(mech_rxn['equation']))
            if pDep and self.normalize(mech_rxn['equation']) in blendRxnNames:
            # rxn is specifically covered either in defaults or user input
                newRxn = {
                    'equation': mech_rxn['equation'],
                    'type': 'linear-Burke'
                }
                if mech_rxn.get('duplicate') is not None:
                    newRxn['duplicate'] = True
                if mech_rxn.get('units') is not None:
                    newRxn['units'] = mech_rxn['units']
                idx = blendRxnNames.index(self.normalize(mech_rxn['equation']))
                blend_rxn = data['blend']['reactions'][idx]
                colliders = self.colliders(data,mech_rxn,blend_rxn=blend_rxn)
                newRxn['colliders'] = [colliderM] + colliders
                newData['reactions'].append(newRxn)
            elif pDep and data['allPdep']:
                # user has opted to have generic 3b effs applied to all p-dep reactions
                # which lack a specification in thirdbodydefaults and testinput
                # print(data['defaults']['allPdep''])
                newRxn = {
                    'equation': mech_rxn['equation'],
                    'type': 'linear-Burke'
                }
                if mech_rxn.get('duplicate') is not None:
                    newRxn['duplicate'] = True
                colliders = self.colliders(data,mech_rxn,generic=True)
                newRxn['colliders'] = [colliderM] + colliders
                newData['reactions'].append(newRxn)
            elif PLOG and data['allPLOG']:
                # user has opted to have generic 3b effs applied to all PLOG reactions
                # which lack a specification in thirdbodydefaults and testinput
                # print(data['defaults']['allPdep''])
                newRxn = {
                    'equation': mech_rxn['equation'],
                    'type': 'linear-Burke'
                }
                if mech_rxn.get('duplicate') is not None:
                    newRxn['duplicate'] = True
                colliders = self.colliders(data,mech_rxn,generic=True)
                newRxn['colliders'] = [colliderM] + colliders
                newData['reactions'].append(newRxn)
            else: # just append it as-is
                newData['reactions'].append(mech_rxn)
        data['output']=newData

    def loadYAML(self, fName):
        with open(fName) as f:
            return yaml.safe_load(f)

    def saveYAML(self, dataSet, fName):
        with open(fName, 'w') as outfile:
            yaml.dump(copy.deepcopy(dataSet), outfile,
                    default_flow_style=None,
                    sort_keys=False)