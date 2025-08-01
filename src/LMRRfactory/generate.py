# datapath = pkg_resources.resource_filename('LMRRfactory', 'data') + "/"
import yaml
import numpy as np
from scipy.optimize import least_squares
import copy
from collections import Counter
import re
import os
# import pkg_resources
from importlib.resources import files
import cantera as ct
import warnings
import io

warnings.filterwarnings("ignore")

class makeYAML:
    def __init__(self,mechInput=None, colliderInput=None, lmrrInput=None,outputPath=".",allPdep=False):
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
        # self.allPdep is option to apply generic 3b-effs to all p-dep reactions in mechanism that haven't already been explicitly specified in either "thirdbodydefaults.yaml" or self.colliderInput
        self.allPdep = False
        os.makedirs(outputPath,exist_ok=True)
        path=outputPath+'/'

        if mechInput:
            if colliderInput:
                self.colliderInput = colliderInput
            self.mechInput = mechInput
            self.foutName = os.path.basename(self.mechInput).replace(".yaml","")
            self.foutName = path + self.foutName + "_LMRR"
            if allPdep:
                self.allPdep = True
                self.foutName = self.foutName + "_allP"
            self.data = self.generateYAML()
        if lmrrInput:
            try:
                with open(lmrrInput) as f:
                    self.data = yaml.safe_load(f)
                self.foutName = path + lmrrInput.replace(".yaml","")
            except FileNotFoundError:
                print(f"Error: The file '{lmrrInput}' was not found.")

    def generateYAML(self):
        # data_path = pkg_resources.resource_filename('LMRRfactory', '/')
        data_path = str(files("LMRRfactory"))
        mech_obj = ct.Solution(self.mechInput)
        self.lookForPdep(mech_obj) # Verify that 'mech' has >=1 relevant p-dep reaction
        data = {
            'mech_obj': mech_obj,
            'mech_pes': self.getPES(mech_obj),
            'defaults': self.loadYAML(data_path+'/'+"thirdbodydefaults.yaml"),
            'input': self.loadYAML(self.colliderInput) if self.colliderInput is not None else None,
            'allPdep': self.allPdep, # True or False
        }
        # yml = self.loadYAML(self.mechInput)
        # data['extras'] = [yml[]]
        # Remove defaults colliders and reactions that were explictly provided by user
        self.deleteDuplicates(data)
        # Blend the user inputs and remaining collider defaults into a single YAML
        self.blendedInput(data)
        # Sub the colliders into their corresponding reactions in the input mechanism
        self.zippedMech(data)
        # self.validate(data['output'])
        self.saveYAML(data['output'], self.foutName+".yaml")
        print(f"LMR-R mechanism successfully generated and stored at "
            f"{self.foutName}.yaml")
        return data['output']

    def lookForPdep(self, gas):
        if not any(
            reaction.reaction_type in ['falloff-Troe','pressure-dependent-Arrhenius', 'Chebyshev', 'three-body-linear-Burke']
            for reaction in gas.reactions()
        ):
            raise ValueError("No pressure-dependent reactions found in mechanism."
                            " Please choose another mechanism.")
    
    def getPES(self,gas): #must input an equation that has already been normalized
        pes=[]
        for reaction in gas.reactions():
            compositions = []
            reactant_species = list(reaction.reactants.keys())
            reactant_coeffs = list(reaction.reactants.values())
            for i, reactant in enumerate(reactant_species):
                spec = gas.species(reactant)
                c = spec.composition
                c_scaled = {k.upper(): v*reactant_coeffs[i] for k, v in c.items()}
                compositions.append(c_scaled)
            counters = [Counter(comp) for comp in compositions]
            pes.append(sum(counters, Counter()))
        return pes

    def deleteDuplicates(self, data): # delete duplicates from thirdBodyDefaults
        newData = {'generic-colliders': data['defaults']['generic-colliders'],
                'reactions': []}
        inputRxnNames = None
        if data.get('input') is not None:
            if data['input'].get('reactions') is not None:
                inputRxnNames = [rxn['pes'] for rxn in data['input']['reactions']]
                inputColliderNames = [[col['composition'] for col in rxn['colliders']]
                                    for rxn in data['input']['reactions']]
        for defaultRxn in data['defaults']['reactions']:
            if inputRxnNames is not None and defaultRxn['pes'] in inputRxnNames:
                idx = inputRxnNames.index(defaultRxn['pes'])
                inputColliders = inputColliderNames[idx]
                newColliderList = [col for col in defaultRxn['colliders']
                                if col['composition'] not in inputColliders]
                if len(newColliderList)>0:
                    newData['reactions'].append({
                        'name': defaultRxn['name'],
                        'pes': defaultRxn['pes'],
                        'reference-collider': defaultRxn['reference-collider'],
                        'colliders': newColliderList
                    })
            else: # reaction isn't in input, so keep the entire default rxn
                newData['reactions'].append(defaultRxn)
        data['defaults']=newData

    def blendedInput(self, data):
        blendData = {'reactions': []}
        # speciesList = data['mech']['phases'][0]['species']
        # speciesList = [sp.composition for sp in data['mech_obj'].species()]
        # speciesList = [sp.name for sp in data['mech_obj'].species()]
        speciesDict = {}
        for sp in data['mech_obj'].species():
            speciesDict[sp.name] = sp.composition
        # print(speciesDict)
        # speciesList = data['mech_obj'].species()
        # print(speciesList)

        # first fill it with all of the default reactions and colliders (which have valid species)
        for defaultRxn in data['defaults']['reactions']:
            newCollList = []
            for col in defaultRxn['colliders']:
                print(col)
                if col['composition'] in list(speciesDict.values()):
                    newCollList.append(col)
            defaultRxn['colliders'] = newCollList
            blendData['reactions'].append(defaultRxn)
        defaultRxnNames = [rxn['pes'] for rxn in blendData['reactions']]
        if data.get('input') is not None:
            if data['input'].get('reactions') is not None:
                for inputRxn in data['input']['reactions']:
                    # Check if input reaction also exists in defaults file, otherwise add the entire input reaction to the blend as-is
                    if inputRxn['pes'] in defaultRxnNames:
                        idx = defaultRxnNames.index(inputRxn['pes'])
                        blendRxn = blendData['reactions'][idx]
                        # If reference colliders match, append new colliders, otherwise override with the user inputs
                        if inputRxn['reference-collider'] == blendRxn['reference-collider']:
                            newColliders = [col for col in inputRxn['colliders']
                                            if col['composition'] in list(speciesDict.values())]
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
                                            if col['composition'] in list(speciesDict.values())]
                            blendRxn['colliders'] = newColliders
                            # blendRxn['colliders'] = inputRxn['colliders']
                    else:
                        if all(col['composition'] in list(speciesDict.values()) for col in inputRxn['colliders']):
                            blendData['reactions'].append(inputRxn)
        data['blend']=blendData
        # print(blendData)

    def arrheniusFit(self, col):
        newCol = copy.deepcopy(col)
        temps=np.array(newCol['temperatures'])
        eps = np.array(newCol['efficiency'])
        def arrhenius_rate(T, A, beta, Ea):
            R = 1.987 # cal/molK
            return A * T**beta * np.exp(-Ea / (R * T))
        def fit_function(params, T, ln_eps):
            A, beta, Ea = params
            return np.log(arrhenius_rate(T, A, beta, Ea)) - ln_eps
        initial_guess = [3, 0.5, 50.0]
        result = least_squares(fit_function, initial_guess, args=(temps, np.log(eps)))
        A_fit, beta_fit, Ea_fit = result.x
        newCol['efficiency'] = {'A': float(round(A_fit.item(),8)),'b': float(round(beta_fit.item(),8)),'Ea': float(round(Ea_fit.item(),8))}
        newCol.pop('temperatures', None)
        newCol.pop('composition')
        return dict(newCol)

    def colliders(self,data,mech_rxn,blend_rxn=None,generic=False):
        # speciesList = [sp.name for sp in data['mech_obj'].species()]
        speciesDict = {}
        for sp in data['mech_obj'].species():
            speciesDict[sp.name] = sp.composition
        divisor = 1
        colliders=[]
        colliderNames=[]
        is_M_N2 = False
        troe_efficiencies={}
        # print(mech_rxn.reaction_type)
        if mech_rxn.reaction_type == 'falloff-Troe':
            troe_efficiencies_raw= mech_rxn.input_data.get('efficiencies', {})
            for sp_name in list(troe_efficiencies_raw.keys()):
                #make the keys the compositions instead of species names
                troe_efficiencies[speciesDict[sp_name]] = troe_efficiencies_raw[sp_name]
                print(troe_efficiencies)

            # print(troe_efficiencies)
        elif mech_rxn.reaction_type == 'three-body-linear-Burke': #case where we've used the linear Burke format so that Troe params can be used alongside a PLOG
            for c, col in enumerate(mech_rxn['colliders']):
                if c>0 and col['efficiency']['b']==0 and col['efficiency']['Ea']==0:
                    troe_efficiencies[speciesDict[col['name']]]=col['efficiency']['A']
        for comp, val in troe_efficiencies.items():
            # Check if N2 is the reference collider instead of Ar
            if comp=={'Ar': 1} and val!=0 and val !=1:
                is_M_N2 = True
                divisor = 1/val #ratio of N2:Ar
            if comp=={'Ar': 1} and val==0 :
                is_M_N2 = True
                divisor = 1 #imperfect solution, doesn't scale colliders, i.e. 1/val, to avoid dividing by zero but still acknowledges that rxn is w.r.t. N2
            # Give warning if both Ar and N2 are non-unity colliders
            if is_M_N2 and comp=={'N': 2} and val!=0 and val !=1:
                print(f"Warning: {mech_rxn} has both Ar and N2 as non-unity colliders!")
        if is_M_N2:
            if blend_rxn:
                divisors=[]
                # Extract T-dependent values for N2 if blend_rxn is provided
                for col in blend_rxn['colliders']:
                    if col['composition']=={'N': 2}:
                        divisors.append(col['efficiency']) #T-dep divisor of length 2 or 3
                # Make reaction-specific colliders wrt N2 and append to collider list 
                for col in blend_rxn['colliders']:
                    #Convert N2:Ar database entry to Ar:N2
                    if col['composition']=={'N': 2}:
                        col['composition']=={'Ar': 1}
                        col['name']=next(k for k, v in speciesDict.items() if v == col['composition'])
                        col['efficiency']=np.divide(1,col['efficiency'])
                        colliders.append(self.arrheniusFit(col))
                        colliderNames.append(col['composition'])
                    elif col['composition'] in list(speciesDict.values()):
                        # print(col['efficiency'])
                        for i in range(len(divisors)):
                            try:
                                col['efficiency'] = np.divide(col['efficiency'], divisors[i])
                                break
                            except:
                                pass
                        # print(col['efficiency'])
                        colliders.append(self.arrheniusFit(col))
                        colliderNames.append(col['composition'])
            # Add troe efficiencies that haven't already been given a value
            for comp, val in troe_efficiencies.items():
                already_given = comp in colliderNames
                if not already_given and not speciesDict[name]=={'N': 2}: #ignores the redundant n2=1 entry
                    colliders.append({
                        'name': next(k for k, v in speciesDict.items() if v == comp),
                        'efficiency': {'A':val,'b':0,'Ea':0 },
                        'note': 'present work',
                    })
                    colliderNames.append(speciesDict[name])
            if generic:
                for col in data['defaults']['generic-colliders']:
                    already_given = col['composition'] in colliderNames
                    if col['composition'] in list(speciesDict.keys()) and not already_given and not col['composition']=={'N': 2}:
                        if col.get('temperatures') is not None:
                            col['efficiency'] = np.divide(col['efficiency'],divisor)
                            colliders.append(self.arrheniusFit(col))
                        else:
                            colliders.append({
                                'name': next(k for k, v in speciesDict.items() if v == col['composition']),
                                'efficiency': {'A': col['efficiency']/divisor,'b':0,'Ea':0},
                                'note': col['note']
                            })
        else:
            if blend_rxn:
                # Make reaction-specific colliders wrt Ar and append to collider list 
                for col in blend_rxn['colliders']:
                    if col['composition'] in list(speciesDict.values()):
                        colliders.append(self.arrheniusFit(col))
                        colliderNames.append(col['composition'])
            # Add troe efficiencies that haven't already been given a value
            for comp, val in troe_efficiencies.items():
                # already_given = any(col['name'] == name for col in colliders)
                already_given = comp in colliderNames
                if not already_given and not comp=={'Ar': 1}:
                    colliders.append({
                        'name': next(k for k, v in speciesDict.items() if v == comp),
                        'efficiency': {'A':val,'b':0,'Ea':0 },
                        'note': 'present work',
                    })
                    colliderNames.append(comp)
            if generic:
                for col in data['defaults']['generic-colliders']:
                    already_given = col['composition'] in colliderNames
                    if col['composition'] in list(speciesDict.values()) and not already_given and not col['composition']=={'Ar': 1}:
                        if col.get('temperatures') is not None:
                            colliders.append(self.arrheniusFit(col))
                        else:
                            colliders.append({
                                'name': next(k for k, v in speciesDict.items() if v == col['composition']),
                                'efficiency': {'A': col['efficiency'],'b':0,'Ea':0},
                                'note': col['note']
                            })
        return colliders

    def to_builtin(self, obj):
        if isinstance(obj, dict):
            return {self.to_builtin(k): self.to_builtin(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self.to_builtin(i) for i in obj]
        elif hasattr(obj, 'as_dict'):
            return self.to_builtin(obj.as_dict())
        elif hasattr(obj, '__dict__'):
            return self.to_builtin(vars(obj))
        elif hasattr(obj, 'tolist'):  # NumPy arrays or similar
            return obj.tolist()
        else:
            return obj
        
    # def validate(self,soln):
    #     for rxn in soln.reactions():
    #         try:
    #             buf = io.StringIO()
    #             test_soln = {
    #                 'species': soln.species(),  # list of ct.Species objects
    #                 'thermo': soln.thermo_model,
    #                 'transport': soln.transport_model,
    #                 'reactions': [rxn],
    #             }
    #             testSoln = ct.Solution(**test_soln)
    #             testSoln.write_yaml(buf)
    #         except Exception as e:
    #             print(f"❌ {rxn.equation} could not be scanned, check the input file for errors.")
    #             print("   Error message:", e)
    #             print(rxn.input_data)


    def zippedMech(self, data):
        # input_data = gas.input_data
        newReactions = []
        blendRxnNames = [rxn['pes'] for rxn in data['blend']['reactions']]
        # for mech_rxn in data['mech']['reactions']:
        for i, mech_rxn in enumerate(data['mech_obj'].reactions()):
            # print(mech_rxn)
            pDep = False
            # Create the M-collider entry for the pressure-dependent reactions
            if mech_rxn.reaction_type in ['falloff-Troe','pressure-dependent-Arrhenius','Chebyshev','three-body-linear-Burke']:
                pDep = True
                if mech_rxn.reaction_type == 'three-body-linear-Burke':
                    d = self.to_builtin(mech_rxn.input_data['colliders'][0]) #use the pdep format given for collider M when rebuilding the reaction
                    d.pop("name")
                else:
                    d = self.to_builtin(mech_rxn.input_data)
                    d.pop("equation")
                    d.pop("efficiencies",None) #only applies to Troe reactions
                d.pop("duplicate", None)
                d.pop("units", None)
                if d.get('Troe') is not None:
                    d['Troe']=dict(d['Troe'])
                if d.get('low-P-rate-constant') is not None:
                    d['low-P-rate-constant']=dict(d['low-P-rate-constant'])
                if d.get('high-P-rate-constant') is not None:
                    d['high-P-rate-constant']=dict(d['high-P-rate-constant'])
                colliderM = {'name': 'M'}
                
                colliderM.update(dict(d))
            if pDep and data['mech_pes'][i] in blendRxnNames:
            # rxn is specifically covered either in defaults or user input
                newRxn = {
                    'equation': mech_rxn.equation,
                    # 'pes': data['mech_pes'][i],
                    'type': 'linear-Burke'
                }
                d = self.to_builtin(mech_rxn.input_data)
                if d.get('duplicate') is not None:
                    newRxn['duplicate'] = True
                if d.get('units') is not None:
                    newRxn['units'] = d['units']
                if 'note' in d and re.fullmatch(r'\n+', d['note']):
                    newRxn['note'] = ''
                idx = blendRxnNames.index(data['mech_pes'][i])
                blend_rxn = data['blend']['reactions'][idx]
                colliders = self.colliders(data,mech_rxn,blend_rxn=blend_rxn)
                newRxn['colliders'] = [colliderM] + colliders
                yaml_str = yaml.dump(newRxn, sort_keys=False)
                newRxn_obj = ct.Reaction.from_yaml(yaml_str,data['mech_obj'])
                newReactions.append(newRxn_obj)
                print(f"{mech_rxn} {dict(data['mech_pes'][i])} converted to LMR-R with ab initio parameters")
            elif pDep and data['allPdep']:
                # user has opted to have generic 3b effs applied to all p-dep reactions which lack a specification in thirdbodydefaults and testinput
                newRxn = {
                    'equation': mech_rxn.equation,
                    'type': 'linear-Burke'
                }
                d = self.to_builtin(mech_rxn.input_data)
                if d.get('duplicate') is not None:
                    newRxn['duplicate'] = True
                if d.get('units') is not None:
                    newRxn['units'] = d['units']
                if 'note' in d and re.fullmatch(r'\n+', d['note']):
                    newRxn['note'] = ''
                colliders = self.colliders(data,mech_rxn,generic=True)
                newRxn['colliders'] = [colliderM] + colliders
                yaml_str = yaml.dump(newRxn, sort_keys=False)
                newRxn_obj = ct.Reaction.from_yaml(yaml_str,data['mech_obj'])
                newReactions.append(newRxn_obj)
                print(f"{mech_rxn} {dict(data['mech_pes'][i])} converted to LMR-R with generic parameters")
            else: # just append it as-is
                d = mech_rxn.input_data
                if 'note' in d and re.fullmatch(r'\n+', d['note']):
                    mech_rxn.update_user_data({'note': ''})
                newReactions.append(mech_rxn)
        # output_data = {
        #     'species': data['mech_obj'].species(),  # list of ct.Species objects
        #     'thermo': data['mech_obj'].thermo_model,
        #     'transport': data['mech_obj'].transport_model,
        #     'reactions': newReactions,
        #     'name': 'outputMech'
        # }
        # print(data['mech_obj'].composite)
        # print(data['mech_obj'].kinetics_model)
        output_data = {
            'thermo': data['mech_obj'].thermo_model,
            'kinetics': data['mech_obj'].kinetics_model,
            'transport': data['mech_obj'].transport_model,
            'species': data['mech_obj'].species(),
            'reactions': newReactions,
            'name': 'outputMech',
        }
        data['output'] = ct.Solution(**output_data)

    def loadYAML(self, fName):
        with open(fName) as f:
            return yaml.safe_load(f)
    
    
    def saveYAML(self, dataSet, fName):
        dataSet.write_yaml(filename=fName,
                units={'length': 'cm', 'time': 's', 'quantity': 'mol', 'activation-energy': 'cal/mol'})