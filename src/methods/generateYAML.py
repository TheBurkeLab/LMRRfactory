# datapath = pkg_resources.resource_filename('LMRRfactory', 'data') + "/"
import yaml
import numpy as np
from scipy.optimize import least_squares
import copy
from collections import Counter
import re
import os

def generateYAML(self):
    data = {
        'input': loadYAML(self.colliderInput), # load input colliders
        'mech': loadYAML(self.mechInput), # load input mechanism}
        'defaults': loadYAML("data\\thirdbodydefaults.yaml"), # load default colliders
        'generic': self.allPdep # True or False
    }
    # foutName = self.foutName.replace(".yaml","")
    cleanMechInput(data) # clean up 'NO' parsing errors in 'mech'
    # saveYAML(data['mech'], f"{foutName}_cleaned.yaml")
    lookForPdep(data) # Verify that 'mech' has >=1 relevant p-dep reaction

    for reaction in data['input']['reactions']:
        reaction['equation'] = normalize(reaction['equation'])
    for reaction in data['defaults']['reactions']:
        reaction['equation'] = normalize(reaction['equation'])

    # Remove defaults colliders and reactions that were explictly provided by user
    deleteDuplicates(data)
    # saveYAML(data['defaults'], f"{foutName}_uniqueDefaults.yaml")
    # Blend the user inputs and remaining collider defaults into a single YAML
    blendedInput(data)
    # saveYAML(data['blend'], f"{foutName}_blended.yaml")
    # Sub the colliders into their corresponding reactions in the input mechanism
    zippedMech(data)
    saveYAML(data['output'], self.foutName+".yaml")
    print(f"LMR-R mechanism successfully generated from {os.path.basename(self.mechInput)}")
    print(f"The new file is stored at {self.foutName+".yaml"}\n")
    return data['output']

def cleanMechInput(data):
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

def lookForPdep(data):
    # Raise an error if the input mech has no Troe, PLOG, or Chebyshev reactions
    if not any(
        reaction.get('type') in ['pressure-dependent-Arrhenius', 'Chebyshev'] or
        (reaction.get('type') == 'falloff' and 'Troe' in reaction)
        for reaction in data['mech']['reactions']
    ):
        raise ValueError("No pressure-dependent reactions found in mechanism."
                         " Please choose another mechanism.")

def normalize(equation):
    # Split the equation into reactants and products
    reactants, products = equation.split('<=>')
    reactants = reactants.strip().replace('(+M)', '').replace(' ', '')
    products = products.strip().replace('(+M)', '').replace(' ', '')
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
    # print(norm_reaction)
    return norm_equation

def deleteDuplicates(data): # delete duplicates from thirdBodyDefaults
    newData = {'generic-colliders': data['defaults']['generic-colliders'],
               'reactions': []}
    inputRxnNames = [rxn['equation'] for rxn in data['input']['reactions']]
    inputColliderNames = [[col['name'] for col in rxn['colliders']]
                          for rxn in data['input']['reactions']]
    for defaultRxn in data['defaults']['reactions']:
        if defaultRxn['equation'] in inputRxnNames:
            # print(defaultRxn['equation'])
            # print(normalize(defaultRxn['equation']))
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

def blendedInput(data):
    blendData = {'reactions': []}
    speciesList = data['mech']['phases'][0]['species']

    # first fill it with all of the default reactions and colliders (which have valid
    # species)
    for defaultRxn in data['defaults']['reactions']:
        if all(col['name'] in speciesList for col in defaultRxn['colliders']):
            blendData['reactions'].append(defaultRxn)

    defaultRxnNames = [rxn['equation'] for rxn in blendData['reactions']]

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

    # Convert collision efficiencies to arrhenius format and save to blended YAML
    for reaction in blendData['reactions']:
        reaction['colliders']=[arrheniusFit(col) for col in reaction['colliders']]


        # colliders = [arrheniusFit(col)
        #             for col in reaction['colliders']
        #             if col['name'] in speciesList]
    data['blend']=blendData

def arrheniusFit(col):
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
    # Update eps values and remove temperatures
    newCol['efficiency'] = {'A': round(float(A_fit),5),
                  'b': round(float(beta_fit),5),
                  'Ea': round(float(Ea_fit),5)}
    newCol.pop('temperatures', None)
    return newCol


def zippedMech(data):
    newData={
        'units': data['mech']['units'],
        'phases': data['mech']['phases'],
        'species': data['mech']['species'],
        'reactions': []
        }
    blendRxnNames = [rxn['equation'] for rxn in data['blend']['reactions']]
    for mech_rxn in data['mech']['reactions']:
        pDep = False
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
        elif mech_rxn.get('type') == 'pressure-dependent-Arrhenius':
            pDep = True
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

        addReaction = True
        for rxn in newData['reactions']:
            if mech_rxn['equation'] == rxn['equation']:
                addReaction = False # Avoids adding a duplicate LMRR rxn that already exists in mech
        if pDep and addReaction:
            # rxn is specifically covered either in defaults or user input
            if normalize(mech_rxn['equation']) in blendRxnNames:
                idx = blendRxnNames.index(normalize(mech_rxn['equation']))
                blend_rxn = data['blend']['reactions'][idx]
                refCol = data['blend']['reactions'][idx]['reference-collider']
                colliders = blend_rxn['colliders']
                # user has opted to have generic 3b effs applied to all p-dep reactions
                # which lack a specification in thirdbodydefaults and testinput
                # print(data['defaults']['generic'])
            elif data['generic'] == True:
                refCol = 'AR' #just assumed, not aiming for perfection
                speciesList = data['mech']['phases'][0]['species']
                colliders = [arrheniusFit(col)
                             for col in data['defaults']['generic-colliders']
                             if col['name'] in speciesList]
            newData['reactions'].append({
                'equation': mech_rxn['equation'],
                'type': 'linear-Burke',
                'reference-collider': refCol,
                'colliders': [colliderM] + colliders,
            })
        elif pDep and not addReaction:
            del mech_rxn['duplicate']
            newData['reactions'].append(mech_rxn)
        else: # not a pressure dependent reaction, so just append it as-is
            newData['reactions'].append(mech_rxn)
    data['output']=newData

def loadYAML(fName):
    with open(fName) as f:
        return yaml.safe_load(f)

def saveYAML(dataSet, fName):
    with open(fName, 'w') as outfile:
        yaml.dump(copy.deepcopy(dataSet), outfile,
                  default_flow_style=None,
                  sort_keys=False)