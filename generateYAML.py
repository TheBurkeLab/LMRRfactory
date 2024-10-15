# datapath = pkg_resources.resource_filename('LMRRfactory', 'data') + "/"
import yaml
import numpy as np
from scipy.optimize import least_squares
import copy

def generateYAML(self):
    data = {
        'input': loadYAML(self.colliderInput), # load input colliders
        'mech': loadYAML(self.mechInput), # load input mechanism}
        'defaults': loadYAML("thirdbodydefaults.yaml") # load default colliders
    }
    fname2 = self.foutName.replace(".yaml","")
    cleanMechInput(data) # clean up 'NO' parsing errors in 'mech'
    saveYAML(data['mech'], fname2+"_cleaned.yaml")
    lookForPdep(data) # Verify that 'mech' has >=1 relevant p-dep reaction
    # Remove defaults colliders and reactions that were explictly provided by user
    deleteDuplicates(data)
    saveYAML(data['defaults'], fname2+"_uniqueDefaults.yaml")
    # Blend the user inputs and remaining collider defaults into a single YAML
    blendedInput(data)
    saveYAML(data['blend'], fname2+"_blended.yaml")
    # Sub the colliders into their corresponding reactions in the input mechanism
    zippedMech(data)
    saveYAML(data['output'], self.foutName)
    return data['output']

def cleanMechInput(data):
    cleanedData = data['mech']
    newMolecList = []
    # Prevent 'NO' from being misinterpreted as bool in species list
    for molec in cleanedData['phases'][0]['species']:
        if str(molec).lower()=="false":
            newMolecList.append("NO")
        else:
            newMolecList.append(molec)
    # Updated the YAML with the corrected species list
    cleanedData['phases'][0]['species'] = newMolecList
    for species in cleanedData['species']:
        name = str(species['name']).lower()
        if name == "false":
            species['name']="NO"
    # Prevent 'NO' from being misinterpreted as bool in efficiencies list found in
    # Troe falloff reactions
    for reaction in cleanedData['reactions']:
        effs = reaction.get('efficiencies')
        if effs is not None:
            for key in list(effs.keys()):
                keyStr = str(key).lower()
                if keyStr == "false":
                    reaction['efficiencies']["NO"] = reaction['efficiencies'].pop(key)
    data['mech']=cleanedData

def lookForPdep(data):
    # Raise an error if the input mech has no Troe, PLOG, or Chebyshev reactions
    if not any(
        reaction.get('type') in ['pressure-dependent-Arrhenius', 'Chebyshev'] or
        (reaction.get('type') == 'falloff' and 'Troe' in reaction)
        for reaction in data['mech']['reactions']
    ):
        raise ValueError("No pressure-dependent reactions found in mechanism."
                         " Please choose another mechanism.")

def deleteDuplicates(data): # delete duplicates from thirdBodyDefaults
    newData = {'reactions': []}
    inputRxnNames = [rxn['equation'] for rxn in data['input']['reactions']]
    inputColliderNames = [[col['name'] for col in rxn['colliders']]
                          for rxn in data['input']['reactions']]
    for defaultRxn in data['defaults']['reactions']:
        if defaultRxn['equation'] in inputRxnNames:
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
                blendData['reactions'][idx]['colliders'] = inputRxn['colliders']
        else:
            if all(col['name'] in speciesList for col in inputRxn['colliders']):
                blendData['reactions'].append(inputRxn)

    # Convert collision efficiencies to arrhenius format and save to blended YAML
    for reaction in blendData['reactions']:
        [arrheniusFit(col) for col in reaction['colliders']]
    data['blend']=blendData

def arrheniusFit(col):
    temps=np.array(col['temperatures'])
    eps = np.array(col['eps'])
    def arrhenius_rate(T, A, beta, Ea):
        # R = 8.314  # Gas constant in J/(mol K)
        R = 1.987 # cal/molK
        return A * T**beta * np.exp(-Ea / (R * T))
    def fit_function(params, T, ln_rate_constants):
        A, beta, Ea = params
        return np.log(arrhenius_rate(T, A, beta, Ea)) - ln_rate_constants
    initial_guess = [3, 0.5, 50.0]
    result = least_squares(fit_function, initial_guess, args=(temps, np.log(eps)))
    A_fit, beta_fit, Ea_fit = result.x
    # Update eps values and remove temperatures
    col['eps'] = {'A': round(float(A_fit),5),
                  'b': round(float(beta_fit),5),
                  'Ea': round(float(Ea_fit),5)}
    del col['temperatures']


def zippedMech(data):
    newData={
        'units': data['mech']['units'],
        'phases': data['mech']['phases'],
        'species': data['mech']['species'],
        'reactions': []
        }
    blendRxnNames = []
    for rxn in data['blend']['reactions']:
        blendRxnNames.append(rxn['equation'])
    for mech_rxn in data['mech']['reactions']:
        if mech_rxn['equation'] in blendRxnNames:
            idx = blendRxnNames.index(mech_rxn['equation'])
            colliderM = data['blend']['reactions'][idx]['colliders'][0] #figure out what this is for
            colliderMlist=[]
            if mech_rxn['type'] == 'falloff' and 'Troe' in mech_rxn:
                colliderMlist.append({
                    'name': 'M',
                    'eps': {'A': 1, 'b': 0, 'Ea': 0},
                    'low-P-rate-constant': mech_rxn['low-P-rate-constant'],
                    'high-P-rate-constant': mech_rxn['high-P-rate-constant'],
                    'Troe': mech_rxn['Troe'],
                })
            elif mech_rxn['type'] == 'pressure-dependent-Arrhenius':
                colliderMlist.append({
                    'name': 'M',
                    'eps': {'A': 1, 'b': 0, 'Ea': 0},
                    'rate-constants': mech_rxn['rate-constants'],
                })
            elif mech_rxn['type'] == 'Chebyshev':
                colliderMlist.append({
                    'name': 'M',
                    'eps': {'A': 1, 'b': 0, 'Ea': 0},
                    'temperature-range': mech_rxn['temperature-range'],
                    'pressure-range': mech_rxn['pressure-range'],
                    'data': mech_rxn['data'],
                })

            # colliderList.append(data['blend']['reactions'][idx]['colliders'])
            newData['reactions'].append({
                'equation': mech_rxn['equation'],
                'type': 'linear-Burke',
                'reference-collider': data['blend']['reactions'][idx]['reference-collider'],
                'colliders': colliderMlist + data['blend']['reactions'][idx]['colliders']
                })
        else:
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