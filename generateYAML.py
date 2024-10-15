# datapath = pkg_resources.resource_filename('LMRRfactory', 'data') + "/"
import yaml
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
import copy

def main(self):
    # input = self.openyaml(inputFile)
    data = {}
    with open(self.colliderInput) as f:
        data['input'] = yaml.safe_load(f) # load input colliders
    with open(self.mechInput) as f:
        data['mech'] = yaml.safe_load(f) # load input mechanism
    cleanMechInput(data) # clean up 'NO' parsing errors in data['mech']
    saveYAML(data['mech'], "Alzueta_cleaned.yaml")
    lookForPdep(data) # Verify that mech has >=1 relevant p-dep reaction

    with open("thirdbodydefaults.yaml") as f:
        data['defaults'] = yaml.safe_load(f) # load default colliders
    # Remove defaults colliders and reactions that were explictly provided by user
    deleteDuplicates(data)
    saveYAML(data['defaults'], "defaults_uniqueOnly.yaml")
    # Blend the user inputs and remaining collider defaults into a single YAML
    blendedInput(data)
    saveYAML(data['blend'], "blendedInput.yaml")
    # Sub the colliders into their corresponding reactions in the input mechanism
    zippedMech(data)
    saveYAML(data['output'], "Alzueta_LMRR.yaml")
    saveYAML(data['output'], "Alzueta_LMRR.yaml")

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
    # # Uncomment if you want to see what the cleaned YAML looks like
    # with open("newAlzueta.yaml", 'w') as outfile:
    #     yaml.dump(copy.deepcopy(cleanedData), outfile, default_flow_style=None,sort_keys=False)
    data['mech']=cleanedData

def lookForPdep(data):
    # Raise an error if the input mech has no Troe, PLOG, or Chebyshev reactions
    if not any(
        reaction.get('type') in ['pressure-dependent-Arrhenius', 'Chebyshev'] or
        (reaction.get('type') == 'falloff' and 'Troe' in reaction)
        for reaction in data['mech']['reactions']
    ):
        raise ValueError("No pressure-dependent reactions found in mechanism. Please choose another mechanism.")

def deleteDuplicates(data): # delete duplicates from thirdBodyDefaults
    # defaults2 = {'reactions': []}
    # defaultRxnNames=[]
    # for defaultRxn in defaults['reactions']:
    #     defaultRxnNames.append(defaultRxn['equation'])
    # for inputRxn in input['reactions']:
    #     if inputRxn['equation'] in defaultRxnNames:
    #         defaultRxnNames.remove(inputRxn['equation'])
    # newReactionList = []
    # for defaultRxn in defaults['reactions']:
    #     if defaultRxn['equation'] in defaultRxnNames:
    #         newReactionList.append(defaultRxn)
    defaults2 = {'reactions': []}
    inputRxnNames=[]
    inputColliderNames=[]
    # bigInput = self.generalizedEquations(input)
    # for inputRxn in bigInput['reactions']:
    for inputRxn in input['reactions']:
        inputRxnNames.append(inputRxn['equation'])
        inputRxnColliderNames=[]
        for inputCol in inputRxn['colliders']:
            inputRxnColliderNames.append(inputCol['name'])
        inputColliderNames.append(inputRxnColliderNames)
    print(inputRxnNames)
    for defaultRxn in data['defaults']['reactions']:
        if defaultRxn['equation'] in inputRxnNames:
            idx = inputRxnNames.index(defaultRxn['equation'])
            defaultColliderNames=[]
            for defaultCol in defaultRxn['colliders']:
                defaultColliderNames.append(defaultCol['name'])
            # print(defaultRxn['equation'])
            for defaultCol in defaultRxn['colliders']:
                if defaultCol['name'] in inputColliderNames[idx]:
                    defaultColliderNames.remove(defaultCol['name'])
            # print(inputColliderNames[idx])
            # print(defaultColliderNames)
            newColliderList=[] #only contains colliders that aren't already in the input
            for defaultCol in defaultRxn['colliders']:
                if defaultCol['name'] in defaultColliderNames:
                    newColliderList.append(defaultCol)
            if len(newColliderList)>0:
                defaults2['reactions'].append({
                    'equation': defaultRxn['equation'],
                    'reference-collider': defaultRxn['reference-collider'],
                    'colliders': newColliderList
                })
        else: # reaction isn't in input, so keep the entire default rxn
            defaults2['reactions'].append(defaultRxn)
    data['defaults']=defaults2

def blendedInput(data):
    # with open("defaults2.yaml", 'w') as outfile:
    #     yaml.dump(defaults2, outfile, default_flow_style=None,sort_keys=False)
    # with open("inputs2_double.yaml", 'w') as outfile:
    #     yaml.dump(self.generalizedEquations(input), outfile, default_flow_style=None,sort_keys=False)
    blend = {'reactions': []}
    speciesList = data['mech']['phases'][0]['species']
    defaultRxnNames = []
    defaultColliderNames = []
    # for defaultRxn in self.generalizedEquations(defaults2)['reactions']:
    for defaultRxn in data['defaults']['reactions']:
        defaultRxnNames.append(defaultRxn['equation'])
        for defaultCol in defaultRxn['colliders']:
            defaultColliderNames.append(defaultCol['name'])
    # first fill it with all of the default reactions and colliders (which have valid species)
    # for defaultRxn in self.generalizedEquations(defaults2)['reactions']:
    for defaultRxn in data['defaults']['reactions']:
        flag = True
        for defaultCol in defaultRxn['colliders']:
            if defaultCol['name'] not in speciesList:
                flag = False
        if flag == True:
            blend['reactions'].append(defaultRxn)

    blendRxnNames = []
    for blendRxn in blend['reactions']:
        blendRxnNames.append(blendRxn['equation'])

    # for inputRxn in self.generalizedEquations(input)['reactions']:
    for inputRxn in data['input']['reactions']:
        if inputRxn['equation'] in blendRxnNames: #input reaction also exists in defaults file
            idx = blendRxnNames.index(inputRxn['equation'])
            # print(inputRxn['reference-collider'])
            if inputRxn['reference-collider'] == blend['reactions'][idx]['reference-collider']: #no blending conflicts bc colliders have same ref
                for inputCol in inputRxn['colliders']:
                    if inputCol['name'] in speciesList:
                        blend['reactions'][idx]['colliders'].append(inputCol)
            else: #blending conflict -> delete all default colliders and override with the user inputs
                print(f"The user-provided reference collider for {inputRxn['equation']}, ({inputRxn['reference-collider']}) does not match the program default ({blend['reactions'][idx]['reference-collider']}).")
                print(f"The default colliders have thus been deleted and the reaction has been completely overrided by (rather than blended with) the user's custom input values.")
                blend['reactions'][idx]['colliders'] = inputRxn['colliders']
        else:
            flag = True
            for inputCol in inputRxn['colliders']:
                if inputCol['name'] not in speciesList:
                    flag = False
            if flag == True:
                blend['reactions'].append(inputRxn)
    for reaction in blend['reactions']:
        for col in reaction['colliders']:
            # print(reaction['equation'])
            temperatures=np.array(col['temperatures'])
            eps = np.array(col['eps'])
            # epsLow=effs['epsLow']['A']
            # epsHigh=effs['epsHigh']['A']
            # rate_constants=np.array([epsLow,epsHigh])
            def arrhenius_rate(T, A, beta, Ea):
                # R = 8.314  # Gas constant in J/(mol K)
                R = 1.987 # cal/molK
                return A * T**beta * np.exp(-Ea / (R * T))
            def fit_function(params, T, ln_rate_constants):
                A, beta, Ea = params
                return np.log(arrhenius_rate(T, A, beta, Ea)) - ln_rate_constants
            initial_guess = [3, 0.5, 50.0]
            result = least_squares(fit_function, initial_guess, args=(temperatures, np.log(eps)))
            A_fit, beta_fit, Ea_fit = result.x
            col['eps'] = {'A': round(float(A_fit),5),'b': round(float(beta_fit),5),'Ea': round(float(Ea_fit),5)}
            del col['temperatures']
    # with open('blend_double.yaml', 'w') as outfile:
    #         yaml.dump(blend, outfile, default_flow_style=None,sort_keys=False)
    data['blend']=blend

def zippedMech(data):
    blend=data['blend']
    shortMechanism={
        'units': data['mech']['units'],
        'phases': data['mech']['phases'],
        'species': data['mech']['species'],
        'reactions': []
        }
    blendRxnNames = []
    for rxn in blend['reactions']:
        blendRxnNames.append(rxn['equation'])
    for mech_rxn in data['mech']['reactions']:
        if mech_rxn['equation'] in blendRxnNames:
            idx = blendRxnNames.index(mech_rxn['equation'])
            colliderM = blend['reactions'][idx]['colliders'][0] #figure out what this is for
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

            # colliderList.append(blend['reactions'][idx]['colliders'])
            shortMechanism['reactions'].append({
                        'equation': mech_rxn['equation'],
                        'type': 'linear-Burke',
                        'reference-collider': blend['reactions'][idx]['reference-collider'],
                        'colliders': colliderMlist + blend['reactions'][idx]['colliders']
                        })
        else:
            shortMechanism['reactions'].append(mech_rxn)
    data['output']=shortMechanism

def saveYAML(dataSet, fOutName):
    with open(fOutName, 'w') as outfile:
        yaml.dump(copy.deepcopy(dataSet['output']), outfile,
                  default_flow_style=None,
                  sort_keys=False)