# import tools.yamlParser as parse
# import sys, os
# sys.path.append("C:\\Users\\pjsin\\Documents\\ext")
import ext.pyyaml.yaml as yaml
# from io import StringIO

#collider_input = "test\\inputs\\testinput.yaml"

class masterFitter:
    def __init__(self, T_ls, P_ls, reactions, collider_input,n_P=7, n_T=7, M_only=False):
        self.T_ls = T_ls
        self.P_ls = P_ls
        self.n_P=n_P
        self.n_T=n_T
        # self.collider = collider
        # self.reaction=reaction
        self.reactions=reactions
        self.collider_input = collider_input
        self.P_min = P_ls[0]
        self.P_max = P_ls[-1]
        self.T_min = T_ls[0]
        self.T_max = T_ls[-1]
        # self.colours = colours
        self.M_only=M_only
        # self.P_ls_Troe=np.logspace(-1,3,num=60)
        # self.P_ls_Troe=P_ls

    def short_mech(self): # creates a version of Alzueta that only has the rxns for which we have eps_i vals
        def openyaml(fname):
            with open(fname) as f:
                return yaml.safe_load(f)
        keyReactions = {}
        collider_input = openyaml(self.collider_input)
        for rxn in collider_input['reactions']:
            newdict={}
            for collider in rxn['colliders']:
                newdict[list(collider.keys())[0]]=collider[list(collider.keys())[0]]
                # val = str(list(collider.values())[0])
                # val2=val.replace("'","")
                # newdict[list(collider.keys())[0]]=val2
            keyReactions[rxn['equation']]=newdict


        chemical_input = openyaml(collider_input['chemical-mechanism'])
        if len(chemical_input['phases'])>1:
            print("Error: multiphase kinetics are not supported.")

        shortMechanism={
            'units': chemical_input['units'],
            'phases': chemical_input['phases'],
            'species': chemical_input['species'],
            'reactions': []
            }

        for i, rxn in enumerate(chemical_input['reactions']):
            if rxn['equation'] in keyReactions.keys():
                colliderList = []
                if rxn['type'] == 'pressure-dependent-Arrhenius': 
                    colliderList.append({
                        'collider': 'M',
                        'eps': {'A': 1, 'b': 0, 'Ea': 0},
                        'rate-constants': rxn['rate-constants'],
                    })
                elif rxn['type'] == 'falloff' and 'Troe' in rxn:
                    colliderList.append({
                        'collider': 'M',
                        'eps': {'A': 1, 'b': 0, 'Ea': 0},
                        'low-P-rate-constant': rxn['low-P-rate-constant'],
                        'high-P-rate-constant': rxn['high-P-rate-constant'],
                        'Troe': rxn['Troe'],
                    })
                for j, col in enumerate(keyReactions[rxn['equation']].keys()):
                    colliderList.append({
                        'collider': col,
                        'eps': keyReactions[rxn['equation']][col]
                    })
                shortMechanism['reactions'].append({
                    'equation': rxn['equation'],
                    'type': 'linear-burke',
                    'collider-list': colliderList
                })
        return shortMechanism
    


        with open('shortdata.yml', 'w') as outfile:
            yaml.dump(shortMechanism, outfile, default_flow_style=None,sort_keys=False)


    def 


    

    for i, rxn in enumerate(chemical_input['reactions']):
        if rxn['equation'] in keyReactions.keys():
            colliderList = []
            if rxn['type'] == 'pressure-dependent-Arrhenius': 
                colliderList.append({
                    'collider': 'M',
                    'eps': {'A': 1, 'b': 0, 'Ea': 0},
                    'rate-constants': rxn['rate-constants'],
                })
            elif rxn['type'] == 'falloff' and 'Troe' in rxn:
                colliderList.append({
                    'collider': 'M',
                    'eps': {'A': 1, 'b': 0, 'Ea': 0},
                    'low-P-rate-constant': rxn['low-P-rate-constant'],
                    'high-P-rate-constant': rxn['high-P-rate-constant'],
                    'Troe': rxn['Troe'],
                })
            for j, col in enumerate(keyReactions[rxn['equation']].keys()):
                colliderList.append({
                    'collider': col,
                    'eps': keyReactions[rxn['equation']][col]
                })
            newMechanism['reactions'].append({
                'equation': rxn['equation'],
                'type': 'linear-burke',
                'collider-list': colliderList
            })
        else:
            newMechanism['reactions'].append(rxn)

    with open('data.yml', 'w') as outfile:
        yaml.dump(newMechanism, outfile, default_flow_style=None,sort_keys=False)
