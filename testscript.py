# import tools.yamlParser as parse
# import sys, os
# sys.path.append("C:\\Users\\pjsin\\Documents\\ext")
import ext.pyyaml.yaml as yaml

def openyaml(fname):
    with open(fname) as f:
        return yaml.load(f,Loader=yaml.SafeLoader)

reactions = {}
collider_input = openyaml("test\\inputs\\testinput.yaml")
for rxn in collider_input['reactions']:
    newdict={}
    for collider in rxn['colliders']:
        val = str(list(collider.values())[0])
        val2=val.replace("'","")
        newdict[list(collider.keys())[0]]=val2
    reactions[rxn['equation']]=newdict


chemical_input = openyaml(collider_input['chemical-mechanism'])
if len(chemical_input['phases'])>1:
    print("Error: multiphase kinetics are not supported.")
newdict = chemical_input['phases'][0]['species']

print(chemical_input['phases'][0]['species'])

