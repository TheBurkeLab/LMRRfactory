# python "test/testinput.py"
from LMRRfactory import makeYAML

models = {
    'Alzueta': 'test/data/alzuetamechanism.yaml'
    }

for m in models.keys():
    makeYAML(mechInput=models[m],
             outputPath='USSCI/factory_mechanisms',
             allPdep=False,date='Oct24')