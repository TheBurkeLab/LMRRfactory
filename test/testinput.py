# python "test/testinput.py"
from LMRRfactory import makeYAML

models = {
    'Alzueta': 'test/data/alzuetamechanism.yaml'
    }

for m in models.keys():
    makeYAML(mechInput=models[m],
             outputPath='test/outputs/Dec16')
    makeYAML(mechInput=models[m],
             outputPath='test/outputs/Dec16',
             allPdep=True)
    makeYAML(mechInput=models[m],
             outputPath='test/outputs/Dec16',
             allPLOG=True)