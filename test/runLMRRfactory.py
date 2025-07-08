import LMRRfactory
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--date', type=str)
parser.add_argument('--allPdep', type=str, default='False')
args = parser.parse_args()

allPdep = 'True'
allPLOG = 'True'
date = 'May07'

models = {
    'AramcoMech30': 'test/data/aramco30.yaml',
    }

for m in models.keys():
    LMRRfactory.makeYAML(mechInput=models[m],
                        outputPath=f"test/outputs/{date}")
    if allPdep == 'True':
        LMRRfactory.makeYAML(mechInput=models[m],
                            outputPath=f"test/outputs/{date}",
                            allPdep=True)