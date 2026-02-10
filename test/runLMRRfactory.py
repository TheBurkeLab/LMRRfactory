import LMRRfactory
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--date', type=str)
parser.add_argument('--allPdep', type=str, default='False')
args = parser.parse_args()

allPdep = 'True'
allPLOG = 'True'
date = 'Feb09'

models = {
    # 'testmech': 'test/data/test_mech.yaml',
    'AramcoMech30': 'test/data/aramco30.yaml',
    # 'Klippenstein-2018': 'test/data/klippenstein-CNF2018-original.yaml',
    # 'Alzueta': 'test/data/alzuetamechanism.yaml',
    # 'Gutierrez-2025': 'test/data/gutierrez-2025.yaml',
    # 'kane-h2nno': 'test/data/glarborg_H2NNO_B.yaml',
    # 'kane-h2nno_uni': 'test/data/glarborg_H2NNO_B_unidirectionalexample.yaml',
    # # 'meng-h2nno2': 'test/data/glarborg_H2NNO2_B.yaml'
    }

for m in models.keys():
    # # LMRRfactory.makeYAML(mechInput=models[m],
    # #                     outputPath=f"test/outputs/{date}")
    # # LMRRfactory.makeYAML(mechInput=models[m],
    # #                     outputPath=f"test/outputs/{date}",reaction='HCO + HCO <=> CO + CH2O')
    # # LMRRfactory.makeYAML(mechInput=models[m],
    # #                     outputPath=f"test/outputs/{date}",reaction='NH2 + H (+M) <=> NH3 (+M)')
    # # LMRRfactory.makeYAML(mechInput=models[m],
    # #                     outputPath=f"test/outputs/{date}",reaction='NH2 + H (+AR) <=> NH3 (+AR)')
    # # LMRRfactory.makeYAML(mechInput=models[m],
    # #                     outputPath=f"test/outputs/{date}",reaction='C$_2$H$_2$+H(+M)=C$_2$H$_3$(+M)')
    # # LMRRfactory.makeYAML(mechInput=models[m],
    # #                     outputPath=f"test/outputs/{date}",reaction='H + O2 (+M) <=> HO2 (+M)') 
    # if allPdep == 'True':
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                     outputPath=f"test/outputs/{date}",
                            # allPdep=True)
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                     outputPath=f"test/outputs/{date}",allPdep=True,reaction='HCO + HCO <=> CO + CH2O')
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                     outputPath=f"test/outputs/{date}",allPdep=True,reaction='NH2 + H (+M) <=> NH3 (+M)')
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                     outputPath=f"test/outputs/{date}",allPdep=True,reaction='NH2 + H (+AR) <=> NH3 (+AR)')
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                     outputPath=f"test/outputs/{date}",allPdep=True,reaction='HOCO (+M) <=> CO2 + H (+M)')
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                     outputPath=f"test/outputs/{date}",allPdep=True,reaction='C2H2+H(+M)=>C2H3(+M)')
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                     outputPath=f"test/outputs/{date}",allPdep=True,reaction='C$_2$H$_2$+H(+M)=C$_2$H$_3$(+M)')
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                     outputPath=f"test/outputs/{date}",allPdep=True,reaction='C2H2 + H (+M) = C2H3 ( + M)')
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                     outputPath=f"test/outputs/{date}",allPdep=True,reaction='C2H2+H(+  M)=C2H3(+M)')
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                     outputPath=f"test/outputs/{date}",allPdep=True,reaction='C2H2 +H( + M)=C2H3(+M)')
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                     outputPath=f"test/outputs/{date}",allPdep=True,reaction='H2O$$$+N2 =>tNHNOH')
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                     outputPath=f"test/outputs/{date}",allPdep=True,reaction='NH2NO <= tNHNOH')
        # LMRRfactory.makeYAML(mechInput=models[m],
        #                 outputPath=f"test/outputs/{date}",allPdep=True,reaction='H + O2 (+M) <=> HO2 (+M)') 
    LMRRfactory.makeYAML(mechInput=models[m],
                         colliderInput='test/testinput.yaml',
                        outputPath=f"test/outputs/{date}",allPdep=True)