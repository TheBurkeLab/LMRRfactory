
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'src')))
from src.LMRRfactory import LMRRfactory

models = {
    'Alzueta': 'test\\data\\alzuetamechanism.yaml',
    'Mei': 'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Mei-2019\\mei-2019.yaml',
    # 'Glarborg': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Glarborg-2018\\glarborg-2018.yaml",
    'Zhang': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Zhang-2017\\zhang-2017.yaml",
    'Otomo': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Otomo-2018\\otomo-2018.yaml",
    'Stagni': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Stagni-2020\\stagni-2020.yaml",
    # 'Shrestha': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Shrestha-2021\\shrestha-2021.yaml",
    'Han': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Han-2021\\han-2021.yaml"
    }

for m in models.keys():
    base = {'mechanism': models[m]}
    base['colliders'] = 'test\\testinput.yaml'
    LMRRfactory(baseInput=base,allPdep=False,date='Oct22')

# import numpy as np
# T_list=np.linspace(200,2000,100)
# P_list=np.logspace(-1,2,num=10)
# mF.convertToTroe(T_list,P_list)