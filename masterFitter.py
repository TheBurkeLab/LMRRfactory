"""
Class that allows for fitting of rate constants at various temperatures and pressures (k(T,P))
"""
# import sys, os
# import pkg_resources
# sys.path.append(pkg_resources.resource_filename('LMRRfactory', 'ext/cantera/build/python'))
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
import yaml
import copy
import re


class masterFitter:
    def __init__(self, T_ls, P_ls, colliderInput,mechInput,n_P=7, n_T=7, M_only=False):
        self.T_ls = T_ls
        self.P_ls = P_ls
        self.n_P=n_P
        self.n_T=n_T
        self.P_min = P_ls[0]
        self.P_max = P_ls[-1]
        self.T_min = T_ls[0]
        self.T_max = T_ls[-1]
        self.M_only=M_only

        # datapath = pkg_resources.resource_filename('LMRRfactory', 'data') + "/"

        # self.input = self.openyaml(inputFile)
        with open(colliderInput) as f:
            self.input = yaml.safe_load(f) # load input colliders
        with open(mechInput) as f:
            self.mech = yaml.safe_load(f) # load input mechanism
        self.mech = self.cleanedYAML(self.mech) # clean up 'NO' parsing errors
        with open("thirdbodydefaults.yaml") as f:
            self.defaults = yaml.safe_load(f) # load default colliders
        # Remove defaults colliders and reactions that were explictly provided by user
        self.defaults2 = self.deleteDuplicates()
        # Blend the user inputs and remaining collider defaults into a single YAML
        self.blend = self.blendedInput()
        # Sub the colliders into their corresponding reactions in the input mechanism
        self.outMech = self.zippedMech()
        



        self.pDepReactionNames=[]
        self.pDepReactions = []
        for reaction in self.mech['reactions']:
            if reaction.get('type') == 'falloff'  and 'Troe' in reaction:
                self.pDepReactions.append(reaction)
                self.pDepReactionNames.append(reaction['equation'])
            elif reaction.get('type') == 'pressure-dependent-Arrhenius':
                self.pDepReactions.append(reaction)
                self.pDepReactionNames.append(reaction['equation'])
            elif reaction.get('type') == 'Chebyshev':
                self.pDepReactions.append(reaction)
                self.pDepReactionNames.append(reaction['equation'])

        if len(self.pDepReactions)==0:
            print("No pressure-dependent reactions found in mechanism. Please choose another mechanism.")
        else:
            
            self.outMech = self.zippedMech()

    

    def cleanedYAML(self,data):
        cleanedData = data
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
        return cleanedData

    def deleteDuplicates(self):
        # defaults2 = {'reactions': []}
        # defaultRxnNames=[]
        # for defaultRxn in self.defaults['reactions']:
        #     defaultRxnNames.append(defaultRxn['equation'])
        # for inputRxn in self.input['reactions']:
        #     if inputRxn['equation'] in defaultRxnNames:
        #         defaultRxnNames.remove(inputRxn['equation'])
        # newReactionList = []
        # for defaultRxn in self.defaults['reactions']:
        #     if defaultRxn['equation'] in defaultRxnNames:
        #         newReactionList.append(defaultRxn)
        defaults2 = {'reactions': []}
        inputRxnNames=[]
        inputColliderNames=[]
        # bigInput = self.generalizedEquations(self.input)
        # for inputRxn in bigInput['reactions']:
        for inputRxn in self.input['reactions']:
            inputRxnNames.append(inputRxn['equation'])
            inputRxnColliderNames=[]
            for inputCol in inputRxn['colliders']:
                inputRxnColliderNames.append(inputCol['name'])
            inputColliderNames.append(inputRxnColliderNames)
        print(inputRxnNames)
        for defaultRxn in self.defaults['reactions']:
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
        return defaults2

    def blendedInput(self):
        # with open("defaults2.yaml", 'w') as outfile:
        #     yaml.dump(defaults2, outfile, default_flow_style=None,sort_keys=False)
        # with open("inputs2_double.yaml", 'w') as outfile:
        #     yaml.dump(self.generalizedEquations(self.input), outfile, default_flow_style=None,sort_keys=False)
        defaults2 = self.defaults2
        blend = {'reactions': []}
        speciesList = self.mech['phases'][0]['species']
        defaultRxnNames = []
        defaultColliderNames = []
        # for defaultRxn in self.generalizedEquations(defaults2)['reactions']:
        for defaultRxn in defaults2['reactions']:
            defaultRxnNames.append(defaultRxn['equation'])
            for defaultCol in defaultRxn['colliders']:
                defaultColliderNames.append(defaultCol['name'])
        # first fill it with all of the default reactions and colliders (which have valid species)
        # for defaultRxn in self.generalizedEquations(defaults2)['reactions']:
        for defaultRxn in defaults2['reactions']:
            flag = True
            for defaultCol in defaultRxn['colliders']:
                if defaultCol['name'] not in speciesList:
                    flag = False
            if flag == True:
                blend['reactions'].append(defaultRxn)

        blendRxnNames = []
        for blendRxn in blend['reactions']:
            blendRxnNames.append(blendRxn['equation'])
        
        # for inputRxn in self.generalizedEquations(self.input)['reactions']:
        for inputRxn in self.input['reactions']:
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
                print(reaction['equation'])
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
        with open('blend_double.yaml', 'w') as outfile:
                yaml.dump(blend, outfile, default_flow_style=None,sort_keys=False)
        return blend

    def zippedMech(self):
        blend=self.blend
        shortMechanism={
            'units': self.mech['units'],
            'phases': self.mech['phases'],
            'species': self.mech['species'],
            'reactions': []
            }
        blendRxnNames = []
        for rxn in blend['reactions']:
            blendRxnNames.append(rxn['equation'])
        for mech_rxn in self.mech['reactions']:
            if mech_rxn['equation'] in blendRxnNames:
                idx = blendRxnNames.index(mech_rxn['equation'])
                colliderM = blend['reactions'][idx]['colliders'][0]
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
        return shortMechanism

    def get_Xvec(self,reaction):
        Prange = self.P_ls
        Xvec=[]
        for i,P in enumerate(Prange):
            for j,T in enumerate(self.T_ls):
                Xvec.append([P,T])
        Xvec=np.array(Xvec)
        return Xvec.T

    def get_Xdict(self,reaction):
        Prange = self.P_ls
        Xdict={}
        for i,P in enumerate(Prange):
            Xdict[P]=self.T_ls
        return Xdict

    def get_YAML_kTP(self,reaction,collider):
        # gas = ct.Solution("shortMech.yaml")
        gas = ct.Solution(yaml=yaml.safe_dump(self.outMech))
        k_TP = []
        for T in self.T_ls:
            temp = []
            for P in self.P_ls:
                gas.TPX = T, P*ct.one_atm, {collider:1.0}
                val=gas.forward_rate_constants[gas.reaction_equations().index(reaction['equation'])]
                temp.append(val)
            k_TP.append(temp)
        return np.array(k_TP)
    
    def first_cheby_poly(self, x, n):
        '''Generate n-th order Chebyshev ploynominals of first kind.'''
        if n == 0: return 1
        elif n == 1: return x
        result = 2. * x * self.first_cheby_poly(x, 1) - self.first_cheby_poly(x, 0)
        m = 0
        while n - m > 2:
            result = 2. * x * result - self.first_cheby_poly(x, m+1)
            m += 1
        # print(result)
        return result
    def reduced_P(self,P):
        '''Calculate the reduced pressure.'''
        P_tilde = 2. * np.log10(P) - np.log10(self.P_min) - np.log10(self.P_max)
        P_tilde /= (np.log10(self.P_max) - np.log10(self.P_min))
        return P_tilde
    def reduced_T(self,T):
        '''Calculate the reduced temperature.'''
        T_tilde = 2. * T ** (-1) - self.T_min ** (-1) - self.T_max ** (-1)
        T_tilde /= (self.T_max ** (-1) - self.T_min ** (-1))
        return T_tilde
    def cheby_poly(self,reaction,collider):
        '''Fit the Chebyshev polynominals to rate constants.
            Input rate constants vector k should be arranged based on pressure.'''
        k_TP = self.get_YAML_kTP(reaction,collider)
        cheb_mat = np.zeros((len(k_TP.flatten()), self.n_T * self.n_P))
        for m, T in enumerate(self.T_ls):
            for n, P in enumerate(self.P_ls):
                for i in range(self.n_T):
                    for j in range(self.n_P):
                        T_tilde = self.reduced_T(T)
                        P_tilde = self.reduced_P(P)
                        T_cheb = self.first_cheby_poly(T_tilde, i)
                        P_cheb = self.first_cheby_poly(P_tilde, j)
                        cheb_mat[m * len(self.P_ls) + n, i * self.n_P + j] = T_cheb * P_cheb
        coef = np.linalg.lstsq(cheb_mat, np.log10(k_TP.flatten()),rcond=None)[0].reshape((self.n_T, self.n_P))
        return coef
    def get_cheb_table(self,reaction,collider,label,epsilon,kTP='off'):
        coef = self.cheby_poly(reaction,collider)
        if kTP=='on':
            colDict = {
                'name': label,
                'eps': epsilon,
                'temperature-range': [float(self.T_min), float(self.T_max)],
                'pressure-range': [f'{self.P_min:.3e} atm', f'{self.P_max:.3e} atm'],
                'data': []
            }
            for i in range(len(coef)):
                row=[]
                for j in range(len(coef[0])):
                    # row.append(f'{coef[i,j]:.4e}')
                    row.append(float(coef[i,j]))
                colDict['data'].append(row)
        else:
            colDict = {
                'name': label,
                'eps': epsilon,
            }
        
        return colDict
    
    def get_PLOG_table(self,reaction,collider,label,epsilon,kTP='off'):
        if kTP=='on':
            colDict = {
                'name': label,
                'eps': epsilon,
                'rate-constants': []
            }
            def arrhenius(T, A, n, Ea):
                return np.log(A) + n*np.log(T)+ (-Ea/(1.987*T))
            # gas = ct.Solution("shortMech.yaml")
            gas = ct.Solution(yaml=yaml.safe_dump(self.outMech))
            Xdict = self.get_Xdict(reaction)
            for i,P in enumerate(Xdict.keys()):
                k_list = []
                for j,T in enumerate(Xdict[P]):
                    gas.TPX = T, P*ct.one_atm, {collider:1}
                    k_T = gas.forward_rate_constants[gas.reaction_equations().index(reaction['equation'])]
                    k_list.append(k_T)
                k_list=np.array(k_list)
                popt, pcov = curve_fit(arrhenius, self.T_ls, np.log(k_list),maxfev = 30000)
                newDict = {
                    'P': f'{P:.3e} atm',
                    'A': float(popt[0]),
                    'b': float(popt[1]),
                    'Ea': float(popt[2]),
                }
                colDict['rate-constants'].append(newDict)
        else:
            colDict = {
                'name': label,
                'eps': epsilon,
            }
        
        return colDict

    def get_Troe_table(self,reaction,collider,label,epsilon,kTP='off'):
        def f(X,a0,n0,ea0,ai,ni,eai,Fcent):
            N= 0.75 - 1.27 * np.log10(Fcent) 
            c= -0.4 - 0.67 * np.log10(Fcent)
            d=0.14
            Rcal=1.987
            Rjoule=8.3145
            M = X[0]*ct.one_atm/Rjoule/X[1]/1000000.0
            k0 = a0 * (X[1] ** n0) * np.exp(-ea0 / (Rcal * X[1]))
            ki = ai * (X[1] ** ni) * np.exp(-eai / (Rcal * X[1]))
            logps = np.log10(k0) + np.log10(M) - np.log10(ki)
            den = logps + c
            den = den / (N - d * den)
            den = np.power(den, 2) + 1.0
            logF = np.log10(Fcent) / den
            logk_fit = np.log10(k0) + np.log10(M) + np.log10(ki) + logF - np.log10(ki + k0 * M)
            return logk_fit
        Xdict=self.get_Xdict(reaction)
        # gas = ct.Solution("shortMech.yaml")
        gas = ct.Solution(yaml=yaml.safe_dump(self.outMech))
        logk_list=[]
        for i,P in enumerate(Xdict.keys()):
            for j,T in enumerate(Xdict[P]):
                gas.TPX=T,P*ct.one_atm,{collider:1.0}
                k_TP=gas.forward_rate_constants[gas.reaction_equations().index(reaction['equation'])]
                logk_list.append(np.log10(k_TP))
        # NEED TO GENERALIZE THE FOLLOWING LINES
        # if "H + OH (+M)" in reaction:
        k0_g = [4.5300E+21, -1.8100E+00, 4.9870E+02]
        ki_g = [2.5100E+13, 0.234, -114.2]
        # # elif "H + O2 (+M)" in reaction:
        #     k0_g = [6.366e+20, -1.72, 524.8]
        #     ki_g = [4.7e+12,0.44,0.0]
        # # elif "H2O2 (+M)" in reaction:
        #     k0_g = [2.5e+24,-2.3, 4.8749e+04]
        #     ki_g = [2.0e+12,0.9,4.8749e+04]
        # # elif "NH2 (+M)" in reaction:
        #     k0_g = [1.6e+34,-5.49,1987.0]
        #     ki_g = [5.6e+14,-0.414,66.0]
        # # elif "NH3 <=>" in reaction:
        #     k0_g = [2.0e+16, 0.0, 9.315e+04]
        #     ki_g = [9.0e+16, -0.39, 1.103e+05]
        # # elif "HNO" in reaction:
        #     k0_g = [2.4e+14, 0.206, -1550.0]
        #     ki_g = [1.5e+15, -0.41, 0.0]
        guess = k0_g+ki_g+[1]
        bounds = (
                [1e-100, -np.inf, -np.inf, 1e-100, -np.inf, -np.inf, 1e-100],  # Lower bounds
                [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 1]         # Upper bounds
            )
        Xvec=self.get_Xvec(reaction)
        popt, pcov = curve_fit(f,Xvec,logk_list,p0=guess,maxfev=1000000,bounds=bounds)
        a0,n0,ea0,ai,ni,eai=popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]
        def numFmt(val):
            return round(float(val),3)
        if kTP=='on':
            colDict = {
                'name': label,
                'eps': epsilon,
                'low-P-rate-constant': {'A':numFmt(a0), 'b': numFmt(n0), 'Ea': numFmt(ea0)},
                'high-P-rate-constant': {'A': numFmt(ai), 'b': numFmt(ni), 'Ea': numFmt(eai)},
                'Troe': {'A': round(float(popt[6]),3), 'T3': 1.0e-30, 'T1': 1.0e+30}
            }
        else:
            colDict = {
                'name': label,
                'eps': epsilon,
            }
        return colDict

    def final_yaml(self,foutName,fit_fxn): # returns PLOG in LMRR YAML format
        newMechanism={
                'units': self.mech['units'],
                'phases': self.mech['phases'],
                'species': self.mech['species'],
                'reactions': []
                }
        for reaction in self.outMech['reactions']:
            if reaction.get('type')=='linear-Burke':
                colliderList=[]
                for i, col in enumerate(reaction['colliders']):
                    if i == 0:
                        colliderList.append(fit_fxn(reaction,reaction['reference-collider'],"M",col['eps'],kTP='on'))
                    elif len(list(reaction['colliders'][i].keys()))>3:
                        colliderList.append(fit_fxn(reaction,col['name'],col['name'],col['eps'],kTP='on'))
                    else:
                        colliderList.append(fit_fxn(reaction,col['name'],col['name'],col['eps'],kTP='off'))
                newMechanism['reactions'].append({
                    'equation': reaction['equation'],
                    'type': 'linear-Burke',
                    'reference-collider': reaction['reference-collider'],
                    'colliders': colliderList
                })
            else:
                newMechanism['reactions'].append(reaction)
        with open(foutName, 'w') as outfile:
            yaml.dump(newMechanism, outfile, default_flow_style=None,sort_keys=False)

    def colliders(self,foutName): # returns PLOG in LMRR YAML format
        with open(foutName, 'w') as outfile:
            yaml.dump(self.outMech, outfile, default_flow_style=None,sort_keys=False)
    def Troe(self,foutName): # returns PLOG in LMRR YAML format
        self.final_yaml(foutName,self.get_Troe_table)
    def PLOG(self,foutName): # returns PLOG in LMRR YAML format
        self.final_yaml(foutName,self.get_PLOG_table)
    def cheb2D(self,foutName): # returns Chebyshev in LMRR YAML format
        self.final_yaml(foutName,self.get_cheb_table)

########################################################################################


models = [
    {'name': 'Alzueta', 'path': 'alzuetamechanism.yaml'},
    # {'name': 'Mei', 'path': 'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Mei-2019\\mei-2019.yaml'},
    # # {'name': 'Glarborg', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Glarborg-2018\\glarborg-2018.yaml"},
    # {'name': 'Zhang', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Zhang-2017\\zhang-2017.yaml"},
    # {'name': 'Otomo', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Otomo-2018\\otomo-2018.yaml"},
    # {'name': 'Stagni', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Stagni-2020\\stagni-2020.yaml"},
    # # {'name': 'Shrestha', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Shrestha-2021\\shrestha-2021.yaml"},
    # {'name': 'Han', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Han-2021\\han-2021.yaml"},
    ]
# colours = ["xkcd:grey","xkcd:purple", "xkcd:teal", "orange", "r", "b", "xkcd:lime green", "xkcd:magenta", "xkcd:navy blue"]

date='Oct15'

for i, model in enumerate (models):
    # # INPUTS
    T_list=np.linspace(200,2000,100)
    # P_list=np.logspace(-12,12,num=120)
    P_list=np.logspace(-1,2,num=10)
    mF = masterFitter(T_list,P_list,colliderInput="testinput.yaml",mechInput=model['path'],n_P=7,n_T=7,M_only=True)
    path=f'outputs\\{date}'
    os.makedirs(path,exist_ok=True)
    mF.colliders(path+f"\\{model['name']}_LMRR.yaml")
    # mF.Troe(path+"\\LMRtest_Troe_M")
    # mF.PLOG(path+"\\LMRtest_PLOG_M")
    # mF.cheb2D(path+"\\LMRtest_cheb_M")








#########################

    # def yamlLoader2(self,fname):
    #     with open(fname) as f:
    #         data = yaml.safe_load(f)
    #     newMolecList = []
    #     # Prevent 'NO' from being misinterpreted as bool in species list
    #     for molec in data['phases'][0]['species']:
    #         if str(molec).lower()=="false":
    #             newMolecList.append("NO")
    #         else:
    #             newMolecList.append(molec)
    #     # Updated the YAML with the corrected species list
    #     data['phases'][0]['species'] = newMolecList
    #     for species in data['species']:
    #         name = str(species['name']).lower()
    #         if name == "false":
    #             species['name']="NO"
    #     # Prevent 'NO' from being misinterpreted as bool in efficiencies list found in
    #     # Troe falloff reactions
    #     for reaction in data['reactions']:
    #         effs = reaction.get('efficiencies')
    #         if effs is not None:
    #             for key in list(effs.keys()):
    #                 keyStr = str(key).lower()
    #                 if keyStr == "false":
    #                     reaction['efficiencies']["NO"] = reaction['efficiencies'].pop(key)
                        
    #         # if reaction.get('type')=='pressure-dependent-Arrhenius' or reaction.get('type')=='Chebyshev' or (reaction.get('type')=='falloff' and reaction.get('Troe') is not None):
    #         #     eqn = reaction['equation']
    #         #     reactants, products = eqn.split("<=>")
    #         #     if "(+M)" in reactants:
    #         #         reactants = reactants.replace("(+M)","").strip()
    #         #         products = products.replace("(+M)","").strip()
    #         #         if "+" in reactants:
    #         #             spec1,spec2 = reactants.split("+")
    #         #             spec1=spec1.strip()
    #         #             spec2=spec2.strip()
    #         #             if spec1==spec2:
    #         #                 reactants = f"2 {spec1}"
    #         #                 # rxn = copy.deepcopy(reaction)
    #         #                 reaction['equation'] = f"{reactants} (+M) <=> {products} (+M)"

    #         #         if "+" in products:
    #         #             spec1,spec2 = products.split("+")
    #         #             spec1=spec1.strip()
    #         #             spec2=spec2.strip()
    #         #             if spec1==spec2:
    #         #                 products = f"2 {spec1}"
    #         #                 # rxn = copy.deepcopy(reaction)
    #         #                 reaction['equation'] = f"{reactants} (+M) <=> {products} (+M)"
    #         #     else:
    #         #         reactants = reactants.strip()
    #         #         products = products.strip()
    #         #         print(reactants)
    #         #         if "+" in reactants:
    #         #             spec1, spec2 = reactants.split("+")
    #         #             spec1=spec1.strip()
    #         #             spec2=spec2.strip()
    #         #             if spec1==spec2:
    #         #                 reactants = f"2 {spec1}"
    #         #                 # rxn = copy.deepcopy(reaction)
    #         #                 reaction['equation'] = f"{reactants} <=> {products}"
    #         #         elif "+" in products:
    #         #             spec1, spec2 = products.split("+")
    #         #             spec1=spec1.strip()
    #         #             spec2=spec2.strip()
    #         #             if spec1==spec2:
    #         #                 products = f"2 {spec1}"
    #         #                 # rxn = copy.deepcopy(reaction)
    #         #                 reaction['equation'] = f"{reactants} <=> {products}"
    #     with open("newAlzueta.yaml", 'w') as outfile:
    #         yaml.dump(copy.deepcopy(data), outfile, default_flow_style=None,sort_keys=False)
    #     return data






    # def generalizedEquations(self,mech):
    #     newMech={'reactions': []}
    #     for rxn in mech['reactions']:
    #         eqn = rxn['equation']
    #         reactants, products = eqn.split('<=>')
    #         reactants = reactants.strip()
    #         products = products.strip()
    #         newMech['reactions'].append(rxn) #append the fwd direction
    #         rxn = copy.deepcopy(rxn)
    #         rxn['equation']=f"{products} <=> {reactants}"
    #         newMech['reactions'].append(rxn) #append the rev direction
    #         rxn = copy.deepcopy(rxn)
    #         products = products.replace("(+M)","").strip()
    #         reactants = reactants.replace("(+M)","").strip()
    #         rxn['equation']=f"{reactants} <=> {products}"
    #         newMech['reactions'].append(rxn)
    #         rxn = copy.deepcopy(rxn)
    #         rxn['equation']=f"{products} <=> {reactants}"
    #         newMech['reactions'].append(rxn)
    #         if "+" in reactants:
    #             spec1, spec2 = reactants.split("+")
    #             spec1=spec1.strip()
    #             spec2=spec2.strip()
    #             if spec1==spec2:
    #                 reactants = f"2 {spec1}"
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} <=> {products}"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} <=> {reactants}"
    #                 newMech['reactions'].append(rxn)
    #             else:
    #                 reactants = f"{spec2} + {spec1}"
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} <=> {products}"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} <=> {reactants}"
    #                 newMech['reactions'].append(rxn)
                    
    #         elif "+" in products:
    #             spec1, spec2 = products.split("+")
    #             spec1=spec1.strip()
    #             spec2=spec2.strip()
    #             if spec1==spec2:
    #                 products = f"2 {spec1}"
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} <=> {products}"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} <=> {reactants}"
    #                 newMech['reactions'].append(rxn)
                    
    #             else:
    #                 products = f"{spec2} + {spec1}"
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{reactants} <=> {products}"
    #                 newMech['reactions'].append(rxn)
    #                 rxn = copy.deepcopy(rxn)
    #                 rxn['equation']=f"{products} <=> {reactants}"
    #                 newMech['reactions'].append(rxn)
                    
    #     #     elif re.search(r'2\b',reactants) == True:
    #     #         spec = reactants.replace("2","")
    #     #         spec = spec.strip()
    #     #         reactants = f"{spec} + {spec}"
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
    #     #         newMech['reactions'].append(rxn)
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
    #     #         newMech['reactions'].append(rxn)
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{reactants} <=> {products}"
    #     #         newMech['reactions'].append(rxn)
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{products} <=> {reactants}"
    #     #         newMech['reactions'].append(rxn)
            
    #     #     elif re.search(r'2\b',products) == True:
    #     #         spec = products.replace("2","")
    #     #         spec = spec.strip()
    #     #         products = f"{spec} + {spec}"
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
    #     #         newMech['reactions'].append(rxn)
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
    #     #         newMech['reactions'].append(rxn)
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{reactants} <=> {products}"
    #     #         newMech['reactions'].append(rxn)
    #     #         rxn = copy.deepcopy(rxn)
    #     #         rxn['equation']=f"{products} <=> {reactants}"
    #     #         newMech['reactions'].append(rxn)
    #     # # with open('blend_double.yaml', 'w') as outfile:
    #     # #         yaml.dump(newMech, outfile, default_flow_style=None,sort_keys=False)
    #     return newMech