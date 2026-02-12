import yaml
import numpy as np
import copy
from collections import Counter
import re
import os
import io
from importlib.resources import files
import cantera as ct
import warnings
import contextlib

warnings.filterwarnings("ignore")

class makeYAML:
    def __init__(self, mechInput, colliderInput=None, outputPath=".", allPdep=False, reaction=None):
        self.T_ls = None
        self.P_ls = None
        self.n_P = None
        self.n_T = None
        self.P_min = None
        self.P_max = None
        self.T_min = None
        self.T_max = None
        self.rxnIdx = None
        self.colliderInput = None
        self.units = self._loadYAML(mechInput).get("units", {})
        self.allPdep = allPdep
        self.input = self._loadYAML(colliderInput) if colliderInput else None
        self.mechInput = mechInput
        self.reaction = reaction
        os.makedirs(outputPath, exist_ok=True)
        self.foutName = f"{outputPath}/{os.path.basename(self.mechInput).replace('.yaml', '_LMRR')}"
        if allPdep:
            self.foutName = f"{self.foutName}_allP"
        self.mech_obj = ct.Solution(mechInput)
        self._lookForPdep() # Verify that 'mech' has >=1 relevant p-dep reaction
        self.mech_pes = self._getPES()
        self.defaults = self._loadYAML(f"{str(files("LMRRfactory"))}/thirdbodydatabase.yaml")
        self._normalizedKeys() # normalize species as uppercase
        self.species_dict = {}
        for sp in self.mech_obj.species():
            self.species_dict[sp.name.upper()] = dict(sp.composition.items())
        # Remove defaults colliders and reactions that were explictly provided by user
        self._deleteDuplicates()
        # Blend the user inputs and remaining collider defaults into a single YAML
        self.blend = self._blendedInput()
        # Sub the colliders into their corresponding reactions in the input mechanism
        self.skipSave = False
        self.output = self._zippedMech()
        # self.validate()
        if not self.skipSave:
            self._saveYAML()
            print(f"New mechanism generated and stored at "
                f"{self.foutName}.yaml")

    def _normalizedKeys(self):
        def capitalize(mapping):
            return {k.capitalize(): v for k, v in mapping.items()}
        for defaultRxn in self.defaults['reactions']:
            for col in defaultRxn['colliders']:
                col['composition'] = capitalize(col['composition'])
            defaultRxn['reference-collider'] = defaultRxn['reference-collider'].upper()
            defaultRxn['pes'] = capitalize(defaultRxn['pes'])
        if self.input and self.input.get('reactions'):
            for inputRxn in self.input['reactions']:
                for col in inputRxn['colliders']:
                    col['composition'] = capitalize(col['composition'])
                inputRxn['reference-collider'] = inputRxn['reference-collider'].upper()
                if inputRxn['reference-collider'] != 'AR':
                    raise ValueError(
                        f"Reference collider '{inputRxn['reference-collider']}' in collider input "
                        f"reaction '{inputRxn['name']}' is not supported. Only parameters scaled "
                        f"with respect to AR as the reference collider are allowed in the collider "
                        f"input file.")
                inputRxn['pes'] = capitalize(inputRxn['pes'])

    def _normalizedUserRxn(self):
        s = self.reaction.strip()
        arrow_match = re.search(r"(<=>|=>|<=|=)", s)
        arrow = arrow_match.group(1) if arrow_match else "="
        s = re.sub(r"\$_\{(\d+)\}\$", r"\1", s)
        s = re.sub(r"\$_(\d+)\$", r"\1", s)
        s = s.replace("$", "")
        s = re.sub(r"\(\s*\+\s*([A-Za-z0-9_]+)\s*\)", r"(+\1)", s)
        s = re.sub(r"(?<=\w)\(\+([A-Za-z0-9_]+)\)", r" (+\1)", s)
        s = re.sub(
            r"\(\+([A-Za-z0-9_]+)\)",
            lambda m: f"__TB_{m.group(1)}__",
            s
        )
        s = re.sub(r"\s*(<=>|=>|<=|=)\s*", f" {arrow} ", s)
        s = re.sub(r"\s*\+\s*", " + ", s)
        s = re.sub(r"__TB_([A-Za-z0-9_]+)__", r"(+\1)", s)
        s = re.sub(r"\s+", " ", s).strip()
        return s


    def _lookForPdep(self):
        if not any(
            reaction.reaction_type in ['falloff-Troe','pressure-dependent-Arrhenius', 'Chebyshev', 'three-body-linear-Burke']
            for reaction in self.mech_obj.reactions()
        ):
            raise ValueError("No pressure-dependent reactions found in mechanism."
                            " Please choose another mechanism.")
    
    def _getPES(self):
        pes = []
        for reaction in self.mech_obj.reactions():
            compositions = []
            reactant_species = list(reaction.reactants.keys())
            reactant_coeffs = list(reaction.reactants.values())
            for i, reactant in enumerate(reactant_species):
                spec = self.mech_obj.species(reactant)
                c = spec.composition
                c_scaled = {k.upper(): v*reactant_coeffs[i] for k, v in c.items()}
                compositions.append(c_scaled)
            counters = [Counter(comp) for comp in compositions]
            pes.append(sum(counters, Counter()))
        return pes

    def _deleteDuplicates(self):
        newData = {'generic-colliders': self.defaults['generic-colliders'],
                   'reactions': []}
        inputRxnNames = None
        if self.input:
            if self.input.get('reactions'):
                inputRxnNames = [rxn['pes'] for rxn in self.input['reactions']]
                inputColliderNames = [[col['composition'] for col in rxn['colliders']]
                                    for rxn in self.input['reactions']]
        for defaultRxn in self.defaults['reactions']:
            if inputRxnNames and defaultRxn['pes'] in inputRxnNames:
                idx = inputRxnNames.index(defaultRxn['pes'])
                inputColliders = inputColliderNames[idx]
                newColliderList = [col for col in defaultRxn['colliders']
                                if col['composition'] not in inputColliders]
                if newColliderList:
                    newData['reactions'].append({
                        'name': defaultRxn['name'],
                        'pes': defaultRxn['pes'],
                        'reference-collider': defaultRxn['reference-collider'],
                        'colliders': newColliderList
                    })
            else: # reaction isn't in input, so keep the entire default rxn
                newData['reactions'].append(defaultRxn)
        self.defaults = newData

    def _blendedInput(self):
        blendData = {'reactions': []}
        # first fill it with all of the default reactions and colliders (which have valid species)
        for defaultRxn in self.defaults['reactions']:
            newCollList = []
            for col in defaultRxn['colliders']:
                if col['composition'] in list(self.species_dict.values()):
                    newCollList.append(col)
            defaultRxn['colliders'] = newCollList
            blendData['reactions'].append(defaultRxn)
        defaultRxnNames = [rxn['pes'] for rxn in blendData['reactions']]
        if self.input:
            if self.input.get('reactions'):
                for inputRxn in self.input['reactions']:
                    # Check if input reaction also exists in defaults file, otherwise add the entire input reaction to the blend as-is
                    if inputRxn['pes'] in defaultRxnNames:
                        idx = defaultRxnNames.index(inputRxn['pes'])
                        blendRxn = blendData['reactions'][idx]
                        # If reference colliders match, append new colliders, otherwise override with the user inputs
                        if inputRxn['reference-collider'] == blendRxn['reference-collider']:
                            newColliders = [col for col in inputRxn['colliders']
                                            if col['composition'] in list(self.species_dict.values())]
                            blendRxn['colliders'].extend(newColliders)
                        else:
                            print(f"User-provided reference collider for {inputRxn['equation']}, "
                                f"({inputRxn['reference-collider']}) does not match the program "
                                f"default ({blendData['reactions'][idx]['reference-collider']})."
                                f"\nThe default colliders have thus been deleted and the reaction"
                                f" has been completely overrided by (rather than blended with) "
                                f"the user's custom input values.")
                            blendRxn['reference-collider'] = inputRxn['reference-collider']
                            newColliders = [col for col in inputRxn['colliders']
                                            if col['composition'] in list(self.species_dict.values())]
                            blendRxn['colliders'] = newColliders
                            # blendRxn['colliders'] = inputRxn['colliders']
                    else:
                        inputRxn['colliders'] = [col for col in inputRxn['colliders']
                                                 if col['composition'] in list(self.species_dict.values())]
                        if inputRxn['colliders']:
                            blendData['reactions'].append(inputRxn)
        return blendData

    def _to_builtin(self, obj):
        if isinstance(obj, dict):
            return {self._to_builtin(k): self._to_builtin(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._to_builtin(i) for i in obj]
        elif hasattr(obj, 'as_dict'):
            return self._to_builtin(obj.as_dict())
        elif hasattr(obj, '__dict__'):
            return self._to_builtin(vars(obj))
        elif hasattr(obj, 'tolist'):  # NumPy arrays or similar
            return obj.tolist()
        else:
            return obj

    def _arrheniusFit(self, temps, eps):
        temps = np.asarray(temps, dtype=float)
        eps = np.asarray(eps, dtype=float)
        y = np.log(eps)
        X = np.column_stack([np.ones_like(temps), np.log(temps), 1 / temps])
        if len(temps) == 2:
            c0, c1 = np.linalg.lstsq(X[:, :2], y, rcond=None)[0]
            A, b, Ea = float(np.exp(c0)), float(c1), 0.0
        else:
            c0, c1, c2 = np.linalg.lstsq(X, y, rcond=None)[0]
            A, b, Ea = round(float(np.exp(c0)), 8), round(float(c1), 8), round(float(-c2 * ct.gas_constant), 8)
        return {'A': A, 'b': b, 'Ea': Ea}

    def _rescaleArrhenius(self, k_ref, k_i):
        A_new = k_i['A'] / k_ref['A']
        b_new = k_i['b'] - k_ref['b']
        Ea_new = k_i['Ea'] - k_ref['Ea']
        return {'A': A_new, 'b': b_new, 'Ea': Ea_new}

    def _colliders(self, mech_rxn, blend_rxn=None, generic=False):
        divisor = 1
        colliders = []
        colliderNames = []
        is_M_N2 = False
        is_M_X = False
        ar_troe_eff = None
        troe_efficiencies = {}
        if mech_rxn.reaction_type == 'falloff-Troe':
            troe_efficiencies = mech_rxn.input_data.get('efficiencies', {})
        elif mech_rxn.reaction_type == 'three-body-linear-Burke':
            for c, col in enumerate(mech_rxn.input_data.get('colliders', {})):
                if c > 0 and col['efficiency']['b'] == 0 and col['efficiency']['Ea'] == 0:
                    troe_efficiencies[col['name']] = col['efficiency']['A']
        for name, val in troe_efficiencies.items():
            comp = self.species_dict[name.upper()]
            # Check if N2 is the reference collider instead of Ar
            if comp == {'Ar': 1} and val != 0 and val != 1:
                is_M_N2 = True
                divisor = 1/val
                ar_troe_eff = val
        for name, val in troe_efficiencies.items():
            comp = self.species_dict[name.upper()]
            if is_M_N2 and comp == {'N': 2} and val != 0 and val != 1:
                is_M_N2 = False
                is_M_X = True
            if comp == {'Ar': 1} and val == 0:
                print(f"> Warning: {mech_rxn} has Ar assumed as reference collider, since params cannot be scaled by the Ar=0 value provided. Please fix.")
        zero_eff_compositions = [self.species_dict[name.upper()] for name, val in troe_efficiencies.items() if val == 0]
        citeStr='Bath gas: ' # collider citations appended to here
        if is_M_N2:
            citeStr += "N2. Citations: "
            if blend_rxn:
                for col in blend_rxn['colliders']:
                    if col['composition']=={'N': 2}:
                        k_ref = self._arrheniusFit(col['temperatures'],col['efficiency'])
                for col in blend_rxn['colliders']:
                    newCol = copy.deepcopy(col)
                    k_i = None
                    #Convert N2:Ar database entry to Ar:N2
                    if newCol['composition'] == {'N': 2}:
                        newCol['composition'] = {'Ar': 1}
                        newCol['name'] = next(k for k, v in self.species_dict.items() if v == newCol['composition'])
                        k_i = self._arrheniusFit([300, 1000, 2000], [1, 1, 1])
                    elif newCol['composition'] in list(self.species_dict.values()):
                        k_i = self._arrheniusFit(col['temperatures'],col['efficiency'])
                    if k_i:
                        if newCol['composition'] in zero_eff_compositions:
                            continue
                        newCol['efficiency'] = self._rescaleArrhenius(k_ref, k_i)
                        citeStr += f"{newCol['name']}: {newCol['note']}; "
                        colliderNames.append(newCol['composition'])
                        colliders.append(newCol)
            # Add troe efficiencies that haven't already been given a value
            for name, val in troe_efficiencies.items():
                comp = self.species_dict[name.upper()]
                already_given = comp in colliderNames
                if not already_given and comp != {'N': 2}:
                    colName = next(k for k, v in self.species_dict.items() if v == comp)
                    citeStr += f"{colName}: present work; "
                    colliderNames.append(self.species_dict[name.upper()])
                    colliders.append({
                        'name': colName,
                        'efficiency': {'A': val, 'b': 0, 'Ea': 0},
                        # 'note': 'present work',
                    })
            if generic:
                for col in self.defaults['generic-colliders']:
                    already_given = col['composition'] in colliderNames
                    if col['composition'] in list(self.species_dict.values()) and not already_given and col['composition'] != {'N': 2} and col['composition'] not in zero_eff_compositions:
                        colName = next(k for k, v in self.species_dict.items() if v == col['composition'])
                        citeStr += f"{colName}: {col['note']}; "
                        colliders.append({
                            'name': colName,
                            'efficiency': {'A': col['efficiency'] / divisor, 'b': 0, 'Ea': 0},
                        })
        elif is_M_X:
            citeStr += "X. Citations: "
            if blend_rxn:
                for col in blend_rxn['colliders']:
                    if col['composition'] in list(self.species_dict.values()) and col['composition'] not in zero_eff_compositions:
                        newCol = copy.deepcopy(col)
                        fitted = self._arrheniusFit(newCol['temperatures'], newCol['efficiency'])
                        fitted['A'] = fitted['A'] / ar_troe_eff
                        newCol['efficiency'] = fitted
                        citeStr += f"{newCol['name']}: {newCol['note']}; "
                        colliderNames.append(newCol['composition'])
                        colliders.append(newCol)
            # Add troe efficiencies that haven't already been given a value
            for name, val in troe_efficiencies.items():
                comp = self.species_dict[name.upper()]
                already_given = comp in colliderNames
                if not already_given:
                    colName = next(k for k, v in self.species_dict.items() if v == comp)
                    citeStr += f"{colName}: present work; "
                    colliderNames.append(comp)
                    colliders.append({
                        'name': colName,
                        'efficiency': {'A': val, 'b': 0, 'Ea': 0},
                    })
            if generic:
                for col in self.defaults['generic-colliders']:
                    already_given = col['composition'] in colliderNames
                    if col['composition'] in list(self.species_dict.values()) and not already_given and col['composition'] not in zero_eff_compositions:
                        colName = next(k for k, v in self.species_dict.items() if v == col['composition'])
                        citeStr += f"{colName}: {col['note']}; "
                        colliders.append({
                            'name': colName,
                            'efficiency': {'A': col['efficiency'] / ar_troe_eff, 'b': 0, 'Ea': 0},
                        })
        else:
            citeStr += "AR. Citations: "
            if blend_rxn:
                # Make reaction-specific colliders wrt Ar and append to collider list
                for col in blend_rxn['colliders']:
                    if col['composition'] in list(self.species_dict.values()) and col['composition'] not in zero_eff_compositions:
                        newCol = copy.deepcopy(col)
                        newCol['efficiency']=self._arrheniusFit(newCol['temperatures'], newCol['efficiency'])
                        citeStr += f"{newCol['name']}: {newCol['note']}; "
                        colliderNames.append(newCol['composition'])
                        colliders.append(newCol)
            # Add troe efficiencies that haven't already been given a value
            for name, val in troe_efficiencies.items():
                comp = self.species_dict[name.upper()]
                # already_given = any(col['name'] == name for col in colliders)
                already_given = comp in colliderNames
                if not already_given and comp != {'Ar': 1}:
                    colName = next(k for k, v in self.species_dict.items() if v == comp)
                    citeStr += f"{colName}: present work; "
                    colliders.append({
                        'name': colName,
                        'efficiency': {'A': val, 'b': 0, 'Ea': 0},
                        # 'note': 'present work',
                    })
                    colliderNames.append(comp)
            if generic:
                for col in self.defaults['generic-colliders']:
                    already_given = col['composition'] in colliderNames
                    if col['composition'] in list(self.species_dict.values()) and not already_given and col['composition'] != {'Ar': 1} and col['composition'] not in zero_eff_compositions:
                        colName = next(k for k, v in self.species_dict.items() if v == col['composition'])
                        citeStr += f"{colName}: {col['note']}; "
                        colliders.append({
                            'name': colName,
                            'efficiency': {'A': col['efficiency'], 'b': 0, 'Ea': 0},
                        })
        return colliders, citeStr

    def _zippedMech(self):
        newReactions = []
        blendRxnNames = [rxn['pes'] for rxn in self.blend['reactions']]
        user_rxn_equation = None
        if self.reaction is not None:
            normalized_eq = self._normalizedUserRxn()
            mech_species = {s.upper() for s in self.mech_obj.species_names}
            eq_no_tb = re.sub(r'\(\+\w+\)', '', normalized_eq)
            tokens = re.split(r'\s*(?:<=>|=>|<=|=|\+)\s*', eq_no_tb)
            eq_species = [t.strip() for t in tokens if t.strip()]
            missing = [s for s in eq_species if s.upper() not in mech_species]
            if missing:
                raise ValueError(
                    f"Reaction '{self.reaction}' contains species not found in the mechanism: "
                    f"{', '.join(missing)}. Please check the equation and try again.")
            user_rxn = {
                "equation": normalized_eq,
                "rate-constant": {"A": 1.0, "b": 0.0, "Ea": 0.0}
            }
            user_rxn_obj = ct.Reaction.from_yaml(yaml.dump(user_rxn, sort_keys=False), self.mech_obj)
            user_rxn_equation = user_rxn_obj.equation
        for i, mech_rxn in enumerate(self.mech_obj.reactions()):
            applyLMRR = False
            if self.reaction is not None:
                if user_rxn_equation == mech_rxn.equation:
                    self.foutName = f"{self.foutName}_{self.reaction}"
                    applyLMRR = True
                else: # just append it as-is
                    d = mech_rxn.input_data
                    if 'note' in d and re.fullmatch(r'\n+', d['note']):
                        mech_rxn.update_user_data({'note': ''})
                    newReactions.append(mech_rxn)
            else:
                applyLMRR = True
            if applyLMRR:
                pDep = False
                # Create the M-collider entry for the pressure-dependent reactions
                if mech_rxn.reaction_type in ['falloff-Troe', 'three-body-pressure-dependent-Arrhenius', 'pressure-dependent-Arrhenius', 'Chebyshev', 'three-body-linear-Burke']:
                    pDep = True
                    if mech_rxn.reaction_type == 'three-body-linear-Burke':
                        d = self._to_builtin(mech_rxn.input_data['colliders'][0]) #use the pdep format given for collider M when rebuilding the reaction
                        d.pop("name")
                    else:
                        d = self._to_builtin(mech_rxn.input_data)
                        d.pop("equation")
                        d.pop("efficiencies", None)
                    d.pop("duplicate", None)
                    d.pop("units", None)
                    if d.get('Troe'):
                        d['Troe'] = dict(d['Troe'])
                    if d.get('low-P-rate-constant'):
                        d['low-P-rate-constant'] = dict(d['low-P-rate-constant'])
                    if d.get('high-P-rate-constant'):
                        d['high-P-rate-constant'] = dict(d['high-P-rate-constant'])
                    colliderM = {'name': 'M'}
                    colliderM.update(dict(d))
                # Leave reactions with explicitly-defined bath gases (e.g., (+AR)) in original format
                if pDep and '(+' in mech_rxn.equation and '(+M)' not in mech_rxn.equation:
                    pDep = False
                if pDep and (self.mech_pes[i] in blendRxnNames or self.allPdep):
                    genericBool = self.allPdep
                    blendRxn = None
                    if self.mech_pes[i] in blendRxnNames:
                        # rxn is specifically covered either in defaults or user input
                        idx = blendRxnNames.index(self.mech_pes[i])
                        blendRxn = self.blend['reactions'][idx]
                    if blendRxn and not genericBool:
                        param_type = "ab initio"
                    elif genericBool and not blendRxn:
                        param_type = "generic"
                    elif blendRxn and genericBool:
                        param_type = "ab initio and generic"
                    colliders, citeStr = self._colliders(mech_rxn, blend_rxn=blendRxn, generic=genericBool)
                    d = self._to_builtin(mech_rxn.input_data)
                    newRxn = {
                        'equation': mech_rxn.equation,
                        **({'duplicate': True} if d.get('duplicate') else {}),
                        **({'units': d['units']} if d.get('units') else {}),
                        'type': 'linear-Burke',
                        'colliders': [colliderM] + colliders,
                    }
                    if 'note' in d and not re.fullmatch(r'\n+', d['note']):
                        newRxn['note'] = d['note'] + " " + citeStr
                    else:
                        newRxn['note'] = citeStr
                    newRxn['note'] = newRxn['note'][:-2] + "."
                    yaml_str = yaml.dump(newRxn, sort_keys=False)
                    newRxn_obj = ct.Reaction.from_yaml(yaml_str, self.mech_obj)
                    newReactions.append(newRxn_obj)
                    print(f"{mech_rxn} {dict(self.mech_pes[i])} converted to LMR-R with {param_type} parameters.")
                else: # just append it as-is
                    d = mech_rxn.input_data
                    if 'note' in d and re.fullmatch(r'\n+', d['note']):
                        mech_rxn.update_user_data({'note': ''})
                    newReactions.append(mech_rxn)
                    if self.reaction is not None:
                        self.skipSave = True
                        if not self.allPdep:
                            print(f"User-provided reaction could not be converted using ab initio parameters. Try enabling generic parameters ('allPdep=True') and rerunning.")
                        else:
                            print(f"User-provided reaction could not be converted.")
        output_data = {
            'thermo': self.mech_obj.thermo_model,
            'kinetics': self.mech_obj.kinetics_model,
            'transport_model': self.mech_obj.transport_model,
            'species': self.mech_obj.species(),
            'reactions': newReactions,
            'name': 'gas',
        }
        return ct.Solution(**output_data)

    def _loadYAML(self, fName):
        with open(fName) as f:
            return yaml.safe_load(f)
        
    def _fixDuplicates(self,fName,mech):
        for _attempt in range(10):
            stderr_buffer = io.StringIO()
            with contextlib.redirect_stderr(stderr_buffer):
                ct.Solution(fName)
                error_msg = stderr_buffer.getvalue()
                if not (error_msg and 'undeclared duplicate' in error_msg.lower()):
                    break
                rxn_numbers = set()
                for line in error_msg.split('\n'):
                    m = re.match(r'\s*Reaction\s+(\d+)', line)
                    if m:
                        rxn_numbers.add(int(m.group(1)) - 1)
                if rxn_numbers:
                    for idx in rxn_numbers:
                        if idx < len(mech['reactions']):
                            mech['reactions'][idx]['duplicate'] = True
                            print(f"  Marked reaction {idx + 1} as duplicate: "
                                f"{mech['reactions'][idx]['equation']}")
                with open(fName, 'w') as outfile:
                    yaml.safe_dump(copy.deepcopy(mech), outfile,
                    default_flow_style=None,
                    sort_keys=False)
    
    def _saveYAML(self):
        fName = f"{self.foutName}.yaml"
        self.output.write_yaml(filename=fName, units=self.units)
        # Resave it to remove formatting inconsistencies
        mech = self._loadYAML(fName)
        # Prevent 'NO' from being misinterpreted as bool in species list
        mech['phases'][0]['species'] = [
            "NO" if str(molec).lower() == "false" else molec
            for molec in mech['phases'][0]['species']
        ]
        for species in mech['species']:
            if str(species['name']).lower() == "false":
                species['name']="NO"
        
        for reaction in mech['reactions']:
            if reaction.get('efficiencies'):
                # Prevent 'NO' from being misinterpreted as bool in efficiencies list found in Troe falloff reactions
                reaction['efficiencies'] = {
                    "NO" if str(key).lower() == "false" else key: reaction['efficiencies'][key]
                    for key in reaction['efficiencies']
                }
            if reaction.get('colliders'):
                for col in reaction['colliders']:
                    if str(col['name']).lower() == "false":
                        col['name']="NO"

        # Prevent 'NO' from being misinterpreted as bool in colliders list for LMRR rxns
        for reaction in mech['reactions']:
            effs = reaction.get('efficiencies')
            if effs:
                reaction['efficiencies'] = {
                    "NO" if str(key).lower() == "false" else key: effs[key]
                    for key in effs
                }
        with open(fName, 'w') as outfile:
            yaml.safe_dump(copy.deepcopy(mech), outfile,
                default_flow_style=None,
                sort_keys=False)
        self._fixDuplicates(fName, mech)