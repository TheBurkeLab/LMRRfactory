import yaml
import copy

def generalized_equations(finName,foutName):
    with open(finName) as f:
        mech = yaml.safe_load(f)
    newMech={'reactions': []}
    for rxn in mech['reactions']:
        eqn = rxn['equation']
        reactants, products = eqn.split('<=>')
        reactants = reactants.strip()
        products = products.strip()
        newMech['reactions'].append(rxn) #append the fwd direction
        rxn = copy.deepcopy(rxn)
        rxn['equation']=f"{products} <=> {reactants}"
        newMech['reactions'].append(rxn) #append the rev direction
        rxn = copy.deepcopy(rxn)
        products = products.replace("(+M)","").strip()
        reactants = reactants.replace("(+M)","").strip()
        rxn['equation']=f"{reactants} <=> {products}"
        newMech['reactions'].append(rxn)
        rxn = copy.deepcopy(rxn)
        rxn['equation']=f"{products} <=> {reactants}"
        newMech['reactions'].append(rxn)
        if "+" in reactants:
            spec1, spec2 = reactants.split("+")
            spec1=spec1.strip()
            spec2=spec2.strip()
            if spec1==spec2:
                reactants = f"2 {spec1}"
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{reactants} <=> {products}"
                newMech['reactions'].append(rxn)
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{products} <=> {reactants}"
                newMech['reactions'].append(rxn)
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
                newMech['reactions'].append(rxn)
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
                newMech['reactions'].append(rxn)
            else:
                reactants = f"{spec2} + {spec1}"
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{reactants} <=> {products}"
                newMech['reactions'].append(rxn)
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{products} <=> {reactants}"
                newMech['reactions'].append(rxn)
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
                newMech['reactions'].append(rxn)
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
                newMech['reactions'].append(rxn)
        elif "+" in products:
            spec1, spec2 = products.split("+")
            spec1=spec1.strip()
            spec2=spec2.strip()
            if spec1==spec2:
                products = f"2 {spec1}"
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{reactants} <=> {products}"
                newMech['reactions'].append(rxn)
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{products} <=> {reactants}"
                newMech['reactions'].append(rxn)
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
                newMech['reactions'].append(rxn)
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
                newMech['reactions'].append(rxn)
            else:
                products = f"{spec2} + {spec1}"
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{reactants} <=> {products}"
                newMech['reactions'].append(rxn)
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{products} <=> {reactants}"
                newMech['reactions'].append(rxn)
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{reactants} (+M) <=> {products} (+M)"
                newMech['reactions'].append(rxn)
                rxn = copy.deepcopy(rxn)
                rxn['equation']=f"{products} (+M) <=> {reactants} (+M)"
                newMech['reactions'].append(rxn)
    with open(foutName, 'w') as outfile:
            yaml.dump(newMech, outfile, default_flow_style=None,sort_keys=False)

bidirectional_input("thirdbodydefaults.yaml","thirdbodydefaults_doubled.yaml")




## IF NON-DUPLICATE SPECIES
# Forwards with (+M): H + OH (+M) <=> H2O (+M)
# Reverse with (+M): H2O (+M) <=> H + OH (+M)
# Forwards without (+M): H + OH <=> H2O
# Reverse without (+M): H2O <=> H + OH
# Swap of Forwards without (+M): OH + H <=> H2O
# Swap of Reverse without (+M): H2O <=> OH + H
# Swap of Forwards with (+M): OH + H (+M) <=> H2O (+M)
# Swap of Reverse with (+M): H2O (+M) <=> OH + H (+M)


## IF DUPLICATE SPECIES
# Forwards with (+M): H2O2 (+M) <=> OH + OH (+M)
# Reverse with (+M): OH + OH (+M) <=> H2O2 (+M)
# Forwards without (+M): H2O2 <=> OH + OH
# Reverse without (+M): OH + OH <=> H2O2
# Merge of Forwards without (+M): H2O2 <=> 2 OH
# Merge of Reverse without (+M): 2 OH <=> H2O2
# Merge of Forwards with (+M): H2O2 (+M) <=> 2 OH (+M)
# Merge of Reverse with (+M): 2 OH (+M) <=> H2O2 (+M)
