'''
i can write the same chemical equation many different ways. 
For example, i) NH3 (+M) <=> NH2 + H (+M), ii) NH3 <=> NH2 + H, 
iii) NH3 (+M) <=> H + NH2 (+M), iv) NH3 <=> H + NH2, 
v) NH2 + H (+M) <=> NH3 (+M), vi) NH2 + H <=> NH3, 
vii) H + NH2(+M) <=> NH3 (+M), viii) H + NH2 <=> NH3. 
Let's say A = 'NH3 (+M) <=> NH2 + H (+M)', and B is a list 
containing 20 different random chemical equations, including one 
of the variations of the equation above. I want to check if A is 
in this list, and return True if A, or any of the equivalent forms 
shown above, are in the list. How to do this with python
'''

import re

def normalize_equation(equation):
    # Remove whitespace and standardize the representation
    equation = equation.replace(" ", "")
    
    # Handle the +M species by ignoring them for normalization
    equation = re.sub(r'\(\+M\)', '', equation)
    
    # Sort the reactants and products
    if '<=>' in equation:
        reactants, products = equation.split('<=>')
        reactants = sorted(reactants.split('+'))
        products = sorted(products.split('+'))
        normalized = '+'.join(reactants) + '<=>' + '+'.join(products)
    else:
        normalized = equation

    return normalized

def check_equation_in_list(A, B):
    normalized_A = normalize_equation(A)
    normalized_B = [normalize_equation(eq) for eq in B]
    
    return normalized_A in normalized_B

# Example usage
A = 'NH3 (+M) <=> NH2 + H (+M)'
B = [
    'NH3 <=> NH2 + H',
    'NH3 (+M) <=> H + NH2 (+M)',
    'NH3 <=> H + NH2',
    # ... more equations ...
]

# Check if A or any equivalent forms are in B
result = check_equation_in_list(A, B)
print(result)  # This will print True or False




'''
i can write the same chemical equation many different ways. 
For example, i) NH3 (+M) <=> NH2 + H (+M), ii) NH3 <=> NH2 + H, 
iii) NH3 (+M) <=> H + NH2 (+M), iv) NH3 <=> H + NH2, 
v) NH2 + H (+M) <=> NH3 (+M), vi) NH2 + H <=> NH3, 
vii) H + NH2(+M) <=> NH3 (+M), viii) H + NH2 <=> NH3. Let's say 
A = 'NH3 (+M) <=> NH2 + H (+M)', and chemical mechanism B.yaml 
contains many different random chemical equations, including one 
of the variations of the equation above. I want to check if A is 
in mechanism B, and return True if A, or any of the equivalent forms 
shown above, are in the mechanism. How to do this with python, making 
using of a Cantera Solution object, i.e. 
gas = ct.Solution("mechanism.yaml"). I want a generalized approach that 
works for equations besides just 'NH3 (+M) <=> NH2 + H (+M)',
'''

import cantera as ct

def normalize_equation(equation):
    # Remove spaces and normalize '(+M)'
    equation = equation.replace(' ', '')
    equation = equation.replace('(+M)', '').replace('M', '')
    
    # Split reactants and products
    reactants, products = equation.split('<=>')
    reactants = sorted(reactants.split('+'))
    products = sorted(products.split('+'))
    
    # Return the normalized form as a tuple for easy comparison
    return (tuple(reactants), tuple(products))

def is_equivalent(equation_A, mechanism_file):
    # Normalize the input equation
    normalized_A = normalize_equation(equation_A)
    
    # Load the mechanism using Cantera
    gas = ct.Solution(mechanism_file)

    # Iterate through the reactions in the mechanism
    for reaction in gas.reactions():
        # Normalize the reaction equation
        normalized_reaction = normalize_equation(reaction.equation)
        
        # Check if the normalized forms are equivalent
        if normalized_A == normalized_reaction:
            return True
    
    return False

# Example usage
A = 'NH3 (+M) <=> NH2 + H (+M)'
mechanism_file = 'B.yaml'
result = is_equivalent(A, mechanism_file)
print(result)  # This will print True or False based on the check
