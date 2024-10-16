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
I can write the same chemical equation many different ways. For example, 
i) NH3 (+M) <=> NH2 + H (+M), ii) NH3 <=> NH2 + H, 
iii) NH3 (+M) <=> H + NH2 (+M), iv) NH3 <=> H + NH2, 
v) NH2 + H (+M) <=> NH3 (+M), vi) NH2 + H <=> NH3, 
vii) H + NH2(+M) <=> NH3 (+M), viii) H + NH2 <=> NH3. As another example, 
the following forms are all equivalent and represent the same equation: 
i) H2 (+M) <=> H + H (+M), ii) H2 <=> H + H, iii) H2 (+M) <=> 2 H (+M), 
iv) H2 <=> 2 H, v) H + H (+M) <=> H2 (+M), vi) H + H <=> H2, 
vii) 2 H (+M) <=> H2 (+M), viii) 2 H <=> H2. Let's say 
A = 'NH3 (+M) <=> NH2 + H (+M)' and B = 'H2 <=> 2 H', and chemical mechanism 
C.yaml contains many different random chemical equations, including one of the 
variations of the equation above. I want to check if A and B are in mechanism 
C, and return True if A, or any of the equivalent forms shown above, are in 
the mechanism. How to do this with python, making using of a Cantera Solution 
object, i.e. gas = ct.Solution("mechanism.yaml"). I want a generalized approach 
that works for equations besides just 'NH3 (+M) <=> NH2 + H (+M)' and 'H2 <=> 2 H'. 
I want the normalization function to automatically know that '2 H' is equivalent 
to 'H + H', for example.
'''

import cantera as ct
import re
from collections import Counter

def normalize_equation(equation):
    # Split reactants and products
    reactants, products = equation.split('<=>')
    reactants = reactants.strip().split('+')
    products = products.strip().split('+')
    
    # Helper function to normalize a single species
    def normalize_species(species):
        # Remove colliders and strip whitespace
        species = species.replace('(+M)', '').strip()
        # Count stoichiometry (e.g., '2 H' becomes 'H + H')
        count = re.findall(r'(\d*)\s*(\S+)', species)
        normalized = []
        for c, s in count:
            c = int(c) if c else 1
            normalized.extend([s] * c)  # Expand to the right number of species
        return normalized

    # Normalize reactants and products
    normalized_reactants = sorted(normalize_species(r) for r in reactants)
    normalized_products = sorted(normalize_species(p) for p in products)
    
    # Combine and return as a tuple (reactants, products)
    return (tuple(normalized_reactants), tuple(normalized_products))

def check_equations_in_mechanism(equation_list, mechanism_file):
    gas = ct.Solution(mechanism_file)
    mechanism_equations = [reaction.equation for reaction in gas.reactions()]
    
    # Normalize the equations from the mechanism
    normalized_mechanism = {normalize_equation(eq) for eq in mechanism_equations}
    
    # Check each equation in the list against the normalized mechanism
    for equation in equation_list:
        normalized_equation = normalize_equation(equation)
        if normalized_equation in normalized_mechanism:
            return True
            
    return False

# Example usage
A = 'NH3 (+M) <=> NH2 + H (+M)'
B = 'H2 <=> 2 H'
equations_to_check = [A, B]

# Replace 'C.yaml' with your actual mechanism file
if check_equations_in_mechanism(equations_to_check, 'C.yaml'):
    print("One or more equations are present in the mechanism.")
else:
    print("No matching equations found.")
