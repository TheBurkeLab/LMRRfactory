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
