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
variations of each equation above. I want to check if A and B are in mechanism 
C, and return True if A or B, or any of their equivalent forms shown above, are 
in the mechanism.. How to do this with python, making using of a Cantera Solution 
object, i.e. gas = ct.Solution("mechanism.yaml"). I want a generalized approach 
that works for equations besides just 'NH3 (+M) <=> NH2 + H (+M)' and 'H2 <=> 2 H'. 
I want the normalization function to automatically know that '2 H' is equivalent 
to 'H + H', for example.
'''

import cantera as ct
import re
from collections import Counter

def normalize_equation(equation):
    """
    Normalize a chemical equation by simplifying coefficients and 
    removing colliders, and ordering species consistently.
    """
    # Split the equation into reactants and products
    reactants, products = equation.split('<=>')
    reactants = reactants.strip().replace('(+M)', '').replace(' ', '')
    products = products.strip().replace('(+M)', '').replace(' ', '')

    # Function to normalize individual side of the equation
    def normalize_side(side):
        # Split into species and their coefficients
        species_list = re.split(r'\s*\+\s*', side)  # Split on '+' with optional whitespace
        species_counter = Counter()

        for species in species_list:
            # Handle cases with coefficients like '2H'
            match = re.match(r'(\d*)(\D+)', species)
            if match:
                coeff, name = match.groups()
                coeff = int(coeff) if coeff else 1  # Default to 1 if no coefficient
                species_counter[name] += coeff
            else:
                species_counter[species] += 1

        # Sort species alphabetically and create a normalized string
        normalized_species = []
        for name in sorted(species_counter.keys()):  # Sort species alphabetically
            count = species_counter[name]
            if count > 1:
                normalized_species.append(f"{count} {name}")
            else:
                normalized_species.append(name)

        return ' + '.join(normalized_species)

    # Normalize both sides and return the formatted equation
    normalized_reactants = normalize_side(reactants)
    normalized_products = normalize_side(products)

    return f"{normalized_reactants} <=> {normalized_products}"

def find_equivalent_reactions(mechanism_file, target_equations):
    """
    Check if any of the target equations or their equivalent forms 
    exist in the Cantera mechanism file.
    """
    # Load the mechanism using Cantera
    gas = ct.Solution(mechanism_file)

    # Extract equations from the mechanism
    mechanism_equations = []
    for rxn in gas.reactions():
        mechanism_equations.append(str(rxn))

    # Normalize the mechanism equations
    normalized_mechanism_equations = {normalize_equation(eqn) for eqn in mechanism_equations}

    # Check for equivalent target equations
    for target in target_equations:
        normalized_target = normalize_equation(target)
        if normalized_target in normalized_mechanism_equations:
            return True  # Found an equivalent reaction

    return False  # No equivalent reactions found

# Example usage:
target_equations = [
    'NH3 (+M) <=> NH2 + H (+M)',
    'H2 <=> 2 H'
]
mechanism_file = 'C.yaml'  # Replace with your mechanism file name

# Call the function to check for equivalents
result = find_equivalent_reactions(mechanism_file, target_equations)
print(result)  # Should return True if any equivalents are found
