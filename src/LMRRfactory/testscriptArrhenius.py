
import copy
import numpy as np
from scipy.optimize import least_squares

M=[2.51,1.62,1.20]

eps=[8.97/M[0],8.55/M[1],8.75/M[2]]

temps=np.array([300,1000,2000])
eps = np.array(eps)
def arrhenius_rate(T, A, beta, Ea):
    # R = 8.314  # Gas constant in J/(mol K)
    R = 1.987 # cal/molK
    return A * T**beta * np.exp(-Ea / (R * T))
def fit_function(params, T, ln_eps):
    A, beta, Ea = params
    return np.log(arrhenius_rate(T, A, beta, Ea)) - ln_eps
initial_guess = [3, 0.5, 50.0]
result = least_squares(fit_function, initial_guess, args=(temps, np.log(eps)))
A_fit, beta_fit, Ea_fit = result.x
newEff = {'A': round(A_fit.item(),8),'b': round(beta_fit.item(),8),'Ea': round(Ea_fit.item(),8)}

print(newEff)
