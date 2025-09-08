
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

temps=[300, 1000, 2000]
eps = [5.22, 4.14, 3.27]
def arrhenius_rate(T, A, n, Ea):
    R=1.987
    # R=8.314
    return np.log(A) + n*np.log(T)+ (-Ea/(R*T))

popt, pcov = curve_fit(arrhenius_rate, temps, np.log(eps),maxfev = 3000000)
A_fit, beta_fit, Ea_fit = popt[0], popt[1], popt[2]

def unlog_arrhenius_rate(T, A, beta, Ea):
    # R = 8.314  # Gas constant in J/(mol K)
    R = 1.987 # cal/molK
    return A * T**beta * np.exp(-Ea / (R * T))


print(A_fit)
print(beta_fit)
print(Ea_fit)

Trange=np.linspace(200, 2000, 100)
plt.plot(Trange,unlog_arrhenius_rate(Trange,A_fit,beta_fit,Ea_fit))
plt.plot(temps,eps,linestyle="None",marker="o")
plt.savefig(f"arrheniusfit_test.png")