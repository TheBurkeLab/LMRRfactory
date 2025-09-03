import matplotlib.pyplot as plt
import numpy as np
import cantera as ct

def compare(params,eps,temps):
    Trange = np.linspace(temps[0],temps[-1])
    A=params['A']
    beta = params['b']
    Ea=params['Ea']
    ylist=[]
    for T in Trange:
        R = ct.gas_constant
        # R=1.987
        yval = A*T**(beta)*np.exp(-Ea/R/T)
        ylist.append(yval)
    plt.figure()
    plt.plot(Trange,ylist)
    plt.plot(temps,eps,linestyle="None",marker="o")
    plt.savefig("comparison.png")
            


# compare({'A': 0.50157049, 'b': 0.45864512, 'Ea': 1.195280919474187},[6.83,12.0,16.3],[300,1000,2000])

compare({'A': 0.64529544, 'b': 0.42626591, 'Ea': 42.89315626},[6.83,12.0,16.3],[300,1000,2000])











# - equation: H2O2 (+M) <=> 2 OH (+M)
#   type: linear-Burke
#   colliders:
#   - name: M
#     type: falloff
#     low-P-rate-constant: {A: 2.49e+24, b: -2.3, Ea: 48749.0}
#     high-P-rate-constant: {A: 2000000000000.0, b: 0.9, Ea: 48749.0}
#     Troe: {A: 0.43, T3: 1.0e-30, T1: 1.0e+30}
#   - name: N2
#     efficiency: {A: 1.14813005, b: 0.04600896, Ea: -2.924426131866635}
#   - name: CO2
#     efficiency: {A: 89.88391622, b: -0.42797367, Ea: 241.4172480999881}
#   - name: H2O2
#     efficiency: {A: 0.50157049, b: 0.45864512, Ea: 1.195280919474187}
#   - name: H2O
#     efficiency: {A: 0.38598987, b: 0.46879325, Ea: 1.195083369452677}
#   - name: HE
#     efficiency: {A: 0.18754314, b: 0.22413586, Ea: 0.0}
#   - name: H2
#     efficiency: {A: 3.53525842, b: -0.07220041, Ea: 0.0}
#   - name: CO
#     efficiency: {A: 2.06652273, b: 0.03887179, Ea: 0.0}
#   - name: CH4
#     efficiency: {A: 160.86525406, b: -0.48030501, Ea: 0.0}
#   - name: O2
#     efficiency: {A: 1.2, b: 0.0, Ea: 0.0}


# - name: H2O2 (+M)
#   pes: {H: 2, O: 2}
#   reference-collider: AR
#   colliders:
#   - name: N2
#     composition: {N: 2}
#     efficiency: [1.50,1.58,1.63]
#     temperatures: [300,1000,2000]
#     note: Jasper-2022
#   - name: CO2
#     composition: {C: 1, O: 2}
#     efficiency: [5.22,4.14,3.27]
#     temperatures: [300,1000,2000]
#     note: Jasper-2022
#   - name: H2O2
#     composition: {H: 2, O: 2}
#     efficiency: [6.83,12.0,16.3]
#     temperatures: [300,1000,2000]
#     note: Jasper-2022
#   - name: H2O
#     composition: {H: 2, O: 1}
#     efficiency: [5.51,10.2,13.3]
#     temperatures: [300,1000,2000]
#     note: Jasper-2022
#   - name: HE
#     composition: {He: 1}
#     efficiency: [0.827,0.966]            #0.43/0.52, 0.57/0.59
#     temperatures: [750,1500]
#     note: Matsugi-2021
#   - name: H2
#     composition: {H: 2}
#     efficiency: [2.192,2.085]            #1.14/0.52, 1.23/0.59
#     temperatures: [750,1500]
#     note: Matsugi-2021
#   - name: CO
#     composition: {C: 1, O: 1}
#     efficiency: [2.673,2.746]            #1.39/0.52, 1.62/0.59
#     temperatures: [750,1500]
#     note: Matsugi-2021
#   - name: CH4
#     composition: {C: 1, H: 4}
#     efficiency: [6.692,4.797]            #3.48/0.52, 2.83/0.59
#     temperatures: [750,1500]
#     note: Matsugi-2021