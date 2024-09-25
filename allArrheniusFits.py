
# import sys, os
# sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import matplotlib.pyplot as plt
import pandas as pd
import time
import scipy
import scipy.optimize
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import gridspec
from scipy.optimize import least_squares
import sys, os
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import matplotlib as mpl
import argparse

# parser = argparse.ArgumentParser()
# parser.add_argument('--figwidth', type=float, help="figwidth = ")
# parser.add_argument('--figheight', type=float, help="figheight = ")
# parser.add_argument('--fsz', type=float, help="mpl.rcParams['font.size'] = ", default=8)
# parser.add_argument('--fszxtick', type=float, help="mpl.rcParams['xtick.labelsize'] = ", default=7)
# parser.add_argument('--fszytick', type=float, help="mpl.rcParams['ytick.labelsize'] = ", default=7)
# parser.add_argument('--fszaxlab', type=float, help="mpl.rcParams['axes.labelsize'] = ", default=8)
# parser.add_argument('--lw', type=float, help="lw = ", default=0.7)
# parser.add_argument('--mw', type=float, help="mw = ", default=0.5)
# parser.add_argument('--msz', type=float, help="msz = ", default=2.5)
# parser.add_argument('--lgdw', type=float, help="lgdw = ", default=0.6)
# parser.add_argument('--lgdfsz', type=float, help="lgdw = ", default=5)
# parser.add_argument('--gridsz', type=int, help="gridsz = ", default=10)
# parser.add_argument('--dpi', type=int, help="dpi = ", default=1000)

# args = parser.parse_args()
mpl.rc('font',family='Times New Roman')
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.size'] = 6
mpl.rcParams['xtick.labelsize'] = 5
mpl.rcParams['ytick.labelsize'] = 5
from matplotlib.legend_handler import HandlerTuple
plt.rcParams['axes.labelsize'] = 6
mpl.rcParams['xtick.major.width'] = 0.5  # Width of major ticks on x-axis
mpl.rcParams['ytick.major.width'] = 0.5  # Width of major ticks on y-axis
mpl.rcParams['xtick.minor.width'] = 0.5  # Width of minor ticks on x-axis
mpl.rcParams['ytick.minor.width'] = 0.5  # Width of minor ticks on y-axis
mpl.rcParams['xtick.major.size'] = 2.5  # Length of major ticks on x-axis
mpl.rcParams['ytick.major.size'] = 2.5  # Length of major ticks on y-axis
mpl.rcParams['xtick.minor.size'] = 1.5  # Length of minor ticks on x-axis
mpl.rcParams['ytick.minor.size'] = 1.5  # Length of minor ticks on y-axis

save_plots = True
fig, ax = plt.subplots(3,2,figsize=(6,9))
# plt.subplots_adjust(wspace=0.4, hspace=1)
import matplotlib.ticker as ticker

lw=0.8
mw=0.5
msz=2.5
dpi=400
lgdw=0.9
lgdfsz=7

plt.subplots_adjust(hspace=0.25)
xrange=[1,2700]

def plot_ratefit(rxn, effs, pltcolour, idx): #3-parameter fit (A, b, Ea)
    
    temperatures=np.array(rxn['temperatures'])
    # print(temperatures)
    labell=effs['collider']
    epsLow=effs['epsLow']['A']
    epsHigh=effs['epsHigh']['A']
    rate_constants=np.array([epsLow,epsHigh])
    def arrhenius_rate(T, A, beta, Ea):
        # R = 8.314  # Gas constant in J/(mol K)
        R = 1.987 # cal/molK
        return A * T**beta * np.exp(-Ea / (R * T))
    def fit_function(params, T, ln_rate_constants):
        A, beta, Ea = params
        return np.log(arrhenius_rate(T, A, beta, Ea)) - ln_rate_constants
    initial_guess = [3, 0.5, 50.0]  
    result = least_squares(fit_function, initial_guess, args=(temperatures, np.log(rate_constants)))
    A_fit, beta_fit, Ea_fit = result.x
    newDict2={
         'collider': labell,
         'epsLowT': effs['epsLow'],
         'epsHighT': effs['epsHigh'],
         'epsFitted': {'A': round(float(A_fit),3),'b': round(float(beta_fit),3),'Ea': round(float(Ea_fit),3)},
    }
    T_range = np.linspace(200, 2000, 100)
    fit_curve = np.exp(np.log(arrhenius_rate(T_range, A_fit, beta_fit, Ea_fit)))
    ax[idx].semilogy(T_range, fit_curve,linewidth=lw, label=labell, color=pltcolour)
    ax[idx].semilogy(temperatures, rate_constants,marker='o',fillstyle='full',markersize=msz,markeredgewidth=mw,linestyle='none',color=pltcolour,label=None)
    ax[idx].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=5))

    return newDict2

import yaml
idx=(0,0)
def openyaml(fname):
            with open(fname) as f:
                return yaml.safe_load(f)
keyReactions = {}
efficiencies = openyaml('efficiencies.yaml')
pltcolours = ["xkcd:grey", 'black', "xkcd:teal", 'r', 'b', 'xkcd:purple','goldenrod','olive', 'brown']

outDict={
        'reactions': []
    }


idxs = [(0,0),(0,1),(1,0),(1,1),(2,0),(2,1)]
for j, rxn in enumerate(efficiencies['reactions']):
    newDict={
        'equation': rxn['equation'],
        'temperatures':rxn['temperatures'],
        'collider-list': []
    }
    for i, col in enumerate(efficiencies['reactions'][j]['collider-list']):
        # print(efficiencies['reactions'][j])
        newDict2 = plot_ratefit(efficiencies['reactions'][j], col, pltcolours[i], idxs[j])
        ax[idxs[j]].legend(fontsize=lgdfsz, frameon=False,loc='right',handlelength=lgdw) 
        ax[idxs[j]].set_title(efficiencies['reactions'][j]['equation'])
        ax[idxs[j]].tick_params(axis='both', direction="in")
        ax[idxs[j]].tick_params(axis='both', which='minor', direction="in")
        ax[idxs[j]].set_xlim(xrange)
        newDict['collider-list'].append(newDict2)
    outDict['reactions'].append(newDict)

with open('efficiencies_fitted.yaml', 'w') as outfile:
    yaml.dump(outDict, outfile, default_flow_style=None,sort_keys=False)


# xrange=[1,2700]

# T_list=[300,1000,2000]
# plot_ratefit(np.array(T_list), np.array([1/2.51,1/1.62,1/1.20]), pltcolours[7], "Ar",idx)
# plot_ratefit(np.array(T_list), np.array([8.97/2.51,8.55/1.62,8.75/1.20]), pltcolours[4], "H$_2$O",idx) 
# # ax[idx].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
# # ax[idx].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
# # ax[idx].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[idx].legend(fontsize=lgdfsz, frameon=False,loc='right') 
# ax[idx].set_title("$\epsilon_{0,i}/\epsilon_{0,N_2}$ for H + OH (+M) $\leftrightharpoons$ H$_2$O (+M)")
# ax[idx].tick_params(axis='both', direction="in")
# ax[idx].tick_params(axis='both', which='minor', direction="in")
# # ax[idx].set_xlim(xrange)





# idx=(0,1)
# T_list=[300,1000,2000]
# plot_ratefit(np.array(T_list), np.array([13.7,8.94,3.03]), pltcolours[2], "CO$_2$",idx)
# plot_ratefit(np.array(T_list), np.array([3.69,3.07,1.71]), pltcolours[5], "H$_2$",idx)
# plot_ratefit(np.array(T_list), np.array([23.3,22.2,21.3]), pltcolours[4], "H$_2$O",idx)
# plot_ratefit(np.array(T_list), np.array([0.90,1.17,1.34]), pltcolours[6], "He",idx)
# plot_ratefit(np.array(T_list), np.array([1.71,1.58,1.20]), pltcolours[0], "N$_2$",idx)
# plot_ratefit(np.array(T_list), np.array([20.4,17.9,18.7]), pltcolours[3], "NH$_3$",idx)

# # ax[idx].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
# # ax[idx].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
# # ax[idx].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[idx].legend(fontsize=lgdfsz, frameon=False,loc='right') 
# ax[idx].set_title("$\epsilon_{0,i}/\epsilon_{0,Ar}$ for H + O$_2$ (+M) $\leftrightharpoons$ HO$_2$ (+M)")
# ax[idx].tick_params(axis='both', direction="in")
# ax[idx].tick_params(axis='both', which='minor', direction="in")
# ax[idx].set_xlim(xrange)

# idx=(1,0)
# T_list=[300,1000,2000]
# plot_ratefit(np.array([300,1000,2000]), np.array([5.22,4.14,3.27]), pltcolours[2], "CO$_2$",idx)
# plot_ratefit(np.array([300,1000,2000]), np.array([5.51,10.2,13.3]), pltcolours[4], "H$_2$O",idx) 
# plot_ratefit(np.array([300,1000,2000]), np.array([6.83,12.0,16.3]), pltcolours[8], "H$_2$O$_2$",idx)
# plot_ratefit(np.array([300,1000,2000]), np.array([1.50,1.58,1.63]), pltcolours[0], "N$_2$",idx)

# # ax[idx].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
# # ax[idx].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
# # ax[idx].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[idx].legend(fontsize=lgdfsz, frameon=False,loc='right') 
# ax[idx].set_title("$\epsilon_{0,i}/\epsilon_{0,Ar}$ for H$_2$O$_2$ (+M) $\leftrightharpoons$ OH + OH (+M)")
# ax[idx].tick_params(axis='both', direction="in")
# ax[idx].tick_params(axis='both', which='minor', direction="in")
# ax[idx].set_xlim(xrange)
# ax[idx].set_ylabel(r'Relative third-body efficiency ($\epsilon_{0,i}/\epsilon_{0,M}$)')

# idx=(1,1)
# T_list=[300,1000,2000]
# plot_ratefit(np.array([300, 1000, 2000]), np.array([11.2,13.4,14.3]), pltcolours[2], "CO$_2$",idx)
# plot_ratefit(np.array([300, 1000, 2000]), np.array([14.0,23.6,27.9]), pltcolours[4], "H$_2$O",idx)
# plot_ratefit(np.array([300, 1000, 2000]), np.array([3.15,2.47,2.25]), pltcolours[0], "N$_2$",idx)
# plot_ratefit(np.array([300, 1000, 2000]), np.array([13.9,20.0,22.2]), pltcolours[3], "NH$_3$",idx)
# plot_ratefit(np.array([300, 1000, 2000]), np.array([1.55,1.48,1.70]), pltcolours[1], "O$_2$",idx)
# # plot_ratefit(np.array([300, 1000, 2000]), np.array([9.94,13.3,14.3]), pltcolours[4], "CH$_4$",idx)

# # ax[idx].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
# # ax[idx].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
# # ax[idx].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[idx].legend(fontsize=lgdfsz, frameon=False,loc='right') 
# ax[idx].set_title("$\epsilon_{0,i}/\epsilon_{0,Ar}$ for NH$_3$ (+M) $\leftrightharpoons$ H + NH$_2$ (+M)")
# ax[idx].tick_params(axis='both', direction="in")
# ax[idx].tick_params(axis='both', which='minor', direction="in")
# # ax[idx].set_xlim(xrange)


# idx=(2,0)
# T_list=[300,1000,2000]
# plot_ratefit(np.array([300,1000,2000]), np.array([5.60,7.98,8.19]), pltcolours[4], "H$_2$O",idx)
# plot_ratefit(np.array([300,1000,2000]), np.array([2.00,1.69,1.43]), pltcolours[0], "N$_2$",idx)
# plot_ratefit(np.array([300,1000,2000]), np.array([5.86,8.25,8.60]), pltcolours[3], "NH$_3$",idx)
# plot_ratefit(np.array([300,1000,2000]), np.array([1.22,1.17,1.14]), pltcolours[1], "O$_2$",idx)

# # ax[idx].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
# # ax[idx].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
# # ax[idx].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[idx].legend(fontsize=lgdfsz, frameon=False,loc='right') 
# ax[idx].set_title("$\epsilon_{0,i}/\epsilon_{0,Ar}$ for NH$_2$ + NH$_2$ (+M) $\leftrightharpoons$ N$_2$H$_4$ (+M)")
# ax[idx].tick_params(axis='both', direction="in")
# ax[idx].tick_params(axis='both', which='minor', direction="in")
# # ax[idx].set_xlim(xrange)

# idx=(2,1)
# T_list=[300,1000,2000]
# plot_ratefit(np.array([300, 1000, 2000]), np.array([6.42,8.15,9.09]), pltcolours[4], "H$_2$O",idx)
# plot_ratefit(np.array([300, 1000, 2000]), np.array([1.61,1.71,1.59]), pltcolours[0], "N$_2$",idx)

# # ax[idx].set_xlim(xrange)

fig.text(0.5, 0.07, r'Temperature [K]', ha='center', va='center',fontsize=6)
name = f'allArrheniusFits'
if save_plots == True:
    plt.savefig(name+'.png', dpi=500, bbox_inches='tight')