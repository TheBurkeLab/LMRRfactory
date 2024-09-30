# import LMRRfactory as LMRR
# from LMRRfactory.masterFitter import masterFitter
import numpy as np
# import sys, os
# sys.path.append("C:/Users/pjsin/Documents/LMRRfactory/build/lib")
from masterFitter import masterFitter

######################################################################

# # INPUTS
T_list=np.linspace(200,2000,100)
# P_list=np.logspace(-12,12,num=120)
P_list=np.logspace(-1,2,num=5)

mF = masterFitter(T_list,P_list,"testinput.yaml",n_P=7,n_T=7,M_only=True)

mF.Troe("LMRtest_Troe_M")
mF.PLOG("LMRtest_PLOG_M")
mF.cheb2D("LMRtest_cheb_M")

# ###################################################
# def makeplot(nom_liste,nom_fig):
#     colours=[(np.random.rand(), np.random.rand(), np.random.rand()) for _ in range(30)]
#     # PLOTTING ACROSS P AT FIXED T 
#     titles=["Reaction_13","Reaction_16","Reaction_25","Reaction_72","Reaction_90","Reaction_167"]
#     indices=[(0,0),(0,1),(1,0),(1,1),(2,0),(2,1)]
#     f, ax = plt.subplots(3, 2, figsize=(10, 10)) 
#     # titles=["Reaction_13","Reaction_16","Reaction_25","Reaction_90"]
#     # indices=[(0,0),(0,1),(1,0),(1,1)]
#     # f, ax = plt.subplots(2, 2, figsize=(10, 10)) 
#     plt.subplots_adjust(wspace=0.2)
#     mpl.rc('font',family='Times New Roman')
#     mpl.rcParams['mathtext.fontset'] = 'stix'
#     mpl.rcParams['font.size'] = 5
#     mpl.rcParams['xtick.labelsize'] = 5
#     mpl.rcParams['ytick.labelsize'] = 5
#     from matplotlib.legend_handler import HandlerTuple
#     plt.rcParams['axes.labelsize'] = 5
#     def get_kTP(fname,P_ls,T_ls,X,reaction,type,linestyle,marker,j,zorder,idx,mkrsz=1.5,reverse=False):
#         gas = ct.Solution(fname)
#         k_TP = []
#         for i,P in enumerate(P_ls):
#             # temp = []
#             # for P in P_ls:
#             gas.TPX = T_ls[0], P*ct.one_atm, {X:1.0}
#             if reverse==True:
#                 k_TP.append(gas.reverse_rate_constants[gas.reaction_equations().index(reaction)])
#             else:
#                 k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
#         ax[idx].loglog(P_ls,k_TP, linestyle=linestyle,linewidth=1,markersize=mkrsz,markeredgewidth=0.6,marker=marker,fillstyle="none",label=f'{X}: {type}',color=colours[j],zorder=zorder)

#     def get_kTPoriginal(fname,P_ls,T_ls,reaction,type,linestyle,marker,j,zorder,idx,mkrsz=1.5,reverse=False):
#         gas = ct.Solution(fname)
#         k_TP = []
#         for i,P in enumerate(P_ls):
#             # temp = []
#             # for P in P_ls:
#             gas.TPX = T_ls[0], P*ct.one_atm, {"AR":1.0}
#             if reverse==True:
#                 k_TP.append(gas.reverse_rate_constants[gas.reaction_equations().index(reaction)])
#             else:
#                 k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
#         ax[idx].loglog(P_ls,k_TP, linestyle=linestyle,linewidth=1,markersize=mkrsz,markeredgewidth=0.6,marker=marker,fillstyle="none",label=f'{type}',color="k",zorder=zorder)

#     T_ls=[1000]

#     nom_PLOG=nom_liste[0]
#     nom_troe=nom_liste[1]
#     nom_cheb=nom_liste[2]
#     # nom_cheb="LMRtest_cheb_M"
#     for j,reaction in enumerate(reactions.keys()):
#         title=titles[j]
#         colliders = reactions[reaction].keys()
#         idx=indices[j]
#         lwr=-1
#         hgr=2
#         fname=f'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\sandbox_substituted.yaml'
#         get_kTPoriginal(fname,np.logspace(lwr,hgr,num=10),T_ls,reaction,"Original","none","o","k",100,idx,mkrsz=4)
#         fname=f'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR.yaml'
#         get_kTPoriginal(fname,np.logspace(lwr,hgr,num=10),T_ls,reaction,"Singal","none","x","k",75,idx,mkrsz=4)
#         for j, X in enumerate(colliders):
#             fname=f'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\{nom_PLOG}.yaml'
#             get_kTP(fname,np.logspace(lwr,hgr,num=60),T_ls,X,reaction.replace(" (+M)",""),"PLOG","-","none",j,1,idx)
#             fname=f'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\{nom_troe}.yaml'
#             get_kTP(fname,np.logspace(lwr,hgr,num=60),T_ls,X,reaction.replace(" (+M)",""),"Troe",":","none",j,50,idx)
#             fname=f'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\{nom_cheb}.yaml'
#             get_kTP(fname,np.logspace(lwr,hgr,num=30),T_ls,X,reaction.replace(" (+M)",""),"Cheb","none","s",j,75,idx)

#         ax[idx].set_title(f"{title}: {reaction} (T=1000K)")
#         ax[idx].legend()
#     plt.savefig('burkelab_SimScripts/rate_constant_plots/'+nom_fig, dpi=1000, bbox_inches='tight')

# makeplot(["LMRtest_PLOG_M","LMRtest_Troe_M","LMRtest_cheb_M"],f'Plog_Troe_Cheb_fixedT.png')
# makeplot(["LMRtest_PLOG_M","LMRtest_Troe_M","LMRtest_cheb_M"],f'Plog_Troe_Cheb_fixedT.svg')
