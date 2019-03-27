import sys
import pandas as pd
import numpy as np
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import seaborn as sns
from sklearn.decomposition import PCA


from autoq import ChemData as CD
import CatO2Df
from pylab import rcParams

font = {'size':14}
mpl.rc('font',**font)
rcParams['figure.figsize'] = 10,6

#catalyst = "tetrid"
#catalyst = "mepyr"
#catalyst = "tetry17"
#catalyst = "tetry20"
catalyst = "all"

csv_dict = {"mepyr":"mepyr_cycle.csv", "tetrid":"tetrid_cycle.csv", "tetry17":"tetry17_cycle.csv", "tetry20":"tetry20_cycle.csv"}
func_dict ={'C=O':0, 'OC':1, 'methylamine':2, 'trifluoromethyl':3, 'C':4, 'N':5, 'O':6, 'F':7, 'Cl':8, 'Br':9, 'cyanide':10}
xyz_dict = {"mepyr":"xyz/mepyr_optsp_a0m2.xyz", "tetry17":"xyz/tetry_optsp_a0m2.xyz", "tetry20":"xyz/tetry_optsp_a0m2.xyz", "tetrid":"xyz/tetrid_optsp_a0m2.xyz"}
AS_dict = {"mepyr":18, "tetry17":28, "tetry20":30, "tetrid":38}
energy_dict = {"tetry17":-0.4793, "tetry20":-0.4793, "tetrid":-0.626, "mepyr":-0.1345}

#func_labels =['CHO', r'OCH$_3$', r'CF$_3$',  r'CH$_3$', r'NH$_2$', 'OH', 'F', 'Cl', 'Br', 'CN', r'CH$_2$NH$_2$']
func_labels =['Br', r'CH$_3$', 'CHO', 'Cl', 'F', r'NH$_2$', 'OH', r'OCH$_3$', 'CN', r'CH$_2$NH$_2$', r'CF$_3$']
mstd_names = ['O2_mean', 'O2H_mean', 'O_mean', 'OH_mean', 'O2_std', 'O2H_std', 'O_std', 'OH_std']

def df_setup(catalyst):
    active_site = AS_dict[catalyst]
    df = pd.read_csv(csv_dict[catalyst])
    print(df.Sub_1.unique())
    #df = df[df['Func_1'] != 'Cl']
    #df = df[df['Func_1'] != 'Br']
    df = df[df['Sub_1'] != active_site]
    unmod_binding_energy = energy_dict[catalyst]
    df['O2_binding_energy'] -= unmod_binding_energy
    return df


if catalyst == "all":
    df_list = []
    for catalyst in AS_dict.keys():
        df_list.append(df_setup(catalyst))
    df = pd.concat(df_list, sort=True)

else:
    df = df_setup(catalyst)

#print(df[['Catalyst_File_Name', 'Sub_1', 'O2H_diff']])

unique_funcs = df.Func_1.unique()
print("Unique Funcs: ", unique_funcs)
func_dict = {}
for func in unique_funcs:
    df_func = df[df.Func_1 == func]
    df_func_mean = df_func.mean(axis=0, numeric_only=True)
    df_func_std = df_func.std(axis=0, numeric_only=True)
    print(func)
    func_mean = df_func_mean[['O2_binding_energy', 'O2H_diff', 'O_diff', 'OH_diff']]
    func_std = df_func_std[['O2_binding_energy', 'O2H_diff', 'O_diff', 'OH_diff']]
    func_dict[func] = list(func_mean) + list(func_std)

df_mean_std = pd.DataFrame.from_dict(func_dict, orient='index', columns=mstd_names)
df_mean_std = df_mean_std.sort_index()

print(df_mean_std)


fig, ax = plt.subplots()

#ax.bar(x=func_labels, height=df_mean_std[mstd_names[1]], yerr = df_mean_std[mstd_names[5]], capsize=3.)

width = 0.2
ind = np.arange(len(func_labels))
for i in range(4):
    #ax.errorbar(ind+i*width, df_mean_std[mstd_names[i]], yerr = df_mean_std[mstd_names[i+4]], capsize=3., fmt='o')
    ax.bar(x=ind+i*width, height=df_mean_std[mstd_names[i]], width=width, yerr = df_mean_std[mstd_names[i+4]], capsize=3.)


ax.set_xticks(ind + 1.5*width)
ax.set_xticklabels(func_labels)


ax.axhline(color='k')
plt.ylabel('Average Effect (eV)')
#ax.bar(x=df_mean_std.index, height=df_mean_std[mstd_names[1]], yerr = df_mean_std[mstd_names[5]])
#ax.errorbar(df_mean_std.index, df_mean_std[mstd_names[1]], yerr = df_mean_std[mstd_names[5]], fmt='o')




#plt.show()
plt.savefig("AllCat_FuncEffect_All.png", transparent=True, bbox_inches='tight', pad_inches=0.05)
