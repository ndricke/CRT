import sys
import pandas as pd
import numpy as np
import re

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns


import CatO2Df

font = {'size':14}
mpl.rc('font',**font)

df_m = pd.read_csv("mepyr1_assembled.csv")
df_tr = pd.read_csv("tetrid1_assembled.csv")
df_ty = pd.read_csv("tetry1_assembled.csv")
df_etr = pd.read_csv("order_tetridEnum_assembled.csv")
df_ety = pd.read_csv("order_tetryEnum_assembled.csv")

### Picking which dataframe
#df = df_m
#df = pd.concat([df_tr, df_etr])
df = pd.concat([df_ty, df_ety])
#df = pd.concat([df_ty, df_ety, df_tr, df_etr, df_m])

diff_names = ['O2H', 'O', 'OH']
energy_names = ['CatalystOOH_Energy', 'CatalystO_Energy', 'CatalystOH_Energy']
#Take difference from Catalyst_Energy
df = CatO2Df.AddEnergyDiffs(df, energy_names, diff_names)
df['O2H_diff'] = df['O2H_diff'] - df['O2H_diff'].mean()
df['O_diff'] = df['O_diff'] - df['O_diff'].mean()
df['OH_diff'] = df['OH_diff'] - df['OH_diff'].mean()
cols = ['O2_binding_energy','O2H_diff', 'O_diff', 'OH_diff']

print(df.Func_1.unique())
print()
print(df.Sub_1.unique())

##Tetry
sub_dict = {15:0, 24:1, 25:2, 28:3, 29:4, 31:5, 32:6}
unmod_binding_energy = -0.4793

##Tetrid
#sub_dict = {30:0, 31:1, 34:2, 35:3, 36:4, 37:5}
#unmod_binding_energy = -0.626

##Mepyr
#sub_dict = {18:0, 19:1, 24:2}
#unmod_binding_energy = -0.1345


func_dict ={'C=O':0, 'OC':1, 'methylamine':2, 'trifluoromethyl':3, 'C':4, 'N':5, 'O':6, 'F':7, 'Cl':8, 'Br':9, 'cyanide':10}

func_labels =['CHO', 'OC', r'CH$_2$NH$_2$', r'CF$_3$', r'CH$_3$', r'NH$_2$', 'OH', 'F', 'Cl', 'Br', 'CN']

fg_array = np.zeros((len(sub_dict), len(func_dict)))

for index, row in df.iterrows():
    i = sub_dict[row['Sub_1']]
    j = func_dict[row['Func_1']]
    #fg_array[i,j] = row['O2_binding_energy'] - unmod_binding_energy
    fg_array[i,j] = row['O2H_diff']

print(fg_array)

##Functional groups are the y-axis, and substitution locations are the x-axis
##bottommost row is the average for the functional group
sns.heatmap(data=fg_array, linewidth=0.5, yticklabels=sub_dict.keys(), xticklabels=func_labels, \
            vmin=-0.2, vmax = 0.2)
plt.show()
#plt.savefig("Tetry_O2HFunc.png", transparent=True, bbox_inches='tight', pad_inches=0.05)
