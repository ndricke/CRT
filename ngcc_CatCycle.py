
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

#df_fname = sys.argv[1]
#df = pd.read_csv(df_fname)

df_m = pd.read_csv("mepyr1_assembled.csv")
df_tr = pd.read_csv("tetrid1_assembled.csv")
df_ty = pd.read_csv("tetry1_assembled.csv")
df_etr = pd.read_csv("order_tetridEnum_assembled.csv")
df_ety = pd.read_csv("order_tetryEnum_assembled.csv")

### Picking which dataframe
#df = df_m
#df = pd.concat([df_tr, df_etr])
#df = pd.concat([df_ty, df_ety])
df = pd.concat([df_ty, df_ety, df_tr, df_etr, df_m])


diff_names = ['O2H', 'O', 'OH']
energy_names = ['CatalystOOH_Energy', 'CatalystO_Energy', 'CatalystOH_Energy']

#Take difference from Catalyst_Energy
df = CatO2Df.AddEnergyDiffs(df, energy_names, diff_names)

df['O2H_diff'] = df['O2H_diff'] - df['O2H_diff'].mean()
df['O_diff'] = df['O_diff'] - df['O_diff'].mean()
df['OH_diff'] = df['OH_diff'] - df['OH_diff'].mean()

print(df)


#plt.scatter(df['O2_binding_energy'], df['O2H_diff'])
#plt.scatter(df['O2_binding_energy'], df['OH_diff'])
#plt.scatter(df['O2_binding_energy'], df['O_diff'])
#plt.show()

cols = ['O2_binding_energy','O2H_diff', 'O_diff', 'OH_diff']
#sns.pairplot(data=df, vars=cols, hue="Sub_1", height=2.5)
#sns.pairplot(data=df, vars=cols, hue="Catalyst_Type", height=2.5)
sns.pairplot(data=df, vars=cols, hue="Func_1", height=2.5)

#sns.pairplot(df[df["Sub_1"] == 19][cols], height=2.5)
#sns.pairplot(df[df["Sub_1"] == 24][cols], height=2.5)
plt.show()
#plt.savefig("TetrMepyr_CatCycleFunc.png", transparent=True, bbox_inches='tight', pad_inches=1.02)
