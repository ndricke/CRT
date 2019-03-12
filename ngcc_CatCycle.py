
import sys
import pandas as pd
import numpy as np
import re
import scipy.stats as stat

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns


import CatO2Df

font = {'size':14}
mpl.rc('font',**font)

#df_fname = sys.argv[1]
#df = pd.read_csv(df_fname)

#df_m = pd.read_csv("mepyr1_assembled.csv")
#df_tr = pd.read_csv("tetrid1_assembled.csv")
#df_ty = pd.read_csv("tetry1_assembled.csv")
#df_etr = pd.read_csv("order_tetridEnum_assembled.csv")
#df_ety = pd.read_csv("order_tetryEnum_assembled.csv")

df_m = pd.read_csv("mepyr_cycle.csv")
df_tr = pd.read_csv("tetrid_cycle.csv")
df_ty17 = pd.read_csv("tetry17_cycle.csv")
df_ty20 = pd.read_csv("tetry20_cycle.csv")

df_ty17["Catalyst_Type"] = "tetry17"
df_ty20["Catalyst_Type"] = "tetry20"

### Picking which dataframe
df = pd.concat([df_ty17, df_ty20, df_tr, df_m], sort=False)


diff_names = ['O2H', 'O', 'OH']
energy_names = ['CatalystOOH_Energy', 'CatalystO_Energy', 'CatalystOH_Energy']

#Take difference from Catalyst_Energy
#df = CatO2Df.AddEnergyDiffs(df, energy_names, diff_names)

#print(df)
df = df[df['Func_1'] != 'Cl']
df = df[df['Func_1'] != 'Br']

cols = ['O2_binding_energy','O2H_diff', 'O_diff', 'OH_diff']

print(df[cols].shape)
#df_cov = np.cov(np.array(df[cols]).T)
df_cov = np.corrcoef(np.array(df[cols]).T)
print(cols)
print(df_cov)



#sns.pairplot(data=df, vars=cols, hue="Catalyst_Type", height=2.5)
#sns.pairplot(data=df, vars=cols, hue="Sub_1", height=2.5)
#sns.pairplot(data=df, vars=cols, hue="Func_1", height=2.5)

#sns.pairplot(df[df["Sub_1"] == 19][cols], height=2.5)
#sns.pairplot(df[df["Sub_1"] == 24][cols], height=2.5)
#plt.show()
#plt.savefig("TetrMepyr_CatCycleFunc.png", transparent=True, bbox_inches='tight', pad_inches=1.02)
