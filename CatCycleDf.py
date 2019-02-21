
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

df1 = pd.read_csv('order_tetrEnum.csv')
df2 = pd.read_csv('order_tetrid1.csv')
df3 = pd.read_csv('order_tetry1.csv')

df1 = CatO2Df.AddTypeName(df1)
df2 = CatO2Df.AddTypeName(df2)
df1 = df1[df1['Catalyst_Type'] == 'tetrid']
df = pd.concat([df1,df2])

#df = df1


diff_names = ['O2H', 'O', 'OH']
energy_names = ['CatalystOOH_Energy', 'CatalystO_Energy', 'CatalystOH_Energy']

#Take difference from Catalyst_Energy


print(df.shape)

df = CatO2Df.AssembleDf(df)
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
sns.pairplot(df[cols], height=2.5)
#plt.show()
plt.savefig("tetrid.png", transparent=True, bbox_inches='tight', pad_inches=0.02)




