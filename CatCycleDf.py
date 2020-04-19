
import sys
import pandas as pd
import numpy as np
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import CatO2Df


def intermediate_differences(df):
    diff_names = ['O2H', 'O', 'OH']
    energy_names = ['CatalystOOH_Energy', 'CatalystO_Energy', 'CatalystOH_Energy']
    df = CatO2Df.AssembleDf(df) # calculates O2 binding energy (plus a few features for ML related to O2 binding)
    df = CatO2Df.AddEnergyDiffs(df, energy_names, diff_names) # calculates intermediate - bare for O2H, O, OH

    # These averages should be calculated within catalyst groups
    df['O2H_diff'] = df['O2H_diff'] - df['O2H_diff'].mean()
    df['O_diff'] = df['O_diff'] - df['O_diff'].mean()
    df['OH_diff'] = df['OH_diff'] - df['OH_diff'].mean()

    print(df.columns)
    print(df.Catalyst_Type.unique())

    return df

font = {'size':14}
mpl.rc('font',**font)
datadir = "~/work/CRT/autoq/"

#df_fname = sys.argv[1]
#df = pd.read_csv(df_fname)

df1 = pd.read_csv(datadir+'order_tetrEnum.csv')
df2 = pd.read_csv(datadir+'order_tetrid1.csv')
df3 = pd.read_csv(datadir+'order_tetry1.csv')
df_mepyr = pd.read_csv(datadir+'order_mepyr.csv')

# Add column Catalyst_Type (mepyr, tetry, tetrid)
df1 = CatO2Df.AddTypeName(df1)
df2 = CatO2Df.AddTypeName(df2)
df3 = CatO2Df.AddTypeName(df3)
df_mepyr = CatO2Df.AddTypeName(df_mepyr)

# split df1 into tetrid and tetry components, and merge parts with tetry1 and tetrid1
df1_tetrid = df1[df1['Catalyst_Type'] == 'tetrid']
df1_tetry = df1[df1['Catalyst_Type'] == 'tetry']
df_tetry = pd.concat([df1_tetry, df3])
df_tetrid = pd.concat([df1_tetrid, df2])

df_tetry_diffs = intermediate_differences(df_tetry)
df_tetrid_diffs = intermediate_differences(df_tetrid)
df_mepyr_diffs = intermediate_differences(df_mepyr)
df = pd.concat([df_tetry_diffs, df_tetrid_diffs, df_mepyr_diffs])

cols = ['O2_binding_energy','O2H_diff', 'O_diff', 'OH_diff', 'Catalyst_Type']
sns.pairplot(df[cols], height=2.5, hue='Catalyst_Type')
plt.show()
#plt.savefig("tetrid.png", transparent=True, bbox_inches='tight', pad_inches=0.02)




