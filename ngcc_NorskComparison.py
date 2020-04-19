
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

indir = "~/work/CRT/autoq/"
df_m = pd.read_csv(indir+"mepyr_cycle.csv")
df_tr = pd.read_csv(indir+"tetrid_cycle.csv")
df_ty17 = pd.read_csv(indir+"tetry17_cycle.csv")
df_ty20 = pd.read_csv(indir+"tetry20_cycle.csv")

df_ty17["Catalyst_Type"] = "Tetry, Active Site 1"
df_ty20["Catalyst_Type"] = "Tetry, Active Site 2"
df_m["Catalyst_Type"] = "Mepyrid"
df_tr["Catalyst_Type"] = "Tetrid"

## Merge chosen dataframes
df = pd.concat([df_ty17, df_ty20, df_tr, df_m], sort=False)

## Cl and Br would sometimes abstract H from a carbon, creating spuriously negative binding energies
df = df[df['Func_1'] != 'Cl']
df = df[df['Func_1'] != 'Br']

# Take difference from Catalyst_Energy
diff_names = ['O2H', 'O', 'OH']
energy_names = ['CatalystOOH_Energy', 'CatalystO_Energy', 'CatalystOH_Energy']
#df = CatO2Df.AddEnergyDiffs(df, energy_names, diff_names)

#df = CatO2Df.AddReactionEnergiesORR(df, pyc_shifts[1:])
df = CatO2Df.AddReactionEnergiesORR(df)

# Calc reaction for merging steps 1 and 2, in the way done in Norskov papers
df = df.assign(bare_to_OOH = df["O2_binding_energy"] + df["O2_to_OOH"])
df["bare_to_OOH"] -= df["bare_to_OOH"].mean()
df["OH_to_bare"] -= df["OH_to_bare"].mean()


df["OH_to_bare"] *= -1. # Reverse reaction, such that it's R --> R-OH

#norsk = np.loadtxt("NorskReview_Scaling.csv", delimiter=",")
norsk = pd.read_csv("NorskReview_Scaling.csv", names=["OH", "OOH"])
ax = sns.scatterplot(data=norsk, x="OH", y="OOH", color="black", s=50)

# plot scaling relation
ax = sns.scatterplot(data=df, x="OH_to_bare", y="bare_to_OOH", hue="Catalyst_Type", legend=False, s=50)
ax.set_ylabel(r"$\Delta G_{OOH}$")
ax.set_xlabel(r"$\Delta G_{OH}$")



plt.show()















#
