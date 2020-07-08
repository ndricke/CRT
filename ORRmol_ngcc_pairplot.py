import sys
import pandas as pd
import numpy as np
import re
import scipy.stats as stat

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns

# TODO this is missing bridge_O (and the other bridge data)

#cols = ["dGrxn_O2", "dGrxn_O2H", "dGrxn_O", "dGrxn_OH", "dGrxn_regen"]
cols = ["dGrxn_O2_and_O2H", "dGrxn_O", "dGrxn_OH", "dGrxn_regen"]

# tetry_AS active site 20, O_intermediate break apart: 2, 4, 5, 6, 7, 

df = pd.read_csv("catdata_dGrxn.csv", index_col=0)

df.loc[(df.Catalyst == "tetry") & (df.Bound_site == 17.0), 'Catalyst'] = "tetry-17"
df.loc[(df.Catalyst == "tetry") & (df.Bound_site == 20.0), 'Catalyst'] = "tetry-20"
df["dGrxn_O2_and_O2H"] = df["dGrxn_O2"] + df["dGrxn_O2H"]


#df_mepyr = df[df["Catalyst"] == "mepyr"]

#ax = sns.kdeplot(df[df["Catalyst_Type"]=="Mepyrid"]["O2_binding_energy"], legend=False, color="r", shade=True)
pair_plt = sns.pairplot(data=df, vars=cols, hue="Catalyst", height=2.5)
plt.show()
