import sys
import pandas as pd
import numpy as np
import re
import scipy.stats as stat

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns

cols = ["dGrxn_O2", "dGrxn_O2H", "dGrxn_O", "dGrxn_OH", "dGrxn_regen"]

df = pd.read_csv("catdata_dGrxn.csv", index_col=0)

for intermediate in ["O2", "O2H", "O", "OH"]:
    df = df[df["GeometryConverged_"+intermediate] == True]

#df_mepyr = df[df["Catalyst"] == "mepyr"]

#ax = sns.kdeplot(df[df["Catalyst_Type"]=="Mepyrid"]["O2_binding_energy"], legend=False, color="r", shade=True)
pair_plt = sns.pairplot(data=df, vars=cols, hue="Catalyst", height=2.5)
plt.show()
