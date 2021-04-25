import sys
import pandas as pd
import numpy as np
import re
import scipy.stats as stat

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns

"""
TODO:
Assess how linear each of these plots are (I think I did this in an earlier instantiation?)
"""

font = {'size':14}
mpl.rc('font',**font)


cols = ["dGrxn_O2", "dGrxn_O2H", "dGrxn_O", "dGrxn_OH", "dGrxn_regen"]
replacements = {"dGrxn_O2":r"$R \rightarrow R-O_2$",
                "dGrxn_O2_and_O2H":r"$R \rightarrow R-O_2H$",
                "dGrxn_O2H":r"$R-O_2 \rightarrow R-O_2H$",
                "dGrxn_O": r"$R-O_2H \rightarrow R-O$",
                "dGrxn_OH": r"$R-O \rightarrow R-OH$",
                "dGrxn_regen": r"$R-OH \rightarrow R$"}

tetry_only = False
combine_O2_O2H = True

# tetry_AS active site 20, O_intermediate break apart: 2, 4, 5, 6, 7, 

df = pd.read_json("catdata_dGrxn.json")


# There are 58 bridge O's at site 20, and only 2 at site 17
if tetry_only:
    df = df[df.Catalyst == "tetry"]
    df.loc[(df.Catalyst == "tetry") & (df.Bound_site == 17.0), 'Catalyst'] = "Tetry, site 1"
    df.loc[(df.Catalyst == "tetry") & (df.Bound_site == 20.0), 'Catalyst'] = "Tetry, site 2"
    df.loc[(df.Catalyst == "Tetry, site 1") & (df.bridge_O == True), 'Catalyst'] = "Tetry, site 1, epoxide bridge"
    df.loc[(df.Catalyst == "Tetry, site 2") & (df.bridge_O == True), 'Catalyst'] = "Tetry, site 2, epoxide bridge"
else:
    df.loc[df.Catalyst == "tetry", 'Catalyst'] = "Tetry"
    df.loc[df.Catalyst == "tetrid", 'Catalyst'] = "Tetrid"
    df.loc[df.Catalyst == "mepyr", 'Catalyst'] = "Mepyr"


if combine_O2_O2H:
    df["dGrxn_O2_and_O2H"] = df["dGrxn_O2"] + df["dGrxn_O2H"]
    cols = ["dGrxn_O2_and_O2H", "dGrxn_O", "dGrxn_OH", "dGrxn_regen"]


#df_mepyr = df[df["Catalyst"] == "mepyr"]

#ax = sns.kdeplot(df[df["Catalyst_Type"]=="Mepyrid"]["O2_binding_energy"], legend=False, color="r", shade=True)
pair_plt = sns.pairplot(data=df, vars=cols, hue="Catalyst", height=2.5)


print(df[cols].shape)
#df_cov = np.cov(np.array(df[cols]).T)
df_cov = np.corrcoef(np.array(df[cols]).T)
print(cols)
print(df_cov)


for i in range(len(cols)):
    for j in range(len(cols)):
        #pair_plt.axes[i,j].set_ylim([-1,1])
        #pair_plt.axes[i,j].set_xlim([-1,1])
        xlabel = pair_plt.axes[i][j].get_xlabel()
        ylabel = pair_plt.axes[i][j].get_ylabel()
        if xlabel in replacements.keys():
            pair_plt.axes[i][j].set_xlabel(replacements[xlabel])
        if ylabel in replacements.keys():
            pair_plt.axes[i][j].set_ylabel(replacements[ylabel])

pair_plt._legend.set_title("Catalysts")
#plt.show()
plt.savefig("ORRmol_ngcc_catalysts4step_pairplot.png", transparent=True, bbox_inches='tight', pad_inches=1.02)
#plt.savefig("ORRmol_ngcc_tetry4step.png", transparent=True, bbox_inches='tight', pad_inches=1.02)
