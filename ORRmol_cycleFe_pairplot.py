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

combine_O2_O2H = False
#cols = ["dGrxn_O2", "dGrxn_O2H", "dGrxn_O", "dGrxn_OH", "dGrxn_regen", "dGrxn_CN", "dGrxn_CO"]
cols = ["dGrxn_tO2H", "dGrxn_O", "dGrxn_OH", "dGrxn_regen", "dGrxn_CN", "dGrxn_CO"]
geo_conv_columns = ["GeometryConverged_bare_O2H", "GeometryConverged_O2H", "GeometryConverged_O", "GeometryConverged_OH", "GeometryConverged_CO", "GeometryConverged_CN"]
replacements = {"dGrxn_O2":r"$R \rightarrow R-O_2$",
                "dGrxn_tO2H":r"$R \rightarrow R-O_2H$",
                "dGrxn_O2H":r"$R-O_2 \rightarrow R-O_2H$",
                "dGrxn_O": r"$R-O_2H \rightarrow R-O$",
                "dGrxn_OH": r"$R-O \rightarrow R-OH$",
                "dGrxn_CO": r"$R \rightarrow R-CO$",
                "dGrxn_CN": r"$R \rightarrow R-CN$",
                "dGrxn_regen": r"$R-OH \rightarrow R$"}

df = pd.read_json("cycleFe_catcycle_GeoConv_CNCO_dGrxn.json")
df = df[df["Filename_O2H"] != "nanFe-functionalized96O2H_optsp_a0m2.out"]

if combine_O2_O2H:
    for conv_column in geo_conv_columns:
        df = df[df[conv_column] == True]
df.dropna(inplace=True)
print(df)

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
plt.show()
#plt.savefig("ORRmol_cycleFes_pairplot.png", transparent=True, bbox_inches='tight', pad_inches=1.02)
