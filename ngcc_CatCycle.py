
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

indir = "~/work/CRT/autoq/"

df_m = pd.read_csv(indir+"mepyr_cycle.csv")
#df_tr = pd.read_csv(indir+"tetrid_cycle.csv")
#df_ty17 = pd.read_csv(indir+"tetry17_cycle.csv")
#df_ty20 = pd.read_csv(indir+"tetry20_cycle.csv")

df_tr = pd.read_csv(indir+"tetrid_msh_cycle.csv")
df_ty17 = pd.read_csv(indir+"tetry17_msh_cycle.csv")
df_ty20 = pd.read_csv(indir+"tetry20_msh_cycle.csv")

#df_ty17["Catalyst_Type"] = "tetry17"
#df_ty20["Catalyst_Type"] = "tetry20"

df_ty17["Catalyst_Type"] = "Tetry, Active Site 1"
df_ty20["Catalyst_Type"] = "Tetry, Active Site 2"
df_m["Catalyst_Type"] = "Mepyrid"
df_tr["Catalyst_Type"] = "Tetrid"

### Picking which dataframe
df = pd.concat([df_ty17, df_ty20, df_tr, df_m], sort=False)


diff_names = ['O2H', 'O', 'OH']
energy_names = ['CatalystOOH_Energy', 'CatalystO_Energy', 'CatalystOH_Energy']

#Take difference from Catalyst_Energy
#df = CatO2Df.AddEnergyDiffs(df, energy_names, diff_names)

#print(df)
df = df[df['Func_1'] != 'Cl']
df = df[df['Func_1'] != 'Br']


pyc_shifts = (0.5127, 16.3444, -2064.4109, 16.4114, -2064.3416)

#df = CatO2Df.AddReactionEnergiesORR(df, pyc_shifts[1:])
df = CatO2Df.AddReactionEnergiesORR(df)
#df['O2_binding_energy'] = df['O2_binding_energy'] + pyc_shifts[0]


df = df.assign(bare_to_OOH = df["O2_binding_energy"] + df["O2_to_OOH"])
df["bare_to_OOH"] -= df["bare_to_OOH"].mean()

#df["O2_binding_energy"] -= df["O2_binding_energy"].mean()
df["O2_to_OOH"] -= df["O2_to_OOH"].mean()
df["OOH_to_O"] -= df["OOH_to_O"].mean()
df["O_to_OH"] -= df["O_to_OH"].mean()
df["OH_to_bare"] -= df["OH_to_bare"].mean()


#cols = ['O2_binding_energy','O2H_diff', 'O_diff', 'OH_diff'] # for direct intermediate energy correlations
#replacements = {"O2_binding_energy":r"O$_2$", "O2H_diff":r"O$_2$H", "O_diff":"O", "OH_diff":"OH"}

"""
cols = ['O2_binding_energy', "O2_to_OOH", "OOH_to_O", "O_to_OH", "OH_to_bare"] # reaction correlations
replacements = {"O2_binding_energy":r"$R \rightarrow R-O_2$",
                "O2_to_OOH":r"$R-O_2 \rightarrow R-O_2H$",
                "OOH_to_O": r"$R-O_2H \rightarrow R-O$",
                "O_to_OH": r"$R-O \rightarrow R-OH$",
                "OH_to_bare": r"$R-OH \rightarrow R$"}
"""
cols = ["bare_to_OOH", "OOH_to_O", "O_to_OH", "OH_to_bare"] # reaction correlations
replacements = {"bare_to_OOH":r"$R \rightarrow R-O_2H$",
                "OOH_to_O": r"$R-O_2H \rightarrow R-O$",
                "O_to_OH": r"$R-O \rightarrow R-OH$",
                "OH_to_bare": r"$R-OH \rightarrow R$"}


print(df[cols].shape)
#df_cov = np.cov(np.array(df[cols]).T)
df_cov = np.corrcoef(np.array(df[cols]).T)
print(cols)
print(df_cov)

"""
pair_plt = sns.pairplot(data=df, vars=cols, hue="Catalyst_Type", height=2.5)
#sns.pairplot(data=df, vars=cols, hue="Sub_1", height=2.5)
#sns.pairplot(data=df, vars=cols, hue="Func_1", height=2.5)

#sns.pairplot(df[df["Sub_1"] == 19][cols], height=2.5)
#sns.pairplot(df[df["Sub_1"] == 24][cols], height=2.5)

for i in range(len(cols)):
    for j in range(len(cols)):
        pair_plt.axes[i,j].set_ylim([-1,1])
        pair_plt.axes[i,j].set_xlim([-1,1])
        xlabel = pair_plt.axes[i][j].get_xlabel()
        ylabel = pair_plt.axes[i][j].get_ylabel()
        if xlabel in replacements.keys():
            pair_plt.axes[i][j].set_xlabel(replacements[xlabel])
        if ylabel in replacements.keys():
            pair_plt.axes[i][j].set_ylabel(replacements[ylabel])

pair_plt._legend.set_title("Catalysts")
"""


df["OH_to_bare"] *= -1.
#sns.jointplot(data=df, x="OH_to_bare", y="bare_to_OOH") #, hue="Catalyst_Type")
#ax.set_ylabel(r"$R \rightarrow R-O_2H$")
#ax.set_xlabel(r"$R \rightarrow R-OH$")
#ax = sns.scatterplot(data=df, x="OOH_to_O", y="O_to_OH", hue="Catalyst_Type", legend=False, s=50)
#ax.set_ylabel(r"$R-O \rightarrow R-OH$")
#ax.set_xlabel(r"$R-O_2H \rightarrow R-O$")
#plt.scatter()

print(df[df["Catalyst_Type"]=="Tetrid"]["O2_binding_energy"])

ax = sns.kdeplot(df[df["Catalyst_Type"]=="Tetrid"]["O2_binding_energy"], legend=False, color="green", shade=True)
ax = sns.kdeplot(df[df["Catalyst_Type"]=="Mepyrid"]["O2_binding_energy"], legend=False, color="r", shade=True)
ax = sns.kdeplot(df[df["Catalyst_Type"]=="Tetry, Active Site 1"]["O2_binding_energy"], legend=False, color="b", shade=True)
ax = sns.kdeplot(df[df["Catalyst_Type"]=="Tetry, Active Site 2"]["O2_binding_energy"], legend=False, color="orange", shade=True)
ax.set_ylabel("Occurrences")
ax.set_xlabel(r"O$_2$ Binding Energy")

plt.show()
#plt.savefig("TetrMepyr_msh_ReactCycleFunc.png", transparent=True, bbox_inches='tight', pad_inches=1.02)
#plt.savefig("TetrMepyr_norsk_ReactCycleFunc.png", transparent=True, bbox_inches='tight', pad_inches=1.02)
#plt.savefig("TetrMepyr_steps34_ReactCycleFunc.png", transparent=True, bbox_inches='tight', pad_inches=1.02)
#plt.savefig("TetrMepyr_ROOHnorsk.png", transparent=True, bbox_inches='tight', pad_inches=1.02)

#plt.savefig("TetrMepyr_O2density.png", transparent=True, bbox_inches='tight', pad_inches=1.02)















#
