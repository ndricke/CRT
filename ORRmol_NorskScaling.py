import sys
import pandas as pd
import numpy as np
import re
import scipy.stats as stat

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns

font = {'size':14}
mpl.rc('font',**font)

save_columns = ["Catalyst", "dGrxn_tO2H", "dGrxn_regen"]

df_ngcc = pd.read_json("~/work/ORRmol/ngcc_func/catdata_dGrxn.json")
df_cycleFe = pd.read_json("~/work/ORRmol/cycleFe_func/cycleFe_catcycle_NorskConv_dGrxn.json")

# Calc reaction for merging steps 1 and 2, in the way done in Norskov papers
df_ngcc["dGrxn_tO2H"] = df_ngcc["dGrxn_O2"] + df_ngcc["dGrxn_O2H"]

df = pd.concat([df_ngcc[save_columns], df_cycleFe[save_columns]])

# Shift free energy by number of electrons yet transferred
df["dGrxn_regen"] *= -1
df["dGrxn_regen"] += -0.27
df["dGrxn_tO2H"] += 1.23*4 - 0.6

print(df)

#norsk = np.loadtxt("NorskReview_Scaling.csv", delimiter=",")
norsk = pd.read_csv("NorskReview_Scaling.csv", names=["OH", "OOH"])
ax = plt.plot([0,2], [3.2,5.2])


# plot scaling relation
ax = sns.scatterplot(data=df, x="dGrxn_regen", y="dGrxn_tO2H", hue="Catalyst", legend="brief", s=50)
ax = sns.scatterplot(data=norsk, x="OH", y="OOH", color="black", s=50, label="111 Metals")
ax.set_ylabel(r"$\Delta G_{OOH}$")
ax.set_xlabel(r"$\Delta G_{OH}$")
plt.legend(loc=2)
plt.xlim([0,2])
plt.ylim([3,5])
plt.show()
