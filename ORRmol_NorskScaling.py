import sys
import pandas as pd
import numpy as np
import re
import scipy.stats as stat

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns


def func2Ea(dGO):
    return 1.8*dGO + 1.52 - 0.13


#def Norsky(dGOH,dGO):
#    dGOH *= -1.0 #reversing x-axis
#    dGa = -1.0*func2Ea(dGO)
#    dGb = -1*(-1.0*dGOH - dGO)
#    dGc = -1.0*dGOH
#    xnum,ynum = dGOH.shape
#    dG = np.zeros((xnum,ynum))
#    for i in range(xnum):
#        for j in range(ynum):
#            dG[i,j] = min(dGa[i,j],dGb[i,j],dGc[i,j])
#    return dG

def Norsky(dGOH_in, dGO2H_in):
    dGOH = dGOH_in  # input is R-OH --> R, which is ideally -1.23eV. Here it's R --> R-OH, and the ideal is 1.23 eV, as in OER
    dGOOH = dGO2H_in  # R --> R-O2H  # thinking of dropping the shift, s.t. it's less confusing.
    dGO = (-4.92 - dGOH - dGOOH)/2.0
    xnum,ynum = dGOH.shape
    dG = np.zeros((xnum,ynum))
    for i in range(xnum):
        for j in range(ynum):
            dG[i,j] = max(dGOH[i,j],dGOOH[i,j],dGO[i,j])
    return dG
    


font = {'size':14}
mpl.rc('font',**font)

#save_columns = ["Catalyst", "dGrxn_O2H", "dGrxn_regen"]
save_columns = ["Catalyst", "dGrxn_O2H", "dGrxn_regen"]

#df_ngcc = pd.read_json("~/work/ORRmol/ngcc_func/catdata_dGrxn.json")  # an older version of free energy calcs, likely errors
#df_ngcc = pd.read_json("~/work/ORRmol/ngcc_func/catdata_ngcc_dGrxn.json")  # free energy solely from DFT output.
df_ngcc = pd.read_json("catdata_ngcc_dGrxn.json")  # free energy solely from DFT output.
#df_cycleFe = pd.read_json("~/work/ORRmol/cycleFe_func/cycleFe_catcycle_NorskConv_dGrxn.json")
df_cycleFe = pd.read_json("~/work/ORRmol/cycleFe_func/catdata_cycleFe_dGrxn.json")

norsk = pd.read_csv("NorskReview_Scaling.csv", names=["OH", "OOH"])
norsk["OOH"] -= 4.92
norsk["OH"] *= -1
#xbound = [-1,3]
#ybound = [1.5,5.5]
xbound = [-2.,0.5]
ybound = [-3,0]

# Calc reaction for merging steps 1 and 2, in the way done in Norskov papers
#df_ngcc["dGrxn_tO2H"] = df_ngcc["dGrxn_O2"] + df_ngcc["dGrxn_O2H"]  # for older versions where ngcc O2, O2H hadn't been combined
#df_ngcc["dGrxn_tO2H"] = df_ngcc["dGrxn_O2H"]
#df = pd.concat([df_ngcc[save_columns], df_cycleFe[save_columns]])
df = df_ngcc[save_columns]

# Shift free energy by number of electrons yet transferred
#df["dGrxn_regen"] *= -1  # to have a positive trend line, the literature uses R --> R-OH, the reverse of the last step
#df["dGrxn_tO2H"] += 1.23*4  # 4 electrons to be transferred # but maybe we don't do this kind of confusing shift?

print(df)
#ax = plt.plot([-2,3], [1.2,6.2])

df.replace({"Catalyst":{"tetry-17":"C-1", "tetry-20":"C-2", "tetrid":"B", "mepyr":"A"}}, inplace=True)

#sns.scatterplot(data=df, x="dGrxn_regen", y="dGrxn_O2H", hue="Catalyst", legend="brief", s=50)
sns.scatterplot(data=df, x="dGrxn_regen", y="dGrxn_O2H", hue_order=["A","B","C-1","C-2"], hue="Catalyst", legend="brief", s=50)
sns.scatterplot(data=norsk, x="OH", y="OOH", color="black", s=50, label="111 Metals")
plt.ylabel(r"$\Delta G_{O_{2}H}$")
plt.xlabel(r"$\Delta G_{OH}$")
# if we would prefer to have access to the plot as a variable:
#ax = sns.scatterplot(data=df, x="dGrxn_regen", y="dGrxn_O2H", hue="Catalyst", legend="brief", s=50)
#ax = sns.scatterplot(data=norsk, x="OH", y="OOH", color="black", s=50, label="111 Metals")
#ax.set_ylabel(r"$\Delta G_{OOH}$")
#ax.set_xlabel(r"$\Delta G_{OH}$")

plt.legend(loc=3)
plt.xlim(xbound)
plt.ylim(ybound)

# Heatmap
#xlow = -3.0; xhigh = 1.0
#ylow = -5.0; yhigh = 1.2
npoints = 1000
x = np.linspace(xbound[0], xbound[1], npoints)
y = np.linspace(ybound[0], ybound[1], npoints)
X, Y = np.meshgrid(x, y)
#Generate Z axis data for heatmap
Z = Norsky(X,Y)
#Create heatmap
im = plt.imshow(Z, interpolation='bilinear', origin='lower',
                 cmap='coolwarm',extent=(xbound[0], xbound[1], ybound[0], ybound[1]))

#CBI = plt.colorbar(im, orientation='horizontal', pad=0.08, shrink=0.7)
#CBI.set_label("Activity")
#CBI.locator = mpl.ticker.MaxNLocator(nbins=7)
#CBI.update_ticks()

plt.show()
