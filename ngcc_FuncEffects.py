import sys
import pandas as pd
import numpy as np
import re

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns


import CatO2Df
from pylab import rcParams

font = {'size':14}
mpl.rc('font',**font)
rcParams['figure.figsize'] = 10,6

#catalyst = "mepyr"
#catalyst = "tetry17"
catalyst = "tetry20"
#catalyst = "tetrid"
intermediate = "O2H"

csv_dict = {"mepyr":"mepyr_cycle.csv", "tetrid":"tetrid_cycle.csv", "tetry17":"tetry17_cycle.csv", "tetry20":"tetry20_cycle.csv"}
diff_dict = {"O":"O_diff", "O2H":"O2H_diff", "OH":"OH_diff", "O2":"O2_binding_energy"}


#df = pd.read_csv("gcc_assembled.csv")
#df = df[df["Catalyst_Type"] == catalyst] #specify which catalyst to plot data for

df = pd.read_csv(csv_dict[catalyst])

print("Unique Funcs: ", df.Func_1.unique())
print("Unique Subs: ", df.Sub_1.unique())

cat_dict = {"tetry17":{15:0, 24:1, 25:2, 28:3, 29:4, 30:5, 31:6, 32:7}, "tetry20":{15:0, 24:1, 25:2, 28:3, 29:4, 30:5, 31:6, 32:7}, "tetrid":{30:0, 31:1, 34:2, 35:3, 36:4, 37:5, 38:6}, "mepyr":{18:0, 19:1, 24:2}}
energy_dict = {"tetry17":-0.4793, "tetry20":-0.4793, "tetrid":-0.626, "mepyr":-0.1345}

sub_dict = cat_dict[catalyst]
unmod_binding_energy = energy_dict[catalyst]

#print(df["O2H_diff"])
#print(df["O_diff"])
#print(df["OH_diff"])

func_dict ={'C=O':0, 'OC':1, 'methylamine':2, 'trifluoromethyl':3, 'C':4, 'N':5, 'O':6, 'F':7, 'Cl':8, 'Br':9, 'cyanide':10}
func_labels =['CHO', r'OCH$_3$', r'CH$_2$NH$_2$', r'CF$_3$', r'CH$_3$', r'NH$_2$', 'OH', 'F', 'Cl', 'Br', 'CN']

fg_array = np.zeros((len(sub_dict), len(func_dict)))

for index, row in df.iterrows():
    i = sub_dict[row['Sub_1']]
    j = func_dict[row['Func_1']]
    if intermediate == "O2":
        fg_array[i,j] = row['O2_binding_energy'] - unmod_binding_energy
    else:
        fg_array[i,j] = row[diff_dict[intermediate]]

print(fg_array)

## Heatmap, mapped onto molecular coordinates

# load catalyst geometry
xyz_dict = {"mepyr":"xyz/mepyr_optsp_a0m2.xyz", "tetry17":"xyz/tetry_optsp_a0m2.xyz", "tetry20":"xyz/tetry_optsp_a0m2.xyz", "tetrid":"xyz/tetrid_optsp_a0m2.xyz"}
atoms, coords = ChemData.loadXYZ(xyz_dict)

tail_dict = {"mepyr":list(range(1,7))+list(range(10,14)), "tetry17":[3,4]+list(range(18,26)), "tetry20":[3,4]+list(range(18,26)), "tetrid":[3,4]+list(range(21,29)}

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

a_atom_color = []
for atom in cat_atoms:
        a_atom_color.append(atom_dict[atom])

a_mark = []
for shift in sitefunc_shifts:
    if np.sign(shift) >= 0:
        a_mark.append('o')
    else:
        a_mark.append('v')

for i in range(func_sites):
    ax.scatter(a_coord[i,0], a_coord[i,1],a_coord[i,2],s=np.abs(sitefunc_shifts[i])*size_scale, \
               c=a_atom_color[i],marker=a_mark[i] )

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()


## Functional groups are the y-axis, and substitution locations are the x-axis
## bottommost row is the average for the functional group
#sns.heatmap(data=fg_array, linewidth=0.5, yticklabels=sub_dict.keys(), xticklabels=func_labels, \
#            vmin=-0.5, vmax = 0.1, annot=True, fmt=".2f")

#plt.show()
#plt.savefig("Tetry20_O2HFunc.png", transparent=True, bbox_inches='tight', pad_inches=0.05)



"""
## I don't think I need these ways of representing the data, but I'll keep them for now
#Columns: Tetry-17, Tetry-20, Tetrid, Mepyr
#Rows: O2H, O, OH
og_diff = np.array([[-151.01969,-151.02064,-151.03000,-150.98593], \
                     [-75.22371,-75.24038,-75.23585,-75.19085], \
                     [-75.88447,-75.89104,-75.89589,-75.85187]])

og_O2H_diff = {"tetry":-151.02064, "tetrid":-151.03000, "mepyr":-150.98593}
og_O_diff = {"tetry":-75.24038, "tetrid":-75.23585, "mepyr":-75.19085}
og_OH_diff = {"tetry":-75.89104, "tetrid":-75.89589, "mepyr":-75.85187}

ty_O2H_diff = {18:-151.01969, 21:-151.02064}
ty_O_diff =   {18:-75.22371,  21:-75.24038}
ty_OH_diff =  {18:-75.88447,  21:-75.89104}
## Assign columns of og_diff to each row, depending on the catalyst type
#for df in [df_m, df_tr, df_etr]:
#    df = df.assign(og_O2H_diff = df["Catalyst_Type"].map(og_O2H_diff))
#    df = df.assign(og_O_diff = df["Catalyst_Type"].map(og_O_diff))
#    df = df.assign(og_OH_diff = df["Catalyst_Type"].map(og_OH_diff))
#
#print(df_m["og_O_diff"])
#print(df_etr["og_O_diff"])
#
#for df in [df_ty, df_ety]:
#    df = df.assign(og_O2H_diff = df["Active_Site"].map(ty_O2H_diff))
#    df = df.assign(og_O_diff = df["Active_Site"].map(ty_O_diff))
#    df = df.assign(og_OH_diff = df["Active_Site"].map(ty_OH_diff))
#
#print(df_ty["og_O2H_diff"])
#print(df_ety["og_O2H_diff"])

## Separate Assign call for each oxygen species
#df = df.assign(og_O2H_diff = df["Catalyst_Type"].map(og_O2H_diff))
#df = df.assign(og_O_diff = df["Catalyst_Type"].map(og_O_diff))
#df = df.assign(og_OH_diff = df["Catalyst_Type"].map(og_OH_diff))


#cols = ['O2_binding_energy','O2H_diff', 'O_diff', 'OH_diff']

#df['O2H_diff'] = df['O2H_diff'] - df['O2H_diff'].mean()
#df['O_diff'] = df['O_diff'] - df['O_diff'].mean()
#df['OH_diff'] = df['OH_diff'] - df['OH_diff'].mean()
"""
