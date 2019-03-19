import sys
import pandas as pd
import numpy as np
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import seaborn as sns
from sklearn.decomposition import PCA


from autoq import ChemData as CD
import CatO2Df
from pylab import rcParams

font = {'size':14}
mpl.rc('font',**font)
rcParams['figure.figsize'] = 10,6

#catalyst = "mepyr"
#catalyst = "tetry17"
#catalyst = "tetry20"
catalyst = "tetrid"
intermediate = "O2"

csv_dict = {"mepyr":"mepyr_cycle.csv", "tetrid":"tetrid_cycle.csv", "tetry17":"tetry17_cycle.csv", "tetry20":"tetry20_cycle.csv"}
diff_dict = {"O":"O_diff", "O2H":"O2H_diff", "OH":"OH_diff", "O2":"O2_binding_energy"}
atom_dict = {'C':'k','N':'b','F':'r','H':'gray'}
energy_dict = {"tetry17":-0.4793, "tetry20":-0.4793, "tetrid":-0.626, "mepyr":-0.1345}
tail_dict = {"mepyr":list(range(2,6))+list(range(10,14)), "tetry17":list(range(18,26)), "tetry20":list(range(18,26)), "tetrid":list(range(21,29))}
## blend tails (0-index):
blend_dict = {"mepyr":[[2,1],[5,6]], "tetry17":[[18,4], [21,3]], "tetry20":[[18,4], [21,3]], "tetrid":[[21,4], [24,3]]}
xyz_dict = {"mepyr":"xyz/mepyr_optsp_a0m2.xyz", "tetry17":"xyz/tetry_optsp_a0m2.xyz", "tetry20":"xyz/tetry_optsp_a0m2.xyz", "tetrid":"xyz/tetrid_optsp_a0m2.xyz"}
AS_dict = {"mepyr":18, "tetry17":31, "tetry20":33, "tetrid":38}

## Including Br, Cl, and tetrid AS 38
func_dict ={'C=O':0, 'OC':1, 'methylamine':2, 'trifluoromethyl':3, 'C':4, 'N':5, 'O':6, 'F':7, 'Cl':8, 'Br':9, 'cyanide':10}
cat_dict = {"tetry17":{15:0, 24:1, 25:2, 28:3, 29:4, 30:5, 31:6, 32:7}, "tetry20":{15:0, 24:1, 25:2, 28:3, 29:4, 30:5, 31:6, 32:7}, "tetrid":{30:0, 31:1, 34:2, 35:3, 36:4, 37:5, 38:6}, "mepyr":{18:0, 19:1, 24:2}}
cat_funcsite_dict = {"mepyr":[18,19,24], "tetry17":[14,27,28,31,32,33,34,35], "tetry20":[14,27,28,31,32,33,34,35], "tetrid":[30,31,34,35,36,37,38]}

## Excluding some problematic sites and functional groups
#func_dict ={'C=O':0, 'OC':1, 'methylamine':2, 'trifluoromethyl':3, 'C':4, 'N':5, 'O':6, 'F':7, 'cyanide':8}
#cat_dict = {"tetry17":{15:0, 24:1, 25:2, 28:3, 29:4, 30:5, 31:6, 32:7}, "tetry20":{15:0, 24:1, 25:2, 28:3, 29:4, 30:5, 31:6, 32:7}, "tetrid":{30:0, 31:1, 34:2, 35:3, 36:4, 37:5}, "mepyr":{18:0, 19:1, 24:2}}
#cat_funcsite_dict = {"mepyr":[18,19,24], "tetry17":[14,27,28,31,32,33,34,35], "tetry20":[14,27,28,31,32,33,34,35], "tetrid":[30,31,34,35,36,37]}

func_labels =['CHO', r'OCH$_3$', r'CH$_2$NH$_2$', r'CF$_3$', r'CH$_3$', r'NH$_2$', 'OH', 'F', 'Cl', 'Br', 'CN']
size_scale = 10000
active_site = AS_dict[catalyst]

#df = pd.read_csv("gcc_assembled.csv")
#df = df[df["Catalyst_Type"] == catalyst] #specify which catalyst to plot data for

df = pd.read_csv(csv_dict[catalyst])
df = df[df['Func_1'] != 'Cl']
df = df[df['Func_1'] != 'Br']
#df = df[(df['Sub_1'] != active_site) | (df['Func_1'] != 'Cl')]
#df = df[(df['Sub_1'] != active_site) | (df['Func_1'] != 'Br')]

print("Unique Funcs: ", df.Func_1.unique())
df_subs = df.Sub_1.unique()
df_subs.sort()
print("Unique Subs: ", df_subs)

sub_dict = cat_dict[catalyst]
unmod_binding_energy = energy_dict[catalyst]
df['O2_binding_energy'] -= unmod_binding_energy
diff_interm = diff_dict[intermediate]

#print(df["O2H_diff"])
#print(df["O_diff"])
#print(df["OH_diff"])

fg_array = np.zeros((len(sub_dict), len(func_dict)))

for index, row in df.iterrows():
    i = sub_dict[row['Sub_1']] #rows --> modified sites
    j = func_dict[row['Func_1']] #columns --> functional group
    #fg_array[i,j] = row['O2_binding_energy'] - unmod_binding_energy #diff from bare catalyst
    fg_array[i,j] = row[diff_interm]

print("FuncGroup Array:")
print(fg_array)

## Heatmap, mapped onto molecular coordinates

# load catalyst geometry
atoms, coords = CD.loadXYZ(xyz_dict[catalyst])

# tail_dict lists the H's and C's that are part of the "graphitic" region of the catalyst
tail = tail_dict[catalyst]

# func_sites are the H indices that get replaced with functional groups
func_sites = cat_funcsite_dict[catalyst]

# cat_atoms are the atoms that don't get functionalized, but aren't in the tail
# remove all indices that are in the union of tail and funcsites
cat_atom_inds = [i for i in range(len(atoms)) if i not in func_sites+tail]

# sitefunc_shifts measures the effect of functional groups on intermediate energies by site
#sitefunc_shifts = fg_array[:,1]
#sitefunc_shifts = np.mean(np.abs(fg_array), axis=1)
#print(sitefunc_shifts)

## Return the average of the absolute energy of all unique sites
shift_list = []
interm_num = len(diff_dict.keys())

for site in df_subs:
    #shift_list.append(np.mean(np.abs(df[df.Sub_1 == site][diff_interm])))

    partial_shift = 0.
    for diff_interm in diff_dict:
        partial_shift += np.mean(np.abs(df[df.Sub_1 == site][diff_dict[diff_interm]]))
    shift_list.append(partial_shift/interm_num)
sitefunc_shifts = np.array(shift_list)

# Do I just directly PCA the FuncGroup-Site array?
#pca = PCA(n_components=1)
#pca.fit(fg_array.T)
#fg_pca = pca.transform(fg_array)
#sitefunc_shifts = fg_pca
#print("Transformed Array:")
#print(fg_pca)
#plt.scatter(fg_pca[:,0], fg_pca[:,1])
#print(pca.components_)
#print(pca.explained_variance_)
#print(np.sum(np.array(pca.explained_variance_)**0.5))

#pca = PCA().fit(fg_array)
#plt.plot(np.cumsum(pca.explained_variance_ratio_))

#plt.show()

#u, s, vh = np.linalg.svd(fg_array)
#print(u)
#print(s)
#print(vh)

bond_cut = 1.55
n = len(cat_atom_inds)
bonds = []
bonded_atoms = cat_atom_inds + func_sites
for index,i in enumerate(bonded_atoms[:-1]):
    #for j in range(i+1,n):
    for j in bonded_atoms[index+1:]:
        if np.linalg.norm(coords[i,:]-coords[j,:]) < bond_cut:
            bonds.append((i,j))


a_atom_color = []
for i in cat_atom_inds:
        a_atom_color.append(atom_dict[atoms[i]])

shift_marks = []
shift_colors = []
for shift in sitefunc_shifts:
    if np.sign(shift) >= 0:
        shift_marks.append('v')
        shift_colors.append('r')
    else:
        shift_marks.append('*')
        shift_colors.append('purple')



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i,ind in enumerate(cat_atom_inds):
    ax.scatter(coords[ind,0], coords[ind,1],coords[ind,2],s=1000, \
               c=a_atom_color[i],marker="o")

for i,ind in enumerate(func_sites):
    if ind == active_site:
        ax.scatter(coords[ind,0], coords[ind,1],coords[ind,2],s=np.abs(sitefunc_shifts[i])*size_scale, \
                   c='purple',marker=shift_marks[i] )
    else:
        ax.scatter(coords[ind,0], coords[ind,1],coords[ind,2],s=np.abs(sitefunc_shifts[i])*size_scale, \
                   c=shift_colors[i],marker=shift_marks[i] )

for bond in bonds:
    r1 = coords[bond[0],:]
    r2 = coords[bond[1],:]
    plt.plot([r1[0], r2[0]], [r1[1],r2[1]], [r1[2],r2[2]], '-', color='k', lw=6, zorder=1)

blend_bonds = blend_dict[catalyst]
for bond in blend_bonds:
    r1 = coords[bond[0],:]
    r2 = coords[bond[1],:]
    plt.plot([r1[0], r2[0]], [r1[1],r2[1]], [r1[2],r2[2]], '--', color='k', lw=6, zorder=1)

ax.set_axis_off()
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

"""
# Functional groups are the y-axis, and substitution locations are the x-axis
# bottommost row is the average for the functional group
sns.heatmap(data=fg_array, linewidth=0.5, yticklabels=sub_dict.keys(), xticklabels=func_labels, \
            vmin=-0.5, vmax = 0.1, annot=True, fmt=".2f")

plt.show()
#plt.savefig("Tetry20_O2HFunc.png", transparent=True, bbox_inches='tight', pad_inches=0.05)
"""



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
