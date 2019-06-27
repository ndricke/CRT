import sys
import pandas as pd
import numpy as np
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import patches

import seaborn as sns
from sklearn.decomposition import PCA


from rotate import rotate2D
from autoq import ChemData as CD
import CatO2Df
from pylab import rcParams



font = {'size':14}
mpl.rc('font',**font)
rcParams['figure.figsize'] = 11,11

#catalyst = "mepyr"
#catalyst = "tetry17"
#catalyst = "tetry20"
catalyst = "tetrid"

csv_dict = {"mepyr":"mepyr_cycle.csv", "tetrid":"tetrid_cycle.csv", "tetry17":"tetry17_cycle.csv", "tetry20":"tetry20_cycle.csv"}
diff_dict = {"O":"O_diff", "O2H":"O2H_diff", "OH":"OH_diff", "O2":"O2_binding_energy"}
atom_dict = {'C':'k','N':'b','F':'r','H':'gray'}
energy_dict = {"tetry17":-0.4793, "tetry20":-0.4793, "tetrid":-0.626, "mepyr":-0.1345}
tail_dict = {"mepyr":list(range(2,6))+list(range(10,14)), "tetry17":list(range(18,26)), "tetry20":list(range(18,26)), "tetrid":list(range(21,29))}
## blend tails (0-index):
blend_dict = {"mepyr":[[2,1],[5,6]], "tetry17":[[18,4], [21,3]], "tetry20":[[18,4], [21,3]], "tetrid":[[21,4], [24,3]]}
xyz_dict = {"mepyr":"xyz/mepyr_optsp_xy_a0m2.xyz", "tetry17":"xyz/tetry_optsp_xy_a0m2.xyz", "tetry20":"xyz/tetry_optsp_xy_a0m2.xyz", "tetrid":"xyz/tetrid_optsp_xy_a0m2.xyz"}
AS_dict = {"mepyr":18, "tetry17":31, "tetry20":33, "tetrid":38}
marker_dict = {1:"v", -1:"D"}

## Including Br, Cl, and tetrid AS 38
#func_dict ={'C=O':0, 'OC':1, 'methylamine':2, 'trifluoromethyl':3, 'C':4, 'N':5, 'O':6, 'F':7, 'Cl':8, 'Br':9, 'cyanide':10}
#cat_dict = {"tetry17":{15:0, 24:1, 25:2, 28:3, 29:4, 30:5, 31:6, 32:7}, "tetry20":{15:0, 24:1, 25:2, 28:3, 29:4, 30:5, 31:6, 32:7}, "tetrid":{30:0, 31:1, 34:2, 35:3, 36:4, 37:5, 38:6}, "mepyr":{18:0, 19:1, 24:2}}
#cat_funcsite_dict = {"mepyr":[18,19,24], "tetry17":[14,27,28,31,32,33,34,35], "tetry20":[14,27,28,31,32,33,34,35], "tetrid":[30,31,34,35,36,37,38]}

## Excluding some problematic sites and functional groups
func_dict ={'C=O':0, 'OC':1, 'methylamine':2, 'trifluoromethyl':3, 'C':4, 'N':5, 'O':6, 'F':7, 'cyanide':8}
cat_dict = {"tetry17":{15:0, 24:1, 25:2, 28:3, 29:4, 30:5, 31:6, 32:7}, "tetry20":{15:0, 24:1, 25:2, 28:3, 29:4, 30:5, 31:6, 32:7}, "tetrid":{30:0, 31:1, 34:2, 35:3, 36:4, 37:5}, "mepyr":{18:0, 19:1, 24:2}}
cat_funcsite_dict = {"mepyr":[18,19,24], "tetry17":[14,27,28,31,32,33,34,35], "tetry20":[14,27,28,31,32,33,34,35], "tetrid":[30,31,34,35,36,37]}

func_labels =['CHO', r'OCH$_3$', r'CH$_2$NH$_2$', r'CF$_3$', r'CH$_3$', r'NH$_2$', 'OH', 'F', 'Cl', 'Br', 'CN']
size_scale = 1000
active_site = AS_dict[catalyst]
print(active_site)
pca_n = 3
color_list = ['green', 'magenta', 'goldenrod']

data_dir = "/home/nricke/work/CRT/autoq/"

df = pd.read_csv(data_dir + csv_dict[catalyst])
df = df[df['Func_1'] != 'Cl']
df = df[df['Func_1'] != 'Br']
df = df[df['Sub_1'] != active_site]
#df = df[(df['Sub_1'] != active_site) | (df['Func_1'] != 'Cl')]
#df = df[(df['Sub_1'] != active_site) | (df['Func_1'] != 'Br')]

df_subs = df.Sub_1.unique()
df_subs.sort()
print("Unique Funcs: ", df.Func_1.unique())
print("Unique Subs: ", df_subs)

sub_dict = cat_dict[catalyst]
unmod_binding_energy = energy_dict[catalyst]
df['O2_binding_energy'] -= unmod_binding_energy # so that all values are differences from unmod catalyst
interm_num = len(diff_dict.keys())

#print(df)


"""
Need a consistent way of handling cases where O2 doesn't bind
I could:
1. measure the energy where they are forced together
2. drop them
3. list the binding energy as 0 (in reality there is likely some slightly negative VdW energy)
I plan to go with 3 for the moment, since dropping them would potentially throw out a lot of good data
Actually, if I just drop the active site, that get's rid of a lot of problematic systems
"""

func_num = len(func_dict)
print(interm_num, func_num)
fg_array = np.zeros((len(sub_dict), interm_num*func_num))
for index, row in df.iterrows():
    i = sub_dict[row['Sub_1']] #rows --> modified sites
    j = func_dict[row['Func_1']] #columns --> functional group
    fg_array[i,j] = row["O2_binding_energy"]
    fg_array[i,j+func_num] = row["O2H_diff"]
    fg_array[i,j+2*func_num] = row["O_diff"]
    fg_array[i,j+3*func_num] = row["OH_diff"]

print("FuncGroup Array:")
print(fg_array)

# Do I just directly PCA the FuncGroup-Site array?
#pca = PCA(n_components=6)
#pca.fit(fg_array)
#fg_pca = pca.transform(fg_array)
#sitefunc_shifts = fg_pca
#print("Transformed Array:")
#print(fg_pca)
##plt.scatter(fg_pca[:,0], fg_pca[:,1])
##print(pca.components_)
#print(pca.explained_variance_)
#print(np.sum(np.array(pca.explained_variance_)**0.5))
#
##pca = PCA().fit(fg_array)
#plt.plot(np.cumsum(pca.explained_variance_ratio_))
#plt.show()

u, s, vh = np.linalg.svd(fg_array)
print("s: ", s)


vh_s = vh[:6,:] # taking non-zero singular values
vh_s = vh_s[:pca_n, :] # taking pca_n largest singular values
print("vh_s:")
print(vh_s)
xvhs = np.dot(fg_array, vh_s.T) # principle components projected onto fg_array
print(xvhs)
for i in range(pca_n):
    print("PCA component %s:" %i, np.linalg.norm(xvhs[:,i]))

#np.savetxt("%s_pca%s_vecs.csv" % (catalyst, pca_n), vh_s.T, delimiter=',')
#sys.exit(-1)

#C = np.dot(fg_array.T, fg_array)
#e, v = np.linalg.eigh(C)

## Heatmap, mapped onto molecular coordinates
# load catalyst geometry
atoms, coords = CD.loadXYZ(data_dir + xyz_dict[catalyst])
tail = tail_dict[catalyst] # tail_dict lists the H's and C's that are part of the "graphitic" region of the catalyst
func_sites = cat_funcsite_dict[catalyst] # func_sites are the H indices that get replaced with functional groups
# cat_atoms are the atoms that don't get functionalized, but aren't in the tail
# remove all indices that are in the union of tail and funcsites
cat_atom_inds = [i for i in range(len(atoms)) if i not in func_sites+tail]

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


fig = plt.figure()
ax = fig.add_subplot(111)
#ax = fig.add_subplot(111, projection='3d')

for bond in bonds:
    r1 = coords[bond[0],:]
    r2 = coords[bond[1],:]
    #plt.plot([r1[0], r2[0]], [r1[1],r2[1]], [r1[2],r2[2]], '-', color='k', lw=6, zorder=1)
    plt.plot([r1[0], r2[0]], [r1[1],r2[1]], '-', color='k', lw=6, zorder=1)

blend_bonds = blend_dict[catalyst]
for bond in blend_bonds:
    r1 = coords[bond[0],:]
    r2 = coords[bond[1],:]
    #plt.plot([r1[0], r2[0]], [r1[1],r2[1]], [r1[2],r2[2]], '--', color='k', lw=6, zorder=1)
    plt.plot([r1[0], r2[0]], [r1[1],r2[1]], '--', color='k', lw=6, zorder=1)

for i,ind in enumerate(cat_atom_inds):
    #ax.scatter(coords[ind,0], coords[ind,1],coords[ind,2],s=1000, \
    #           c=a_atom_color[i],marker="o")
    ax.scatter(coords[ind,0], coords[ind,1], s=1000, \
               c=a_atom_color[i], marker="o", zorder=2)



#v1 = np.array([0.3,0,0]).T
v1 = np.array([0.2,0.]).T
#shift = np.zeros((pca_n, 3))
shift = np.zeros((pca_n, 2))
#shift[0,:] = v1
for i in range(pca_n):
    shift[i,:] = rotate2D(v1, i*(2*np.pi/pca_n)).T

#arc_angle = 2*np.pi/pca_n
arc_angle = 360./pca_n
# Plot all the aggregate functional group effects
for i,ind in enumerate(func_sites):
    arcsum = 0.
    for j in range(pca_n):
        #ax.scatter(coords[ind,0]+shift[j,0], coords[ind,1]+shift[j,1], coords[ind,2]+shift[j,2], s=np.abs(xvhs[i,j])*size_scale, \
        #           c='purple' ,marker='v' )
        #ax.scatter(coords[ind,0]+shift[j,0], coords[ind,1]+shift[j,1], s=np.abs(xvhs[i,j])*size_scale, \
        #           c=color_list[j] ,marker=marker_dict[np.sign(xvhs[i,j])], zorder=2 )

        wsize = np.sum((xvhs[i,:])**2)**0.5
        arc = 360. * (xvhs[i,j]**2)/wsize**2
        th1 = arcsum
        th2 = arcsum + arc
        arcsum += arc
        wed = patches.Wedge(coords[ind,:2], r=wsize*0.4, theta1=th1, theta2=th2, zorder=2, color=color_list[j])
        ax.add_patch(wed)

        #(coords[ind,0]+shift[j,0], coords[ind,1]+shift[j,1], s=np.abs(xvhs[i,j])*size_scale, \
        #           c=color_list[j] ,marker=marker_dict[np.sign(xvhs[i,j])], zorder=2 )


ax.set_axis_off()
ax.set_xlabel('X')
ax.set_ylabel('Y')
#ax.set_zlabel('Z')
#plt.show()
plt.savefig("%s_Wedgepca%s.png" % (catalyst, pca_n), transparent=True, bbox_inches='tight', pad_inches=0.05)
