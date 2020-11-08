"""
This script generates a plot showing the PCA for the effect of modifying each site, intending to evaluate whether
certain sites are more relevant to modify than others

I think I've already filtered out all the cases where O2 doesn't bind or the catalyst doesn't fall apart
This was nice for the volcano and pair plots, but it probably means there are some gaps for a functional group at a particular site
I can take the set of functional groups that succeeded for all sites, or the set of sites that succeeded for all groups
or (I like this best) I can eliminate either in such a way as to maximize data for the PCA plot itself
I still want to see if I can measure the effect of each functional group in an even-handed way

The atom mapping for tetry is all out of order. 
"""

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


def reform_gform(df_in, catalyst, bound_site=None):
    """ Reorganize gform df to get energy change in reaction intermediates"""
    bound_cycle = ["None", "O2", "O2H", "O", "OH"]
    rxn_steps = ["None_O2", "O2_O2H", "O2H_O", "O_OH", "OH_None"]
    df_c = df_in[df_in["Catalyst"] == catalyst]
    if bound_site != None:
        df_c = df_in[(df_in["Bound_site"] == bound_site) | (df_in["Bound"] == "None")]

    print(df_c)

    rxn_Es = {}
    for i, bound in enumerate(bound_cycle):
        ipp = i+1
        if ipp >= len(bound_cycle):
            ipp -= len(bound_cycle)
        b1 = df_c[df_c["Bound"] == bound_cycle[i]].E.values[0]
        b2 = df_c[df_c["Bound"] == bound_cycle[ipp]].E.values[0]
        rxn_Es[rxn_steps[i]] = b2 - b1
    return rxn_Es


font = {'size':14}
mpl.rc('font',**font)
#rcParams['figure.figsize'] = 20,11
rcParams['figure.figsize'] = 15,11

#catalyst = "mepyr"
#catalyst = "tetrid"
#bound_site = None
#catalyst_site = catalyst

catalyst = "tetry"
bound_site = "17"
#bound_site = "20"
catalyst_site = catalyst+bound_site

Ht2eV = 27.211
pca_n = 3
size_scale = 1000
func_labels =['CHO', r'OCH$_3$', r'CH$_2$NH$_2$', r'CF$_3$', r'CH$_3$', r'NH$_2$', 'OH', 'F', 'Cl', 'Br', 'CN']
diff_dict = {"O":"O_diff", "O2H":"O2H_diff", "OH":"OH_diff", "O2":"O2_binding_energy"}
atom_dict = {'C':'k','N':'b','F':'r','H':'gray'}
tail_dict = {"mepyr":list(range(2,6))+list(range(10,14)), "tetry17":list(range(18,26)), "tetry20":list(range(18,26)), 
                "tetrid":list(range(21,29))}


# blend tails (0-index):
blend_dict = {"mepyr":[[2,1],[5,6]], "tetry17":[[18,4], [21,3]], "tetry20":[[18,4], [21,3]], "tetrid":[[21,4], [24,3]]}
# Location of xyz file for base catalyst
xyz_path = "/home/nricke/work/ORRmol/base_xyz_xy/"
xyz_dict = {"mepyr":xyz_path+"mepyr_optsp_xy_a0m2.xyz", "tetry17":xyz_path+"tetry_optsp_xy_a0m2.xyz", 
            "tetry20":xyz_path+"tetry_optsp_xy_a0m2.xyz", "tetrid":xyz_path+"tetrid_optsp_xy_a0m2.xyz"}
AS_dict = {"mepyr":18, "tetry17":31, "tetry20":33, "tetrid":38}
remove_sites = {"tetrid":{38}, "mepyr":{18}, "tetry17":{31}, "tetry20":{33}}
remove_funcs = {"tetrid":set(), "mepyr":{"N", "Cl", "nitroso"}, "tetry17":set(["C=O"]), "tetry20":set(["C=O"])}

# Plotting settings
marker_dict = {1:"v", -1:"D"}
color_list = ['green', 'magenta', 'goldenrod']

df_unmod = pd.read_csv("~/work/ORRmol/dGform_catalysts/ngcc_gform.csv", index_col=0) # load unmodified catalyst data
df_unmod = df_unmod[df_unmod["Catalyst"] == catalyst]
rxn_energies = reform_gform(df_unmod, catalyst, bound_site)
df = pd.read_json("~/work/ORRmol/ngcc_func/catdata_dGrxn.json")
df = df[df["Catalyst"] == catalyst]
df["None_O2"] = (df["Esolv_O2"] - df["Esolv_bare"] - rxn_energies["None_O2"])*Ht2eV
df["O2_O2H"] = (df["Esolv_O2H"] - df["Esolv_O2"] - rxn_energies["O2_O2H"])*Ht2eV
df["O2H_O"] = (df["Esolv_O"] - df["Esolv_O2H"] - rxn_energies["O2H_O"])*Ht2eV
df["O_OH"] = (df["Esolv_OH"] - df["Esolv_O"] - rxn_energies["O_OH"])*Ht2eV
df["OH_None"] = (df["Esolv_bare"] - df["Esolv_OH"] - rxn_energies["OH_None"])*Ht2eV

# these are all 1-site modifications, so convert the lists to values for ease of use
df.loc[:, "func_loc"] = df.loc_O2.map(lambda x: x[0])
df.loc[:, "func"] = df.func_O2.map(lambda x: x[0])

# re-mapping tetry atoms, as these are different between the base and the derivitized molecules
tetry_remap = {32:35, 24:27, 28:31, 25:28, 15:14, 31:34, 30:33, 29:32} # XXX still need to complete
if catalyst == "tetry":
    df.loc[:, "func_loc"] = df.func_loc.map(tetry_remap)

print("Unique Funcs: ", df.func.unique())
print("Unique Subs: ", df.func_loc.unique())
print(rxn_energies)
print(df[["None_O2", "O2_O2H", "O2H_O", "O_OH", "OH_None"]])


# For the grid of functional groups and substitutions, check where data is missing
# Each sub should have all the funcs, each func should have all the subs
# For each sub, get all the funcs it has
# Compare all the func lists, and print the set of funcs that exist for all subs
# removing func_loc 38, the active site I believe, should resolve this for tetrid
sub_list = []
for func in list(df.func.unique()):
    fl = list(df[df["func"] == func].func_loc.unique())
    fl.sort()
    print(func, fl)
    sub_list.append(fl)

res = list(set.intersection(*map(set, sub_list)))
print()
print(res)
print()
func_list = []
for func_loc in list(df.func_loc.unique()):
    fl = list(df[df["func_loc"] == func_loc].func.unique())
    fl.sort()
    print(func_loc, fl)
    func_list.append(fl)

res = list(set.intersection(*map(set, func_list)))
print()
print(res)

# remove difficult functional groups and substitutions
sub_set = set(df.func_loc.unique())
sub_set -= remove_sites[catalyst_site]
func_set = set(df.func.unique())

print(remove_funcs[catalyst_site])
func_set -= remove_funcs[catalyst_site]

print()
print(sub_set)
print(func_set)

# Generate matrix of functional group effects on reaction step
func_num = len(func_set)
subreact_count = 5  # 5 intermediate steps in the ORR cycle
print(subreact_count, func_num)
fg_array = np.zeros((len(sub_set), subreact_count*func_num))  # rows: substitution sites, columns: (func groups)*(intermediates)

sub_sites = list(sub_set)
for i, sub in enumerate(sub_sites):
    for j, func in enumerate(list(func_set)):
        print(sub, func)
        print(df[(df["func"] == func) & (df["func_loc"] == sub)])
        row = df[(df["func"] == func) & (df["func_loc"] == sub)].iloc[0]
        fg_array[i,j] = row["None_O2"]
        fg_array[i,j+func_num] = row["O2_O2H"]
        fg_array[i,j+2*func_num] = row["O2H_O"]
        fg_array[i,j+3*func_num] = row["O_OH"]
        fg_array[i,j+4*func_num] = row["OH_None"]


print("FuncGroup Array:")
print(fg_array)

test_idx = 1
print()
print(fg_array[test_idx,:])
print(df[df["func_loc"] == list(sub_set)[test_idx]][["None_O2", "O2_O2H", "O2H_O", "O_OH", "OH_None"]])


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
atoms, coords = CD.loadXYZ(xyz_dict[catalyst_site])
tail = tail_dict[catalyst_site] # tail_dict lists the H's and C's that are part of the "graphitic" region of the catalyst
# cat_atoms are the atoms that don't get functionalized, but aren't in the tail
# remove all indices that are in the union of tail and the substituted sites
cat_atom_inds = [i for i in range(len(atoms)) if i not in sub_sites+tail]

bond_cut = 1.55
n = len(cat_atom_inds)
bonds = []
bonded_atoms = cat_atom_inds + sub_sites
print(coords)
print(type(coords))
for index,i in enumerate(bonded_atoms[:-1]):
    #for j in range(i+1,n):
    for j in bonded_atoms[index+1:]:
        print(i, j)
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

blend_bonds = blend_dict[catalyst_site]
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
for i,ind in enumerate(sub_sites):
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
        wed = patches.Wedge(coords[ind,:2], r=wsize*0.3, theta1=th1, theta2=th2, zorder=2, color=color_list[j])
        ax.add_patch(wed)

        #(coords[ind,0]+shift[j,0], coords[ind,1]+shift[j,1], s=np.abs(xvhs[i,j])*size_scale, \
        #           c=color_list[j] ,marker=marker_dict[np.sign(xvhs[i,j])], zorder=2 )


fc = 0.6
ax.set_axis_off()
ax.set_xlabel('X')
ax.set_ylabel('Y')
plt.xlim([np.min(coords[:,0])-fc, np.max(coords[:,0])+fc])
plt.ylim([np.min(coords[:,1])-fc, np.max(coords[:,1])+fc])
#plt.show()
if catalyst == "tetry":
    plt.savefig("%s%s_Wedgepca%s.png" % (catalyst, bound_site, pca_n), transparent=True, bbox_inches='tight', pad_inches=0.05)
else:
    plt.savefig("%s_Wedgepca%s.png" % (catalyst, pca_n), transparent=True, bbox_inches='tight', pad_inches=0.05)
