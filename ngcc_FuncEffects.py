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

catalyst = "tetrid"

csv_dict = {"mepyr":"mepyr_cycle.csv", "tetrid":"tetrid_cycle.csv", "tetry17":"tetry17_cycle.csv", "tetry20":"tetry20_cycle.csv"}
func_dict ={'C=O':0, 'OC':1, 'methylamine':2, 'trifluoromethyl':3, 'C':4, 'N':5, 'O':6, 'F':7, 'Cl':8, 'Br':9, 'cyanide':10}
xyz_dict = {"mepyr":"xyz/mepyr_optsp_a0m2.xyz", "tetry17":"xyz/tetry_optsp_a0m2.xyz", "tetry20":"xyz/tetry_optsp_a0m2.xyz", "tetrid":"xyz/tetrid_optsp_a0m2.xyz"}
AS_dict = {"mepyr":18, "tetry17":31, "tetry20":33, "tetrid":38}

df = pd.read_csv(csv_dict[catalyst])
#df = df[df['Func_1'] != 'Cl']
#df = df[df['Func_1'] != 'Br']
df = df[df['Sub_1'] != active_site]
#df = df[(df['Sub_1'] != active_site) | (df['Func_1'] != 'Br')]

unique_funcs = df.Func_1.unique()
print("Unique Funcs: ", unique_funcs)

for func in unique_funcs:
    df_func = df[df.Func_1 == func]
    print(func)
    print(df_func)
