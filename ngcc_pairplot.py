import sys
import pandas as pd
import numpy as np
import re
import scipy.stats as stat

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns

merge_method = "inner"
infile = "catdata_bindE.json"
df = pd.read_json(infile)
df.drop_duplicates(inplace=True)

"""
How to get the free energy of reaction for each of these?
I already have E_binding for bare --> bound
By comparing this relative value to the appropriate base catalyst, that should be enough for dGrxn(bare --> bound)
But for the shift, we would need to groupby on the appropriate bound species
Also need to match all of the catalysts together by matching the columns data_dir and Funcnum

"""

df_O2 = df[df["Bound"] == "O2"].reset_index()
df_O2H = df[df["Bound"] == "O2H"].reset_index()
df_O = df[df["Bound"] == "O"]
df_OH = df[df["Bound"] == "OH"]

for df_i in [df_O2, df_O2H, df_O, df_OH, df]:
    print(df_i.shape)

#df_merge = df_O2.join(df_O2H, on=["data_dir", "Funcnum"], how=merge_method, lsuffix="_O2", rsuffix="_O2H")
df_merge = df_O2.merge(df_O2H, left_on=["data_dir", "Funcnum"], right_on=["data_dir", "Funcnum"], how=merge_method, suffixes=("_O2", "_O2H"))
print(df_merge.shape)


"""
df_merge = df_merge.merge(df_O, on=["data_dir", "Funcnum"], how=merge_method, rsuffix="_O")
print(df_merge.shape)
df_merge = df_merge.merge(df_OH, on=["data_dir", "Funcnum"], how=merge_method, rsuffix="_OH")
print(df_merge.shape)
df_merge.drop_duplicates(inplace=True)
print(df_merge.shape)
"""
