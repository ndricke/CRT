import sys
import numpy as np
import pandas as pd


merge_method = "inner"
infile = "catdata_bindE.json"
df = pd.read_json(infile)
df.drop_duplicates(inplace=True)
df.drop(columns=['0', 'index', 'Charge', 'Multiplicity'], inplace=True)

"""
How to get the free energy of reaction for each of these?
I already have E_binding for bare --> bound
By comparing this relative value to the appropriate base catalyst, that should be enough for dGrxn(bare --> bound)
But for the shift, we would need to groupby on the appropriate bound species
Also need to match all of the catalysts together by matching the columns data_dir and Funcnum

"""

df_O2 = df[df["Bound"] == "O2"]
df_O2H = df[df["Bound"] == "O2H"]
df_O = df[df["Bound"] == "O"]
df_OH = df[df["Bound"] == "OH"]

df_O2.loc[df_O2["data_dir"] == "mepyr", "Bound_site"] = np.NaN
df_O2.loc[(df_O2["Catalyst"] == "tetry") & (df_O2["Bound_site"] == 26), "Bound_site"] = 20
df_O2.loc[(df_O2["Catalyst"] == "tetry") & (df_O2["Bound_site"] == 18), "Bound_site"] = 17

for df_i in [df_O2, df_O2H, df_O, df_OH, df]:
    print(df_i.shape)

# add suffix to all columns except ["data_dir", "Funcnum", "Bound_site"] which are preserved for merging
suffix_list = ["_O2", "_O2H", "_O", "_OH"]
keep_same = ["data_dir", "Funcnum", "Bound_site"]
for i, df_intermediate in enumerate([df_O2, df_O2H, df_O, df_OH]):
    df_intermediate.columns = ['{}{}'.format(c, '' if c in keep_same else suffix_list[i]) for c in df_intermediate.columns]


# TODO for some reason there was no _O data in the final database. Figure out why it's not there
df_merge = df_O2.merge(df_O2H, on=["data_dir", "Funcnum", "Bound_site"], how=merge_method)
print(df_merge.shape)
print(df_merge.columns)
df_merge = df_merge.merge(df_O, on=["data_dir", "Funcnum", "Bound_site"], how=merge_method)
print(df_merge.shape)
print(df_merge.columns)
df_merge = df_merge.merge(df_OH, on=["data_dir", "Funcnum", "Bound_site"], how=merge_method)
print(df_merge.shape)
df_merge.drop_duplicates(inplace=True)
print(df_merge.shape)
print(df_merge.columns)
df_merge.to_csv("catdata_bindE_IntMerge.csv")


