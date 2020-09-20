import sys
import numpy as np
import pandas as pd


merge_method = "inner"
infile = "cycleFe_analyzed.json"
df = pd.read_json(infile)
#df.drop(columns=["Species", 'JobType', 'Charge', 'Multiplicity', '0'], inplace=True)
df = df[df.GeometryConverged == True]
df.drop(columns=["Species", 'JobType', 'Charge', 'Multiplicity', '0', "GeometryConverged"], inplace=True)

df_O2 = df[df["Bound"] == "O2"]
df_O2H = df[df["Bound"] == "O2H"]
df_O = df[df["Bound"] == "O"]
df_OH = df[df["Bound"] == "OH"]
df_CN = df[df["Bound"] == "CN"]
df_CO = df[df["Bound"] == "CO"]

for df_i in [df_O2, df_O2H, df_O, df_OH, df]:
    print(df_i.shape)

# add suffix to all columns except ["data_dir", "Funcnum", "Bound_site"] which are preserved for merging
suffix_list = ["_O2", "_O2H", "_O", "_OH", "_CN", "_CO"]
keep_same = ["Catalyst", "Funcnum"]
for i, df_intermediate in enumerate([df_O2, df_O2H, df_O, df_OH]):
    df_intermediate.columns = ['{}{}'.format(c, '' if c in keep_same else suffix_list[i]) for c in df_intermediate.columns]

df_merge = df_O2.merge(df_O2H, on=keep_same, how=merge_method)
print(df_merge.shape)
print(df_merge.columns)
df_merge = df_merge.merge(df_O, on=keep_same, how=merge_method)
print(df_merge.shape)
print(df_merge.columns)
df_merge = df_merge.merge(df_OH, on=keep_same, how=merge_method)
print(df_merge.shape)
#df_merge.drop_duplicates(inplace=True)
print(df_merge.shape)
print(df_merge.columns)
df_merge.to_json("cycleFe_GeoConv_IntMerge.json")


