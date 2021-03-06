import sys
import numpy as np
import pandas as pd

dir_dict = {"tetr_enum/":"tetrEnum", "tetry1/":"tetry1", "mepyr/":"mepyr", "tetrid1/":"tetrid1", "tetry_AS/":"tetryAs", 
    "tetrid_AS/":"tetridAs"}

merge_method = "inner"
infile = "catdata_all_bindE.json"
df = pd.read_json(infile)
df.drop_duplicates(subset=["Species", "Filename", "data_dir"], inplace=True)
df.drop(columns=['Charge', 'Multiplicity'], inplace=True)
df["fprefix"] = df.Filename.str.split(".").str[0]
print(df.shape)

df_if = pd.read_csv("ngcc_catalyst_init_final_opt.csv")
df_if["fprefix"] = df_if.filename_final.str.split(".").str[0]
df_if["data_dir"] = df_if.parent_dir.replace(dir_dict) 

# drop all rows from df with fprefix and data_dir that match rows in df_break
df = df.merge(df_if[["fprefix", "data_dir", "unchanged", "bridge"]], how="left", on=["fprefix", "data_dir"])
print(df.shape)
df = df[~((df["unchanged"] == False) & (df["bridge"] == False))]
print(df.shape)
df = df[df["GeometryConverged"] == True]
print(df.shape)
print(df["bridge"])

df_O2 = df[df["Bound"] == "O2"]
df_O2H = df[df["Bound"] == "O2H"]
df_O = df[df["Bound"] == "O"]
df_OH = df[df["Bound"] == "OH"]

df_O2.loc[df_O2["data_dir"] == "mepyr", "Bound_site"] = np.NaN  # binding site 16 never actually bound anything, so now it's only 14

for df_i in [df_O2, df_O2H, df_O, df_OH, df]:
    print(df_i.shape)

# add suffix to all columns except ["data_dir", "Funcnum", "Bound_site"] which are preserved for merging
suffix_list = ["_O2", "_O2H", "_O", "_OH"]
keep_same = ["data_dir", "Funcnum", "Bound_site"]
for i, df_intermediate in enumerate([df_O2, df_O2H, df_O, df_OH]):
    df_intermediate.columns = ['{}{}'.format(c, '' if c in keep_same else suffix_list[i]) for c in df_intermediate.columns]

df_merge = df_O2.merge(df_O2H, on=["data_dir", "Funcnum", "Bound_site"], how=merge_method)
print(df_merge.shape)
print(df_merge.columns)
df_merge = df_merge.merge(df_O, on=["data_dir", "Funcnum", "Bound_site"], how=merge_method)
print(df_merge.shape)
print(df_merge.columns)
df_merge = df_merge.merge(df_OH, on=["data_dir", "Funcnum", "Bound_site"], how=merge_method)
print(df_merge.shape)
#df_merge.drop_duplicates(inplace=True)
print(df_merge.shape)
print(df_merge.columns)
df_merge.to_json("catdata_bindE_IntMerge.json")


