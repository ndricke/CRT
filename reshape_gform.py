import pandas as pd



#df = pd.read_csv("gform.csv")
df = pd.read_csv("ngcc_gform_tdat-631gp_dGrxn.csv")
df.dropna(inplace=True)

suffix_list = ["_O2", "_O2H", "_O", "_OH", "_CN", "_CO"]

merge_method = "outer"
#df.drop(columns=["H(kcal/mol)", "S(cal/mol.K)", "E (eV)", "H (eV)", "ST (eV)", "G"], inplace=True)
sel_cols = ["Catalyst", "Bound", "E", "dGrxn"]
df = df[sel_cols]

df_bare = df[df["Bound"] == "bare"].drop(columns=["Bound"])
df_O2 = df[df["Bound"] == "O2"].drop(columns=["Bound"])
df_O2H = df[df["Bound"] == "O2H"].drop(columns=["Bound"])
df_O = df[df["Bound"] == "O"].drop(columns=["Bound"])
df_OH = df[df["Bound"] == "OH"].drop(columns=["Bound"])
#df_CN = df[df["Bound"] == "CN"].drop(columns=["Bound"])
#df_CO = df[df["Bound"] == "CO"].drop(columns=["Bound"])

print("Initial shapes:")
for df_i in [df_O2, df_O2H, df_O, df_OH, df]: # df_CN, df_CO]:
    print(df_i.shape)

# add suffix to all columns except ["data_dir", "Funcnum", "Bound_site"] which are preserved for merging
keep_same = ["Catalyst"]
for i, df_intermediate in enumerate([df_O2, df_O2H, df_O, df_OH]):  #, df_CN, df_CO]):
    df_intermediate.columns = ['{}{}'.format(c, '' if c in keep_same else suffix_list[i]) for c in df_intermediate.columns]

df_merge = df_bare.merge(df_O2, on=keep_same, how=merge_method)
df_merge = df_merge.merge(df_O2H, on=keep_same, how=merge_method)
print(df_merge.shape)
print(df_merge.columns)
df_merge = df_merge.merge(df_O, on=keep_same, how=merge_method)
print(df_merge.shape)
print(df_merge.columns)
df_merge = df_merge.merge(df_OH, on=keep_same, how=merge_method)
print(df_merge.shape)
#df_merge = df_merge.merge(df_CN, on=keep_same, how=merge_method)
#df_merge = df_merge.merge(df_CO, on=keep_same, how=merge_method)
#print(df_merge.shape)
print(df_merge.columns)
#df_merge.to_csv("gform_reshape.csv", index=False)
df_merge.to_csv("gform_reshape_631gp.csv", index=False)


