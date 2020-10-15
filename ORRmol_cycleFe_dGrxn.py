import sys
import pandas as pd
import numpy as np

"""
df is organized with E_binding_<intermediate>. These need to be shifted and converted to dGrxn
for each bare species, calculate a dGrxn_correction that can just be added to the E_binding values
we can get the dGrxn(bare --> intermediate) = dG(i) - dG(0); to go to dGrxn(i --> j) = dG(j) - dG(i) = dG(j) - dG(0) - dG(i) + dG(0)
so just take the difference between dGrxn(0 --> i) values, as the dG(0) will each cancel
"""

Ht2eV = 27.211399
kC2eV = 1/23.06054195 
CN_shift = -1.292 + 4.42
CO_shift = -1.422

df = pd.read_json("~/work/ORRmol/cycleFe_func/cycleFe_GeoConv_CNCO_IntMerge.json")
df = df[df["Filename_O2"] != "porphyrinFe-functionalized95O2_optcdftsp_a0m1.out"]  # bare didn't converge
df.rename(columns={"Esolv_bare_O":"Esolv_bare"}, inplace=True)
df.Catalyst = df.Catalyst.str.replace("-functionalized", "")
#df.columns = ['{}{}'.format(c, '' if c in keep_same else suffix_list[i]) for c in df_intermediate.columns]
df.columns = [c.replace("Binding_Energy", "E_binding") for c in df.columns]
#df.dropna(inplace=True)
#df = df[df.GeometryConverged_bare_OH == True]
#df = df[df.GeometryConverged_O2H == True]
#df = df[df.GeometryConverged_O == True]
#df = df[df.GeometryConverged_OH == True]
print(df.columns)


df_gform = pd.read_csv("~/work/ORRmol/dGform_catalysts/cycleFe_gform_CNCO_IntMerge.csv")
df_gform.replace("None", np.nan, inplace=True)

print(df)

df_gform["dGrxn_corr_O2"] = kC2eV*(df_gform["dG_O2"] - df_gform["dG_bare"]) - Ht2eV*(df_gform["E_O2"] - df_gform["E_bare"])
df_gform["dGrxn_corr_O2H"] = kC2eV*(df_gform["dG_O2H"]-df_gform["dG_O2"]) - Ht2eV*(df_gform["E_O2H"]-df_gform["E_O2"])
df_gform["dGrxn_corr_O"] = kC2eV*(df_gform["dG_O"] - df_gform["dG_O2H"]) - Ht2eV*(df_gform["E_O"] - df_gform["E_O2H"])
df_gform["dGrxn_corr_OH"] = kC2eV*(df_gform["dG_OH"] - df_gform["dG_O"]) - Ht2eV*(df_gform["E_OH"] - df_gform["E_O"])
df_gform["dGrxn_corr_regen"] = kC2eV*(df_gform["dG_bare"] - df_gform["dG_OH"]) - Ht2eV*(df_gform["E_bare"] - df_gform["E_OH"])
df_gform["dGrxn_corr_CN"] = kC2eV*(df_gform["dG_CN"] - df_gform["dG_bare"]) - Ht2eV*(df_gform["E_CN"] - df_gform["E_bare"])
df_gform["dGrxn_corr_CO"] = kC2eV*(df_gform["dG_CO"] - df_gform["dG_bare"]) - Ht2eV*(df_gform["E_CO"] - df_gform["E_bare"])
df_gform["dGrxn_corr_tO2H"] = kC2eV*(df_gform["dG_O2H"]-df_gform["dG_bare"]) - Ht2eV*(df_gform["E_O2H"]-df_gform["E_bare"])

print(df_gform)
print(df_gform.shape)

dfm = df.merge(df_gform[["Catalyst", "dGrxn_corr_O2","dGrxn_corr_O2H","dGrxn_corr_O","dGrxn_corr_OH",
    "dGrxn_corr_regen", "dGrxn_corr_CN", "dGrxn_corr_CO", "dGrxn_corr_tO2H"]], 
    on="Catalyst", how="left")
dfm["dGrxn_O2"] = Ht2eV*dfm["E_binding_O2"] + dfm["dGrxn_corr_O2"]
dfm["dGrxn_O2H"] = Ht2eV*(dfm["E_binding_O2H"]-dfm["E_binding_O2"]) + dfm["dGrxn_corr_O2H"]
dfm["dGrxn_O"] = Ht2eV*(dfm["E_binding_O"]-dfm["E_binding_O2H"]) + dfm["dGrxn_corr_O"] - 2.46
dfm["dGrxn_OH"] = Ht2eV*(dfm["E_binding_OH"]-dfm["E_binding_O"]) + dfm["dGrxn_corr_OH"]
dfm["dGrxn_regen"] = Ht2eV*(dfm["Esolv_bare"]-dfm["Esolv_OH"]) + dfm["dGrxn_corr_regen"] - 2.46
dfm["dGrxn_CN"] = Ht2eV*(dfm["E_binding_CN"]) + dfm["dGrxn_corr_CN"] + CN_shift
dfm["dGrxn_CO"] = Ht2eV*(dfm["E_binding_CO"]) + dfm["dGrxn_corr_CO"] + CO_shift
dfm["dGrxn_tO2H"] = Ht2eV*(dfm["E_binding_O2H"]) + dfm["dGrxn_corr_tO2H"]

print(dfm[["dGrxn_O2","dGrxn_O2H","dGrxn_O","dGrxn_OH", "dGrxn_regen", "dGrxn_CN", "dGrxn_CO", "dGrxn_tO2H"]])
print("5-step cycle checksum:")
print(dfm[["dGrxn_O2","dGrxn_O2H","dGrxn_O","dGrxn_OH", "dGrxn_regen"]].sum(axis=1))
print("4-step cycle checksum:")
print(dfm[["dGrxn_tO2H","dGrxn_O","dGrxn_OH", "dGrxn_regen"]].sum(axis=1))
print(dfm.columns)
print(dfm[["dGrxn_corr_O2","dGrxn_corr_O2H","dGrxn_corr_O","dGrxn_corr_OH", "dGrxn_corr_regen", "dGrxn_corr_CN", "dGrxn_corr_CO"]])
dfm.to_json("cycleFe_catcycle_GeoConv_CNCO_dGrxn.json")
