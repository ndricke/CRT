import sys
import pandas as pd
import numpy as np

"""
Create pairplots for all of the reaction energies for the catalytic cycle

TODO:
1. Get adjustments for base catalysts to get free energies of reaction
2. Apply adjustments
3. Calculate free energy of reactions
4. Make pairplot
"""

Ht2eV = 27.211399
kC2eV = 1/23.06054195 

df = pd.read_json("~/work/ORRmol/ngcc_func/catdata_bindE_IntMerge.json")
df.rename(columns={"Catalyst_O2":"Catalyst", "Esolv_bare_O":"Esolv_bare"}, inplace=True)

(df["Bound_site"] == 20.0) & (df["Catalyst"] == "tetry")
df.loc[(df["Bound_site"] == 17.0) & (df["Catalyst"] == "tetry"), ["Catalyst"]] = "tetry-17"
df.loc[(df["Bound_site"] == 20.0) & (df["Catalyst"] == "tetry"), ["Catalyst"]] = "tetry-20"

#df_gform = pd.read_csv("~/work/ORRmol/dGform_catalysts/gform_IntermediateColumns.csv")
df_gform = pd.read_csv("~/work/ORRmol/dGform_catalysts/gform_reshape.csv")
df_gform.replace("None", np.nan, inplace=True)
print(df)

# df is organized with E_binding_<intermediate>. These need to be shifted and converted to dGrxn

# for each bare species, calculate a dGrxn_correction that can just be added to the E_binding values
# we can get the dGrxn(bare --> intermediate) = dG(i) - dG(0); to go to dGrxn(i --> j) = dG(j) - dG(i) = dG(j) - dG(0) - dG(i) + dG(0)
# so just take the difference between dGrxn(0 --> i) values, as the dG(0) will each cancel

print(df_gform)

df_gform["dGrxn_corr_O2"] = df_gform["dGrxn_O2"] - df_gform["E_O2"] - df_gform["E"]
df_gform["dGrxn_corr_O2H"] = df_gform["dGrxn_O2H"] - Ht2eV*(df_gform["E_O2H"]-df_gform["E"])
df_gform["dGrxn_corr_O"] = df_gform["dGrxn_O"] - Ht2eV*(df_gform["E_O"] - df_gform["E_O2H"])
df_gform["dGrxn_corr_OH"] = df_gform["dGrxn_OH"] - Ht2eV*(df_gform["E_OH"] - df_gform["E_O"])
df_gform["dGrxn_corr_regen"] = df_gform["dGrxn"] - Ht2eV*(df_gform["E"] - df_gform["E_OH"])

dfm = df.merge(df_gform[["Catalyst", "dGrxn_corr_O2","dGrxn_corr_O2H","dGrxn_corr_O","dGrxn_corr_OH","dGrxn_corr_regen"]], on="Catalyst", how="left")
dfm["dGrxn_O2"] = Ht2eV*dfm["E_binding_O2"] + dfm["dGrxn_corr_O2"]
dfm["dGrxn_O2H"] = Ht2eV*(dfm["E_binding_O2H"]) + dfm["dGrxn_corr_O2H"]
dfm["dGrxn_O"] = Ht2eV*(dfm["E_binding_O"]-dfm["E_binding_O2H"]) + dfm["dGrxn_corr_O"] # - 2.46
dfm["dGrxn_OH"] = Ht2eV*(dfm["E_binding_OH"]-dfm["E_binding_O"]) + dfm["dGrxn_corr_OH"]
dfm["dGrxn_regen"] = Ht2eV*(dfm["Esolv_bare"]-dfm["Esolv_OH"]) + dfm["dGrxn_corr_regen"] # - 2.46

print(dfm[["dGrxn_O2H","dGrxn_O","dGrxn_OH", "dGrxn_regen"]])
print(dfm[["dGrxn_O2H","dGrxn_O","dGrxn_OH", "dGrxn_regen"]].sum(axis=1))
print(dfm.columns)
dfm.to_json("catdata_ngcc_dGrxn.json")
