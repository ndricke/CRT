import sys
import pandas as pd
import numpy as np
import re
import CatO2Df

indir = "~/work/CRT/autoq/"

df_m = pd.read_csv(indir+"mepyr1_assembled.csv")
df_tr = pd.read_csv(indir+"tetrid1_assembled.csv")
df_ty = pd.read_csv(indir+"tetry1_assembled.csv")
df_etr = pd.read_csv(indir+"order_tetridEnum_assembled.csv")
df_ety = pd.read_csv(indir+"order_tetryEnum_assembled.csv")
df_trAS = pd.read_csv(indir+"order_tetridAS_assembled.csv")
df_tyAS = pd.read_csv(indir+"order_tetryAS_assembled.csv")


### Picking which dataframe

#catalyst = "mepyr"
#catalyst = "tetry"
#catalyst = "tetrid"
catalyst = "all"

df_tr = pd.concat([df_tr, df_etr, df_trAS], sort=True)
df_ty = pd.concat([df_ty, df_ety, df_tyAS], sort=True)
df_ty_17 = df_ty[df_ty["Active_Site"] == 18]
df_ty_20 = df_ty[df_ty["Active_Site"] == 21]


diff_names = ['O2H', 'O', 'OH']
energy_names = ['CatalystOOH_Energy', 'CatalystO_Energy', 'CatalystOH_Energy']


tetrid_shift =(-150.97535,-75.18120,-75.84124)
mepyr_shift =(-150.96147,-75.16639,-75.82740)
ty20_shift = (-150.96852,-75.18826,-75.83892)
ty17_shift = (-150.96756,-75.17159,-75.83235)

df_m = CatO2Df.AddEnergyDiffs(df_m, energy_names, diff_names, mepyr_shift)
#df_tr = CatO2Df.AddEnergyDiffs(df_tr, energy_names, diff_names, tetrid_shift)
#df_ty_17 = CatO2Df.AddEnergyDiffs(df_ty_17, energy_names, diff_names, ty17_shift)
#df_ty_20 = CatO2Df.AddEnergyDiffs(df_ty_20, energy_names, diff_names, ty20_shift)

df_tr = CatO2Df.AddEnergyDiffs(df_tr, energy_names, diff_names, mepyr_shift)
df_ty_17 = CatO2Df.AddEnergyDiffs(df_ty_17, energy_names, diff_names, mepyr_shift)
df_ty_20 = CatO2Df.AddEnergyDiffs(df_ty_20, energy_names, diff_names, mepyr_shift)

#df_m_unshift = CatO2Df.AddEnergyDiffs(df_m, energy_names, diff_names)
#df_m_shift = CatO2Df.AddEnergyDiffs(df_m, energy_names, diff_names, mepyr_shift)
#df_tr_unshift = CatO2Df.AddEnergyDiffs(df_tr, energy_names, diff_names)
#df_ty_17 = CatO2Df.AddEnergyDiffs(df_ty_17, energy_names, diff_names)
#df_ty_20 = CatO2Df.AddEnergyDiffs(df_ty_20, energy_names, diff_names)

#print(df_m_unshift[["Tag", "O2H_diff"]]/27.211)
#print(df_m_unshift[["Tag", "O_diff"]]/27.211)
#print(df_m_unshift[["Tag", "OH_diff"]]/27.211)
#print()
#print(df_m_shift[["Tag", "O2H_diff"]])
#print(df_m_shift[["Tag", "O_diff"]])
#print(df_m_shift[["Tag", "OH_diff"]])

#print(df_tr_unshift[["Tag", "O2H_diff"]])
#print(df_tr_unshift[["Tag", "O_diff"]])
#print(df_tr_unshift[["Tag", "OH_diff"]])
#print()
#print(df_tr_shift[["Tag", "O2H_diff"]])
#print(df_tr_shift[["Tag", "O_diff"]])
#print(df_tr_shift[["Tag", "OH_diff"]])

catalyst_dict = {"mepyr":[df_m], "tetry":[df_ty_17, df_ty_20], "tetrid":[df_tr], "all":[df_m, df_tr, df_ty_17, df_ty_20]}
df = pd.concat(catalyst_dict[catalyst], sort=True)

#df_name_list = ["mepyr_cycle.csv", "tetrid_cycle.csv", "tetry17_cycle.csv", "tetry20_cycle.csv"]
df_name_list = ["mepyr_msh_cycle.csv", "tetrid_msh_cycle.csv", "tetry17_msh_cycle.csv", "tetry20_msh_cycle.csv"]
for i, df in enumerate(catalyst_dict[catalyst]):
    df = df.drop("Unnamed: 0", 1)
    df = df.drop("Unnamed: 0.1", 1)
    df = df[df["Cat-O2_Bond_Length"] < 1.7]
    df = df[df["Cat-O2_Bond_Length"] > 1.3]
    #df = df.drop("Unnamed: 0.1.1", 1)
    #df = df.drop("Unopt_Cat-O2_Energy", 1)
    df.to_csv(df_name_list[i], index=False)

## This dataframe appears to have duplicates in tetrid, and they aren't all the same format
## This may or may not fix the issue of energies being super off from the bare catalyst
#df.to_csv("gcc_assembled.csv", index=False)
