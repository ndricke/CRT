import pandas as pd
from CRT import find_substructures
import autoq.scandir_analysis as scana
from rdkit import Chem
from rdkit.Chem import inchi

"""
Check the substructures for a list of directories containing xyz's
1. make sure the main catalyst structure is intact
    a. we have a before and after; we could make sure the structure itself hasn't changed
2. make sure the correct O2 species is bound (and that it hasn't rearranged)
3. for O, categorize based on bridge or non-bridge binding

Iterate over a list of directories, applying slightly different functions to each
"""

def analyze_xyz_dir(file_smiles):
        df_init = pd.Series(file_smiles).to_frame().reset_index()
        df_init.rename(columns={"index": "filename", 0: "SMILES"}, inplace=True)
        df_aug = df_init.filename.apply(scana.map_autoq_catalysts).rename(columns={0:"filename"})
        df_aug = df_aug.merge(df_init, on="filename")
        df_aug["mol"] = df_aug.SMILES.apply(Chem.MolFromSmiles)
        df_aug["inchikey"] = df_aug.mol.apply(inchi.MolToInchiKey)
        return df_aug


xyz_matches = {"tetry1/":{"O/": "cat-O-out/xyz/", "OH/": "cat-OH-out/xyz/", "O2H": "cat-O2H-out/xyz/", "cat-xyz/": "cat-out/xyz/",
    "O2/": "cat-O2-out/xyz/"}}

df_list = []
for parent_dir, subpair in xyz_matches.items():
    for init_dir, final_dir in subpair.items():
        print(init_dir, final_dir)
        init_smiles = find_substructures.smiles_from_dir(parent_dir+init_dir)
        final_smiles = find_substructures.smiles_from_dir(parent_dir+final_dir)

        df_init_aug = analyze_xyz_dir(init_smiles)
        df_final_aug = analyze_xyz_dir(final_smiles)

        df_merge = df_init_aug.merge(df_final_aug[["filename", "Funcnum", "Bound_site", "SMILES", "inchikey"]], 
            on=["Funcnum", "Bound_site"], how="outer", suffixes=("_init","_final"))

        #print(df_merge[["filename_x", "SMILES_x", "filename_y", "SMILES_y"]])
        df_merge["unchanged"] = df_merge["inchikey_init"] == df_merge["inchikey_final"]
        df_merge["bridge"] = df_merge.SMILES_final.apply(find_substructures.smiles_substruct, substruct="C1OC1")
        df_merge["init_dir"] = init_dir
        df_merge["parent_dir"] = parent_dir
        print(df_merge.head(5))
        print(df_merge.shape)
        df_list.append(df_merge)

df_all = pd.concat(df_list)
df_all.to_csv("ngcc_catalyst_init_final_opt.csv")

    
