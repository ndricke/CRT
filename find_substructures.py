import os
from rdkit import Chem
import pandas as pd
from openbabel import pybel


def smiles_from_dir(indir):
    # read all xyz files in a directory, and infer their smiles structures
    smiles_dict = {}
    for infile in os.listdir(indir):
        if infile.split(".")[-1] == "xyz":
            smiles = infer_smiles(indir+'/'+infile)
            smiles_dict[infile] = smiles
    return smiles_dict


def dict_substructs(substruct, smiles_dict):
    """
    check if each smiles in list has substruct
    Input:
    substruct (rdkit mol): substructure to check exists in the smiles in smiles_dict
    smiles_dict (dictionary): filename-smiles 
    Output:
    (dictionary): filename: presence or absence
    """
    match_dict = {}
    for key, value in smiles_dict.items():
        m = Chem.MolFromSmiles(value)
        # check if substruct is a substructure of m
        match_dict[key] = m.HasSubstructMatch(substruct)
    return match_dict


def smiles_substruct(smiles, substruct):
        print("I'm trying ", smiles)
        m = Chem.MolFromSmiles(smiles)
        m_sub = Chem.MolFromSmiles(substruct)
        print("Got the big m")
        return m.HasSubstructMatch(m_sub)


def infer_smiles(xyz_filename):
    mol = next(pybel.readfile("xyz", xyz_filename))
    return mol.write('can').split("\t")[0].strip()


"""
>>> infer_smiles.infer_smiles("tetryO-17_fq_a0m2.xyz")
'[O][C@H]1C=CN2c3c1c1nc4ccccc4nc1c1c3c(C=C2)ccc1'
>>> infer_smiles.infer_smiles("tetryO-20_fq_a0m2.xyz")
'c1cc2C=CN3c4c2c(c1)c1nc2ccccc2nc1c4[CH][C@H]1[C@@H]3O1'
"""

if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--substructure", help="substructure to search for")
    parser.add_argument("-d", "--directory", help="search for substructure in the xyz's in this directory")
    args = parser.parse_args()

    smiles_list = []
    inferred_smiles = smiles_from_dir(args.directory)

    m_sub = Chem.MolFromSmiles(args.substructure)
    substructure_match_dict = dict_substructs(m_sub, inferred_smiles)
    s = pd.Series(substructure_match_dict)
    print(s)




