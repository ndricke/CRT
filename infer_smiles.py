import pandas as pd
import pybel


def smiles_from_dir(indir):
    smiles_dict = {}
    for infile in os.listdir(indir):
        if infile.split(".")[-1] == "xyz":
            smiles = infer_smiles(indir+'/'+infile)
            smiles_dict[infile] = smiles
    return smiles_dict


def check_substructs(substruct, smiles_dict):
    # check if each smiles in list has substruct


def infer_smiles(xyz_filename):
    mol = next(pybel.readfile("xyz", xyz_filename))
    return mol.write('can').split("\t")[0].strip()



