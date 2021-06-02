import os
import sys

import numpy as np
import pandas as pd

from autoq import ChemData as CD


def distance(c1, c2):
    return np.linalg.norm(c1 - c2)


def shortest_distance(coords, base_idx, skip=()):
    shortest_dist = sys.float_info.max
    ignore_idxs = [base_idx] + list(skip)
    for i in range(coords.shape[0]):
        if i not in ignore_idxs:
            dist = distance(coords[base_idx,:], coords[i,:])
            if dist < shortest_dist:
                shortest_dist = dist
                bound_idx = i
    return bound_idx, shortest_dist


def locate_O2(atoms, coords):
    # atoms (list) list of atom symbols
    # coords (np array) n_atoms x 3 array of 3d coordinates
    max_O2_bond_length = 1.6  # O-O bond length should be 1.3-1.4. M-O bond lengths can be longer
    
    O_idxs = []  # get index of all O atoms
    for i, sym in enumerate(atoms):
        if sym == "O":
            O_idxs.append(i)
    if len(O_idxs) < 2:  # no O2 if fewer than 2 O atoms
        print("No O2 found")
        return None, None, None, None
            
    O2_found = False
    for ss, O_i1 in enumerate(O_idxs[:-1]):  # iterate up to second-to-last
        for O_i2 in O_idxs[ss+1:]:  # start iteration one past O_i1 to not pick the same
            if np.linalg.norm(coords[O_i1,:] - coords[O_i2,:]) < max_O2_bond_length:
                if O2_found == True:
                    print("Multiple O2 bound")
                    return None, None, None, None
                else:
                    O2_found = True
                    oxy1, oxy2 = O_i1, O_i2 # O_coords index values for oxygens identified as O2

    # find binding site by distance
    closest_O1_idx, dist1 = shortest_distance(coords, oxy1, [oxy2])
    closest_O2_idx, dist2 = shortest_distance(coords, oxy2, [oxy1])
    if dist1 < dist2:
        return oxy1, oxy2, closest_O1_idx, dist1
    else:
        return oxy1, oxy2, closest_O2_idx, dist2


def xyz_locate_O2(xyz_filename):
    atoms, coords = CD.loadXYZ(xyz_filename)
    return locate_O2(atoms, coords)


def dir_locate_O2(xyz_dir):
    sites, bond_lengths, xyz_files = [], [], []
    for xyz_filename in os.listdir(xyz_dir):
        print(xyz_filename)
        if xyz_filename.split(".")[-1] == "xyz":
            oxy1, oxy2, site_idx, bond_length = xyz_locate_O2(os.path.join(xyz_dir, xyz_filename))
            sites.append(site_idx)
            bond_lengths.append(bond_length)
            xyz_files.append(xyz_filename)
    df = pd.DataFrame({"xyz":xyz_files, "site":sites, "bond_length":bond_lengths})
    return df


if __name__ == "__main__":
    indir = sys.argv[1]
    outfile = sys.argv[2]
    df = dir_locate_O2(indir)
    df.to_csv(outfile)
