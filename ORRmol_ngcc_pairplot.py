import sys
import pandas as pd
import numpy as np
import re
import scipy.stats as stat

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns


"""
Create pairplots for all of the reaction energies for the catalytic cycle

TODO:
1. Get adjustments for base catalysts to get free energies of reaction
2. Apply adjustments
3. Calculate free energy of reactions
4. Make pairplot



"""

df = pd.read_csv("catdata_bindE_IntMerge.csv", index_col=0)
df_gform = pd.read_csv("~/work/ORRmol/dGform_catalysts/ngcc_gform.csv", index_col=0)
print(df)


# from df_gform, need the columns mol_tag to know what to map to, then E and dGfaq
# it maybe easier to do this parsing in a separate python script


