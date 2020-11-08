import sys
import pandas as pd
import numpy as np
import re
import scipy.stats as stat

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import seaborn as sns


save_columns = ["Catalyst", "dGrxn_tO2H", "dGrxn_regen"]

df_ngcc = pd.read_json("~/work/ORRmol/ngcc_func/catdata_dGrxn.json")
df_cycleFe = pd.read_json("~/work/ORRmol/cycleFe_func/cycleFe_catcycle_NorskConv_dGrxn.json")

# Calc reaction for merging steps 1 and 2, in the way done in Norskov papers
df_ngcc["dGrxn_tO2H"] = df_ngcc["dGrxn_O2"] + df_ngcc["dGrxn_O2H"]
df = pd.concat([df_ngcc[save_columns], df_cycleFe[save_columns]])

# Shift free energy by number of electrons yet transferred
df["dGrxn_regen"] *= -1
df["dGrxn_regen"] += -0.27
df["dGrxn_tO2H"] += 1.23*4 - 0.6


