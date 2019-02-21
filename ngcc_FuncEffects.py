import sys
import pandas as pd
import numpy as np
import re

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns


import CatO2Df

font = {'size':14}
mpl.rc('font',**font)

#df_fname = sys.argv[1]
#df = pd.read_csv(df_fname)

df_m = pd.read_csv("mepyr1_assembled.csv")
df_tr = pd.read_csv("tetrid1_assembled.csv")
df_ty = pd.read_csv("tetry1_assembled.csv")
df_etr = pd.read_csv("order_tetridEnum_assembled.csv")
df_ety = pd.read_csv("order_tetryEnum_assembled.csv")

### Picking which dataframe
#df = df_m
#df = pd.concat([df_tr, df_etr])
#df = pd.concat([df_ty, df_ety])
df = pd.concat([df_ty, df_ety, df_tr, df_etr, df_m])




#Functional groups are the y-axis, and substitution locations are the x-axis
#bottommost row is the average for the functional group
#create this for each of the 3 catalysts, stack them, and have 1 final total average

#would be good to have a plot like this for variance at each site as well


sns.heatmap(data=fg_arr)
