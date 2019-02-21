import sys
import pandas as pd
import numpy as np
import re

import CatO2Df

"""
Scripts for munging tetry, tetrid, tetr_enum, and mepyr dataframes
"""

"""
df = pd.read_csv('order_tetrEnum.csv')
df = CatO2Df.AssembleDf(df)

df_tetrid = df[df['Catalyst_Type'] == 'tetrid']
df_tetrid = CatO2Df.AssignSubsts(df_tetrid, "/home/nricke/work/autoq/order_funcs/catfunc_list/tetrid_enum.txt")
df_tetrid.to_csv("order_tetridEnum_assembled.csv")

df_tetry = df[df['Catalyst_Type'] == 'tetry']
df_tetry = CatO2Df.AssignSubsts(df_tetry, "/home/nricke/work/autoq/order_funcs/catfunc_list/tetry_enum.txt")
df_tetry.to_csv("order_tetryEnum_assembled.csv")
"""

"""
CatO2Df.SaveAssembledSubsts('/home/nricke/work/CRT/autoq/order_tetrid1.csv', \
                    '/home/nricke/work/autoq/order_funcs/catfunc_list/tetrid_1.txt', \
                    'tetrid1_assembled.csv')

CatO2Df.SaveAssembledSubsts('/home/nricke/work/CRT/autoq/order_tetry1.csv', \
                    '/home/nricke/work/autoq/order_funcs/catfunc_list/tetry_1.txt', \
                    'tetry1_assembled.csv')
"""

"""
df_mepyr = pd.read_csv('order_mepyr.csv')
df_mepyr = CatO2Df.AssembleDf(df_mepyr)

print(df_mepyr['Tag'])
print(df_mepyr.shape)

# drop mepyr funcs 9,10,12, 23,24,26, 37,38,40
#print(df_mepyr['Tag'] in [9,10,12,23,24,26,37,38,40])

tag_list = [9,10,12,23,24,26,37,38,40]
index_list = []
for row in df_mepyr.itertuples(index=True, name='Pandas'):
    if row.Tag in tag_list:
        index_list.append(row[0])

    # de-increment all the other values so they number 1-33
    if row.Tag > 9:
        df_mepyr['Tag'][row[0]] -= 2
    if row.Tag > 12:
        df_mepyr['Tag'][row[0]] -= 1
    if row.Tag > 23:
        df_mepyr['Tag'][row[0]] -= 2
    if row.Tag > 26:
        df_mepyr['Tag'][row[0]] -= 1
    if row.Tag > 37:
        df_mepyr['Tag'][row[0]] -= 2
    if row.Tag > 40:
        df_mepyr['Tag'][row[0]] -= 1
df_mepyr = df_mepyr.drop(index_list)
df_mepyr = df_mepyr.sort_values(by=["Tag"])

print(df_mepyr['Tag'])
print(df_mepyr.shape)

df_mepyr = CatO2Df.AssignSubsts(df_mepyr, '/home/nricke/work/autoq/order_funcs/catfunc_list/mepyr_1.txt')
df_mepyr.to_csv("mepyr1_assembled.csv")
"""
