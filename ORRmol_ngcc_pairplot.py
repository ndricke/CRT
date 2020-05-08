import sys
import pandas as pd
import numpy as np
import re
import scipy.stats as stat

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns

df = pd.read_csv("catdata_bindE_IntMerge.csv", index_col=0)
print(df)


