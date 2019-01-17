import sys

import pandas as pd
import numpy as np
from scipy import stats
import sklearn.linear_model as lm
from sklearn import preprocessing
from sklearn import metrics

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

font = {'size'   : 22}
matplotlib.rc('font', **font)
O2_tpssh = -150.3285416
O2_31gp = -150.32246450 
Ht2eV = 27.211


def AssembleDf(df, active_sites=[], O2_bl = [1.3,1.6], active_site_id='C'):

    df = df.dropna()
#    df = df[df['Cat-O2_Bond_Length'] > O2_bl[0]]
#    df = df[df['Cat-O2_Bond_Length'] < O2_bl[1]]
#    df = df[df['Active_Site_ID'] == active_site_id]

    #df['O2_binding_energy'] = (df['CatalystO2_Energy'] - df['Catalyst_Energy'] - O2_31gp)*Ht2eV
    #df['O2_CHELPG'] = df['Oxygen_1_CHELPG']+df['Oxygen_2_CHELPG']
    #df['IE'] = (df['Catalyst_c1_Energy'] - df['Catalyst_Energy'])*Ht2eV
    #df['Ch_Diff'] = df['Catalyst_Active_Site_CHELPG'] - df['Catalyst_c1_Active_Site_CHELPG']

    df = df.assign(O2_binding_energy = (df['CatalystO2_Energy'] - df['Catalyst_Energy'] - O2_31gp)*Ht2eV)
    df = df.assign(O2_CHELPG = df['Oxygen_1_CHELPG']+df['Oxygen_2_CHELPG'])
    df = df.assign(IE = (df['Catalyst_c1_Energy'] - df['Catalyst_Energy'])*Ht2eV)
    df = df.assign(Ch_Diff = df['Catalyst_Active_Site_CHELPG'] - df['Catalyst_c1_Active_Site_CHELPG'])

    df = df[df['O2_binding_energy'] < 0.5]
    df = df[df['O2_binding_energy'] > -2.0]

#    df = df[df['IE'] > -3.0]
#    df = df[df['IE'] < 3.0]

    #df_list = []
    #print(active_sites)
    #for site in active_sites:
    #    print(df[df['Active_Site'] == site])
    #    df_list.append(df[df['Active_Site'] == site])
    #df = pd.concat(df_list)

    return df

#"""
#in_csv = sys.argv[1]
df1_csv = 'tetrid1_matched.csv'
df2_csv = 'tetrid2_matched.csv'
df3_csv = 'tetry_matched.csv'
df4_csv = 'mepyr_matched.csv'
df_name_list = [df1_csv, df2_csv, df3_csv, df4_csv]
#"""


df_list = []
for df_name in df_name_list:
    df_piece = pd.DataFrame.from_csv(df_name)
    df_piece = AssembleDf(df_piece)
    df_list.append(df_piece)


"""
df = pd.concat(df_list) #in case we want to include data points we would want to combine here

df.to_csv('order_func.csv')
"""


#df = pd.read_csv('order_func.csv')

#slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(df1['IE'], df1['O2_binding_energy'])
#slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(df2['IE'], df2['O2_binding_energy'])
#print(slope1, intercept1, r_value1, p_value1, std_err1)
#print(slope2, intercept2, r_value2, p_value2, std_err2)

"""
fig, ax = plt.subplots()

#plt.scatter(df1['Catalyst_Active_Site_CHELPG'], df1['O2_binding_energy'], s=80, c='blue', label='Cation Single Point') 
#plt.scatter(df2['Catalyst_Active_Site_CHELPG'], df2['O2_binding_energy'], s=80, c='red', label='Optimized Cation') 
#plt.xlabel(r'Active Site CHELPG')

#plt.scatter(df1['Catalyst_c1_Active_Site_CHELPG'], df1['O2_binding_energy'], s=80, c='blue', label='Cation Single Point') 
#plt.scatter(df2['Catalyst_c1_Active_Site_CHELPG'], df2['O2_binding_energy'], s=80, c='red', label='Optimized Cation') 
#plt.xlabel(r'Active Site CHELPG')

x_term = 'O2_CHELPG'
plt.scatter(df_list[0][x_term], df_list[0]['O2_binding_energy'], s=80, c='blue') 
plt.scatter(df_list[1][x_term], df_list[1]['O2_binding_energy'], s=80, c='red') 
plt.scatter(df_list[2][x_term], df_list[2]['O2_binding_energy'], s=80, c='green') 
plt.scatter(df_list[3][x_term], df_list[3]['O2_binding_energy'], s=80, c='orange') 
plt.xlabel(r'Active Site CHELPG Difference')

#plt.scatter(df['O2_CHELPG'], df['O2_binding_energy'], s=80, c='blue') 
#plt.xlabel(r'O$_2$ CHELPG')

#plt.scatter(df1['Cat-O2_Bond_Length'], df1['O2_binding_energy'], s=80, c='blue') 
#plt.xlabel(r'O$_2$ Bond Length (Angstrom)')

#plt.scatter(df['IE'], df['O2_binding_energy'], s=80, c='blue') #possibly nonlinear, possibly noise
#plt.xlabel(r'Ionization Energy (eV)')

#plt.scatter(df['Ch_Diff'], df['O2_binding_energy'], s=80, c='red') #possibly nonlinear, possibly noise

"""

"""
x_term = 'O2_CHELPG'
y_term = 'O2_binding_energy'
cat_list = ['mepyr', 'tetrid1', 'tetrid2', 'tetry']
color_list = ['red', 'blue', 'green', 'orange']

# the filenames are structured like "tetrid-func12..."
# need to select all the rows that start with tetrid, then all the rows that start with tetry, ect.

for i, cat in enumerate(cat_list):
    df_cat = df.loc[df['Catalyst_File_Name'].isin(]
    plt.scatter(df_cat[x_term], df_cat[y_term], s=80, c=color_list[i], label=cat)
plt.xlabel(x_term)
"""

"""
#plt.scatter(df['Cat-O2_Bond_Length'], df['O2_CHELPG'], s=80, c='red') #possibly nonlinear, possibly noise
#plt.scatter(df['AS_Bare_CHELPG'], df['O2_CHELPG'], s=80, c='red') #possibly nonlinear, possibly noise

plt.ylabel(r'O$_2$ Binding Energy (eV)')
#plt.legend(loc=3)
fig.set_size_inches(8.,8.,forward=True)
plt.show()
#plt.savefig('OrderFuncs_tetrid2_O2BLVsO2.png', transparent=True, bbox_inches='tight', pad_inches=0.02)
"""



#"""

df = df_list[3]
#df = pd.concat(df_list[:2])

#df = df.loc[df['Catalyst_File_Name']]
cat_type_list = []
tag_list = []
for name in df['Catalyst_File_Name']:
    name_split = name.split('-')
    cat_type_list.append(name_split[0])
    tag_list.append(int(name_split[1].split('_')[0][4:]))


df = df.assign(Tag = tag_list)
df = df.assign(Catalyst_Type = cat_type_list)

#print(df.loc((df['Tag'] >= 122) & (df['Tag'] <= 242)))
#print(df[(df['Tag'] >= 122) & (df['Tag'] <= 242) & (df['Active_Site'] == 27) ])
#df = df[(df['Tag'] >= 122) & (df['Tag'] <= 242) & (df['Active_Site'] == 27) ]
#df = df[(df['Tag'] >= 122) & (df['Tag'] <= 242) & (df['Active_Site'] == 19) ]
#df = df[(df['Tag'] >= 122) & (df['Tag'] <= 242) & (df['Active_Site'] == (19 or 27)) ]
#df = pd.concat([df, df_list[0], df_list[1]])
df = df[df['Tag'] >= 29]

#sys.exit(-1)
#print(df)
#predictors = ['Ch_Diff', 'IE', 'Catalyst_Active_Site_CHELPG']
#predictors = ['IE', 'Catalyst_Active_Site_CHELPG']
predictors = ['Ch_Diff', 'IE', ]
#predictors = ['IE' ]


a = 0.001

#reg = lm.Ridge(alpha=a)
reg = lm.LinearRegression()
reg.fit(df[predictors], df['O2_binding_energy'])
y_pred = reg.predict(df[predictors])
print(y_pred)

print(metrics.r2_score(df['O2_binding_energy'], y_pred))

#df_pred = preprocessing.scale(df[predictors])
#df_O2 = preprocessing.scale(df['O2_binding_energy'])
#reg = lm.LinearRegression()
#reg.fit(df_pred, df_O2)
#y_pred = reg.predict(df_pred)
#print(y_pred)

fig, ax = plt.subplots()

plt.scatter(df['O2_binding_energy'], y_pred, s=80, c='red')
plot_range = [-0.7, 0.0]
plt.plot(plot_range,plot_range, linewidth=3)
plt.ylim(plot_range)
plt.xlim(plot_range)
plt.ylabel('Fit (eV)')
plt.xlabel(r'DFT O$_2$ Bonding Energy (eV)')

#plt.show()
tick_spacing = 0.05
#ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
fig.set_size_inches(8.,8.,forward=True)
#plt.savefig('OrderFuncs_tetrid-tetry_FitO2-IE-ChDiff.png', transparent=True, bbox_inches='tight', pad_inches=0.02)
plt.show()

#lassoreg = lm.Lasso(alpha=a, normalize=True, max_iter=1e5)
#lassoreg.fit(df[predictors], df['O2bonding'])
#y_pred = lassoreg.predict(df[predictors])
#print(y_pred)
#"""


