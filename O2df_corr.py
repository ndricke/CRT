import sys

import pandas as pd
import numpy as np
import sklearn.linear_model as lm

import matplotlib
import matplotlib.pyplot as plt

font = {'size'   : 18}
matplotlib.rc('font', **font)
O2_tpssh = -150.3285416
O2_31gp = -150.32246450 
Ht2eV = 27.211

in_csv = sys.argv[1]

df = pd.DataFrame.from_csv(in_csv)
#df = df.dropna()
df = df[df['Cat-O2_Bond_Length'] > 1.4]
df = df[df['Cat-O2_Bond_Length'] < 1.9]
#df = df[df['Active_Site_ID'] == 'C']

#df['O2_binding_energy'] = (df['CatalystO2_Energy'] - df['Catalyst_Energy'] - O2_tpssh)*Ht2eV
df['O2_binding_energy'] = (df['CatalystO2_Energy'] - df['Catalyst_Energy'] - O2_31gp)*Ht2eV
#df['IE'] = (df['Energy'] - df['Cation_Energy'])*Ht2eV + CHE_elec_E
df['O2_CHELPG'] = df['Oxygen_1_CHELPG']+df['Oxygen_2_CHELPG']
df['IE'] = df['Catalyst_Energy'] - df['Catalyst_c1_Energy']
df['Ch_Diff'] = df['Catalyst_c1_Active_Site_CHELPG'] - df['Catalyst_c1_Active_Site_CHELPG']


df = df[df['O2_binding_energy'] < 0.0]
df = df[df['O2_binding_energy'] > -1.0]
#df = df[df['O2_CHELPG'] > -0.6]

print(df)

plt.scatter(df['O2_CHELPG'], df['O2_binding_energy'], s=80, c='red') #possibly nonlinear, possibly noise
#plt.xlabel(r'O$_2$ Charge')
#plt.ylabel(r'O$_2$ Binding Strength (eV)')

#plt.scatter(df['AS_Bare_CHELPG'], df['O2_binding_energy'], s=80, c='red') #possibly nonlinear, possibly noise
#plt.xlabel(r'Active Site Charge')
#plt.ylabel(r'O$_2$ Binding Strength (eV)')

#plt.scatter(df['Cat-O2_Bond_Length'], df['O2_binding_energy'], s=80, c='red') #possibly nonlinear, possibly noise
#plt.xlabel(r'C-O$_2$ Bond Length')
#plt.ylabel(r'O$_2$ Binding Strength (eV)')

#plt.scatter(df['IE'], df['O2_binding_energy'], s=80, c='red') #possibly nonlinear, possibly noise
#plt.scatter(df['Ch_Diff'], df['O2_binding_energy'], s=80, c='red') #possibly nonlinear, possibly noise


#plt.scatter(df['Cat-O2_Bond_Length'], df['O2_CHELPG'], s=80, c='red') #possibly nonlinear, possibly noise
#plt.scatter(df['AS_Bare_CHELPG'], df['O2_CHELPG'], s=80, c='red') #possibly nonlinear, possibly noise

plt.show()

"""
predictors = ['O2_CHELPG', 'AS_Bare_CHELPG', 'Cat-O2_Bond_Length']
#predictors = ['Redox(eV)', 'CHELPG(a0m3)', 'Diff(a0m3-c1m4)']
#predictors = ['Diff(a0m3-c1m4)']
#predictors = ['Redox(eV)']
#predictors = ['Redox(eV)', 'CHELPG(a0m3)']
#predictors = ['Redox(eV)', 'Diff(a0m3-c1m4)']
#predictors = ['CHELPG(a0m3)', 'Diff(a0m3-c1m4)']

a = 0.001

reg = lm.Ridge(alpha=a)
reg.fit(df[predictors], df['O2_binding_energy'])
y_pred = reg.predict(df[predictors])
print(y_pred)

plt.scatter(df['O2_binding_energy'],y_pred, s=80, c='red')
plt.plot([-1,-0.2],[-1,-0.2])
plt.ylim([-1,-0.2])
plt.xlim([-1,-0.2])
plt.ylabel('Fit (eV)')
plt.xlabel(r'DFT O$_2$ Bonding Energy (eV)')

plt.show()
#plt.savefig('FeO2_Lf.png')

#lassoreg = lm.Lasso(alpha=a, normalize=True, max_iter=1e5)
#lassoreg.fit(df[predictors], df['O2bonding'])
#y_pred = lassoreg.predict(df[predictors])
#print(y_pred)
"""


