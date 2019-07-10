"""Visualize the convergence of potential differences.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from smamp.tools import find

print('Load data (wide format) ...')
data_paths = find(path='.', folder_keyword='/data', file_keyword='csv', nr_occ=None)

collect_df = pd.DataFrame()
for data_file in data_paths:
   df = pd.read_csv(data_file)
   collect_df = collect_df.append(df)

print('Manipulating data ...')
df = collect_df.reset_index(drop=True)
df.columns = df.columns.map(str.strip)
df['lnrho'] = pd.to_numeric(df['lnrho'])
df['sigma'] = df['sigma'].apply(lambda x: str(x))
df = df.loc[df.lnrho < -3]
df_0 = df.loc[df.charge == 0]
df_1 = df.loc[df.charge == 1]
df_2 = df.loc[df.charge == 2]
print(df)

print('Plotting ...')
p = sns.pointplot(data=df_0, x='lnrho', y='rrmsd', hue='sigma', palette='GnBu_d')
p.figure.savefig('0_lnrho_vs_rrmsd.png')
plt.clf()
p = sns.pointplot(data=df_1, x='lnrho', y='rrmsd', hue='sigma', palette='GnBu_d')
p.figure.savefig('1_lnrho_vs_rrmsd.png')
plt.clf()
p = sns.pointplot(data=df_2, x='lnrho', y='rrmsd', hue='sigma', palette='GnBu_d')
p.figure.savefig('2_lnrho_vs_rrmsd.png')
print('Done.')
