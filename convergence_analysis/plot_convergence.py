"""Visualize the convergence of potential differences.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Data Manipulation
print('Loading data (in wide format) ...')
df = pd.read_csv('convergence.csv')

print('Melting data into long format.')
error_df = df.melt(id_vars=['n_samples'], value_vars=['rrmsd'])

# Plotting
sns.set_context("talk", font_scale=0.9)

print('Making boxplot of convergence vs samplesize ...')
fig = plt.figure(figsize=(16,10))
bp = sns.boxplot(x='n_samples', y='value', data=error_df)
bp.figure.savefig('error_vs_samples.png')
plt.clf()

# Make plot of sample size vs time
fig = plt.figure(figsize=(16,10))
time_df = df.melt(id_vars=['n_samples'], value_vars=['time'])
time_df['n_samples'] = pd.to_numeric(time_df.n_samples)
time_df['time [h]'] = time_df['value'] / 60 / 60
pp = sns.pointplot(x='n_samples', y='time [h]', data=time_df, join=False)
pp.figure.savefig('pp_time_vs_samples.png')
plt.clf()

lp = sns.lineplot(x='n_samples', y='time [h]', data=time_df, markers=True, ci=95)
lp.figure.savefig('time_vs_samples.png')

