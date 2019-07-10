"""Visualize the convergence of potential differences.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load data (wide format)
df = pd.read_csv('convergence.csv')

# Melt data into long format
df = df.melt(id_vars=['n_samples'], value_vars=['rrmsd'])

# Make boxplot of convergence vs samplesize
sns.boxplot(x='n_samples', y='value', data=df)
plt.show()

# Make plot of sample size vs time

