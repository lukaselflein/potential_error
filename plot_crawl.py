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

# Load data (wide format)
data_paths = find(path='.', folder_keyword='/data', file_keyword='csv', nr_occ=None)

collect_df = pd.DataFrame()
for data_file in data_paths:
   df = pd.read_csv(data_file)
   collect_df = collect_df.append(df)

df = collect_df.reset_index(drop=True)
df.columns = df.columns.map(str.strip)
lp = sns.lineplot(data=df, x='lnrho', y='rrmsd', hue='sigma')


