"""Calculate error between ESP for different sample sizes to analysze convergence.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ase
import parmed as pmd
import warnings
import time
from smamp.tools import read_atom_numbers
from compare_potentials import create_structure, get_dimensions, parse_cubefile
from compare_potentials import combine_data, get_esp, parse_charges, extract_dft_esp
from compare_potentials import line_to_xyz, check_distance, reject_sample


def main(sweep=True):
   setup_start = time.time()
   atom_positions, dft_esp, grid_vectors = parse_cubefile(path='esp.cube')
   upper_bound = 7
   lower_bound = 1.8 

   # Calculate HORTON potential
   charge_df = parse_charges()
   pmd2ase, ase2pmd = create_structure()

   for n_samples in [100, 500, 1000, 2500, 5000, 7500, 10000, 20000, 30000, 40000, 50000, 60000]:
      print(n_samples)
      start = time.time()

      # Get the non-rejected probe spots in units of cubefiles line-numbers 
      probe_positions = reject_sample(atom_positions, dft_esp, grid_vectors, 
                                      upper_bound, lower_bound, n_samples)
      # Set up a table 
      df = pd.DataFrame() 
      df['dft_lines'] = probe_positions
      xyz = [line_to_xyz(dft_esp, n, grid_vectors) for n in probe_positions]
      df['x'] = [pos[0] for pos in xyz]
      df['y'] = [pos[1] for pos in xyz]
      df['z'] = [pos[2] for pos in xyz]

      dft_esp_at_probe = []
      for line in probe_positions:
         dft_esp_at_probe += [float(dft_esp[line].strip()) * -1]
      df['dft_esp'] = dft_esp_at_probe

      charge_xyz = combine_data(charge_df, atom_positions, pmd2ase)
      horton_esp_at_probe = []
      for line in probe_positions:
         horton_esp_at_probe += [get_esp(charge_xyz, line_to_xyz(dft_esp, line, grid_vectors))]
      df['horton_esp'] = horton_esp_at_probe

      df['square_dev'] = (df['horton_esp'] - df['dft_esp']).pow(2)
      rrmsd = np.sqrt(df.square_dev.mean())

      with open('convergence.csv', 'a') as outfile:
         correlation = df['horton_esp'].corr(df['dft_esp'])
         run_time = time.time() - start
         outfile.write('{}, {}, {}, {}\n'.format(n_samples, rrmsd, correlation, run_time))

if __name__ == '__main__':
   main()
