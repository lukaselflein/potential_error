"""Crawl structure, determine the difference between DFT and ESP potential.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ase
import os
import parmed as pmd
import warnings
import random
from smamp.tools import read_atom_numbers
from smamp.tools import find
from compare_potentials import create_structure, get_dimensions, parse_cubefile
from compare_potentials import combine_data, get_esp, parse_charges, extract_dft_esp
from compare_potentials import line_to_xyz, check_distance, reject_sample


def calc_error(cube_path='esp.cube', pdb_path='snapshot_600.pdb', top_path='topol.top', 
               hydrogen_path='hydrogen_per_atom.csv', charge_path='fitted_point_charges.csv',
               upper=7, lower=1.8, n_samples=40000, pmd2ase=None, out_path=None, time=None,
               lnrho=None, sigma=None):
   atom_positions, dft_esp, grid_vectors = parse_cubefile(path=cube_path)

   n_samples = 400

   # Calculate HORTON potential
   charge_df = parse_charges(charge_path)


   if pmd2ase is None:
      pmd2ase, ase2pmd = create_structure(infile_pdb=pdb_path, infile_top=top_path,
                                          hydrogen_file=hydrogen_path)

   # Get the non-rejected probe spots in units of cubefiles line-numbers
   probe_positions = reject_sample(atom_positions, dft_esp, grid_vectors, 
                                   upper, lower, n_samples)
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
   correlation = df['horton_esp'].corr(df['dft_esp'])
   
   with open(out_path, 'w') as outfile:
      outfile.write('time, lnrho, sigma, rrmsd\n')
      outfile.write('{}, {}, {}, {}\n'.format(time, lnrho, sigma, rrmsd))

   return rrmsd, pmd2ase

def main():
   pmd2ase = None
   root = '../0_charge/1_charge_cycle'
   pdb_path = find(path=root, folder_keyword='0_initial_structure', file_keyword='.pdb')[0]
   top_path = find(path=root, folder_keyword='0_initial_structure', file_keyword='.top')[0]
   hyd_path = find(path='../0_charge', folder_keyword='fitting', file_keyword='hydrogen_per_atom.csv')[0]
   charge_paths = find(path=root, folder_keyword='4_horton_cost_function/lnrho', file_keyword='charges',
                      nr_occ=None)
   lnrho_range = []
   sigma_range = []
   time_range = []
   random.shuffle(charge_paths)

   for charge_path in charge_paths:
      work_dir = os.path.split(os.path.split(os.path.split(charge_path)[0])[0])[0]
      cube_path = find(path=work_dir, folder_keyword='2_dft_calculations', file_keyword='esp.cube')[0]
      # Parse parameters from filename
      lnrho, sigma = charge_path[-15:-4].split('_')[-2:]
      time = charge_path.split('/')[3].split('_')[0]

      if not sigma in ['0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4', '1.6']:
         print('skipping sigma {}'.format(sigma))
         continue

      out_path = 'data/pot_err_{}_{}_{}.csv'.format(time, lnrho, sigma)

      if os.path.exists(out_path):
         print('exists, skipping.')
         continue

      rrmsd, pmd2ase = calc_error(cube_path, pdb_path, top_path, hyd_path, charge_path, 
                                  pmd2ase=pmd2ase, n_samples=4000, out_path=out_path,
                                  time=time, lnrho=lnrho, sigma=sigma)
   print('Done.')
if __name__ == '__main__':
   main()
