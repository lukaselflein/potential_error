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
import parmed as pmd
import warnings
from smamp.tools import read_atom_numbers
from compare_potentials import create_structure, get_dimensions, parse_cubefile
from compare_potentials import combine_data, get_esp, parse_charges, extract_dft_esp
from compare_potentials import line_to_xyz, check_distance, reject_sample


def calc_error(cube_path='esp.cube', pdb_path='snapshot_600.pdb', top_path='topol.top', 
               hydrogen_path='hydrogen_per_atom.csv',
               upper=7, lower=1.8, n_samples=40000):
   print('Parsing cubefile ...')
   atom_positions, dft_esp, grid_vectors = parse_cubefile(path=cube_path)

   n_samples = 400

   # Calculate HORTON potential
   print('Parsing HORTON charges ...')
   charge_df = parse_charges()
   print('Getting ase - pmd conversion ...')


   pmd2ase, ase2pmd = create_structure(infile_pdb=pdb_path, infile_top=top_path,
                                       hydrogen_file=hydrogen_path)

   print('Rejection sampling of esp locations, n = {} ...'.format(n_samples))
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

   #print('Calculating HORTON ESP ...')
   charge_xyz = combine_data(charge_df, atom_positions, pmd2ase)
   horton_esp_at_probe = []
   for line in probe_positions:
      horton_esp_at_probe += [get_esp(charge_xyz, line_to_xyz(dft_esp, line, grid_vectors))]
   df['horton_esp'] = horton_esp_at_probe

   print('Calculating statistics ...')
   df['square_dev'] = (df['horton_esp'] - df['dft_esp']).pow(2)
   rrmsd = np.sqrt(df.square_dev.mean())
   correlation = df['horton_esp'].corr(df['dft_esp'])

   print('n = {}, rrmsd = {}, corr = {}'.format(n_samples, rrmsd, correlation))
   return rrmsd

def main():
   calc_error()

if __name__ == '__main__':
   main()
