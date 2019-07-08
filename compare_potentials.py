"""Determine the difference between DFT and ESP potential.
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


def create_structure(infile_pdb='snapshot_100.pdb', infile_top='topol.top', hydrogen_file='hydrogen_per_atom.csv', strip_string=':SOL,CL'):
    implicitHbondingPartners = read_atom_numbers(hydrogen_file)

    ase_struct = ase.io.read(infile_pdb)
    pmd_struct = pmd.load_file(infile_pdb)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pmd_top = pmd.gromacs.GromacsTopologyFile(infile_top, parametrize=False)

    # strip water and electrolyte from system (if not yet done in .top)
    pmd_top.strip(strip_string)
    pmd_top.box = pmd_struct.box # Needed because .pdb contains box info
    pmd_top.positions = pmd_struct.positions

    names = [ a.name for a in pmd_top.atoms ]
    residues = [ a.residue.name for a in pmd_top.atoms ]

    ase_index = np.arange(len(ase_struct))

    atom_residue_list = list(zip(names, residues))
    ase2pmd = dict(zip(ase_index, atom_residue_list))
    pmd2ase = dict(zip(atom_residue_list, ase_index))

    return pmd2ase, ase2pmd

def get_dimensions(esp_lines, verbose=True):
      # The number of gridpoints in x/y/z direction is the cuberoot of the total
      grid_dimension_float = np.power(len(esp_lines), 1/3)
      grid_dimension_int = int(round(grid_dimension_float))
      if verbose:
          print('Read grid of dimensions {}, cast to {}.'.format(grid_dimension_float, 
                                                                 grid_dimension_int))
      return grid_dimension_int

def parse_cubefile(path='esp.cube'):
   """Read atom positions from cubfile."""
   with open(path, 'r') as cubefile:
      
      # The first two lines are header
      header = cubefile.readline()
      loop_def = cubefile.readline()
     
      # First entry in third line is the number of atoms
      nr_atoms = cubefile.readline()
      nr_atoms = int([float(s) for s in nr_atoms.split() if s.isdigit][0])

      # Line four to seven are the grid basis vectors
      e_x = float(cubefile.readline().split()[1])
      e_y = float(cubefile.readline().split()[2])
      e_z = float(cubefile.readline().split()[3])
      grid_vectors = (e_x, e_y, e_z)

      # The next 128 lines are the positions of the atoms
      atoms = []
      for i in range(0, nr_atoms):
         line = cubefile.readline()
         q, _, x, y, z = [float(s) for s in line.split() if s.isdigit]
         atoms += [np.array([x, y, z])]
      print('header read.')
      dft_esp = cubefile.readlines()

   return atoms, dft_esp, grid_vectors


def combine_data(df, atom_positions, pmd2ase):
   """Write positional data into charge data frame."""
   def lookup(row, xyz):
      atom_name, residue = row['atom'], row['residue']
      ase_index =  pmd2ase[(atom_name, residue)]
      return atom_positions[ase_index][xyz]
     
   df['x'] = df.apply(lookup, args=(0,), axis=1)
   df['y'] = df.apply(lookup, args=(1,), axis=1)
   df['z'] = df.apply(lookup, args=(2,), axis=1)
   return df


def get_energy(charge, charge_position, probe_position):
   """Calculate electrostatic energy of probe wrt a point charge."""
   # Cast to np.array
   charge_position = np.array(charge_position)
   probe_position = np.array(probe_position)
   
   # Calculate ESP energy
   energy = charge / np.linalg.norm(charge_position - probe_position)

   # Fix units to Hartree
   energy *= ase.units.Hartree / 4 / np.pi 

   return energy


def get_esp(df, probe_position):
   """Add up the energy contribution of point charges.""" 
   esp = 0
   for index, row in df.iterrows():
      esp += get_energy(charge=row.q,
                        charge_position=(row.x, row.y, row.z),
                        probe_position=probe_position)
   return esp


def parse_charges(path='fitted_point_charges.csv'):
   """Read HORTON point charge .csv file."""
   df = pd.read_csv(path)
   df = df[['atom', 'residue', 'q']]
   return df

def extract_dft_esp(path, probe_position, esp_start_pos):
   with open(path, 'r') as cubefile:
      cubefile.seek(esp_start_pos)

def line_to_xyz(dft_esp, line_number, grid_vectors, test=False):
   """Convert cartesian coordinates to line number in cube file."""
   grid_dimension = get_dimensions(dft_esp, verbose=False)

   x = line_number // (grid_dimension ** 2) * grid_vectors[0]
   y = (line_number % (grid_dimension **2)) // grid_dimension * grid_vectors[1]
   z = line_number % grid_dimension * grid_vectors[2]

   if test:
      i = 0
      for x_i in range(0, grid_dimension):
         for y_i in range(0, grid_dimension):
            for z_i in range(0, grid_dimension):
               if i == line_number:
                  xyz = np.array([x_i, y_i, z_i]) * grid_vectors[0]
                  assert all(np.isclose(xyz, np.array([x, y, z])))
                  return np.array([x, y, z])
               i += 1
      raise RuntimeError('Return should have happened in loop.')
 
   return np.array([x, y, z])

def check_distance(xyz, atom_positions, upper_bound, lower_bound):
   distances = []
   for i in range(0, len(atom_positions)):
      distance = np.linalg.norm(xyz - atom_positions[i])
      distances += [distance]
   distances = np.array(distances)
   if any(distances < lower_bound) or all(distances > upper_bound):
      return False
   else:
      return True

def reject_sample(atom_positions, dft_esp, grid_vectors, upper_bound, lower_bound, target_n):
   probe_positions = []
   for i in range(0, target_n * 10000):
      # Sample random positions
      n = int(np.random.uniform(0, len(dft_esp)))
      # convert to cartesian coordinates
      xyz = line_to_xyz(dft_esp, n, grid_vectors)
      # Reject the point if it too close or far from the atom positions
      if check_distance(xyz, atom_positions, upper_bound, lower_bound) is True:
         probe_positions += [n]
         n_positions = len(probe_positions) 
         if n_positions >= target_n:
            break
   if n_positions < target_n:
      raise RuntimeError('Not enough sample points produced: {} of {}'.format(n_positions, target_n))

   return probe_positions
   

def main():
   print('Parsing cubefile ...')
   atom_positions, dft_esp, grid_vectors = parse_cubefile(path='esp.cube')
   upper_bound = 7
   lower_bound = 1.8 

   n_samples = 1000000
   #n_samples = 10

   # Look up DFT eletrostatic potential at probe positions

   # Calculate HORTON potential
   #print('Parsing HORTON charges ...')
   charge_df = parse_charges()
   #print('Getting ase - pmd conversion ...')
   pmd2ase, ase2pmd = create_structure()


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
      dft_esp_at_probe += [float(dft_esp[line].strip())]
   df['dft_esp'] = dft_esp_at_probe

   #print('Calculating HORTON ESP')
   charge_xyz = combine_data(charge_df, atom_positions, pmd2ase)
   horton_esp_at_probe = []
   for line in probe_positions:
      horton_esp_at_probe += [get_esp(charge_xyz, line_to_xyz(dft_esp, line, grid_vectors))]
   df['horton_esp'] = horton_esp_at_probe


   df['square_dev'] = (df['horton_esp'] - df['dft_esp']).pow(2)
   rrmsd = np.sqrt(df.square_dev.mean())

   with open('quality.csv', 'a') as outfile:
      outfile.write('{}, {}\n'.format(rrmsd, df['horton_esp'].corr(df['dft_esp'])))

   print('{} samples done in {}'.format(n_samples, time.time() - start))
   print('Done')

   # TODO: 
   # * make units consistent

if __name__ == '__main__':
   main()
