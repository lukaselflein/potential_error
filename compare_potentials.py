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
      e_x = cubefile.readline()
      e_y = cubefile.readline()
      e_z = cubefile.readline()

      # The next 128 lines are the positions of the atoms
      atoms = []
      for i in range(0, nr_atoms):
         line = cubefile.readline()
         q, _, x, y, z = [float(s) for s in line.split() if s.isdigit]
         atoms += [(x, y, z)]
   # esp_start_pos = cubefile.tell()
      print('header read.')
      esp_lines = cubefile.readlines()
      print('readlines() done')
      #dft_esp = np.loadtxt(esp_lines, delimiter=',')
      print(esp_lines[50000])
      exit()

   return atoms, dft_esp# esp_start_pos


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
   charge_position = np.array(charge_position)
   probe_position = np.array(probe_position)
   energy = charge / np.linalg.norm(charge_position - probe_position)
   # TODO: fix units
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


def main():
   print('Parsing cubefile ...')
   atom_positions, esp_start_pos = parse_cubefile(path='esp.cube')
   print('Parsing HORTON charges ...')
   charge_df = parse_charges()
   print('Getting ase - pmd conversion ...')
   pmd2ase, ase2pmd = create_structure()
   print('Calculating HORTON ESP')
 
   extract_dft_esp('esp.cube', (0, 0, 0), esp_start_pos=esp_start_pos)

   charge_xyz = combine_data(charge_df, atom_positions, pmd2ase)
   horton_esp = get_esp(charge_xyz, (0, 0, 0))

   # dft_esp = extract_dft_esp(path='esp.cube', probe_position = (0, 0, 0))

if __name__ == '__main__':
   main()
