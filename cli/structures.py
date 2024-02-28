import argparse
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter

def xyz_to_cif(filename, a, b, c, alpha, beta, gamma, output):
  """
  Converts XYZ file format into CIF file format.

  filename: Name of XYZ file (e.g. structure.xyz)
  a, b, c: Lattice parameters of cell
  alpha, beta, gamma: Angles of cell
  output: Name of CIF output (e.g. structure.cif)
  """
  # Convert Cartesian coordinates to fractional coordinates
  alpha_rad = np.deg2rad(alpha)
  beta_rad = np.deg2rad(beta)
  gamma_rad = np.deg2rad(gamma)
    
  # Calculate the lattice vectors
  av = [a, 0, 0]
  bv = [b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0]
  cvx = c * np.cos(beta_rad)
  cvy = (b * c * np.cos(alpha_rad) - bv[0] * cvx) / bv[1]
  cvz = np.sqrt(c**2 - cvx**2 - cvy**2)
  cv = [cvx, cvy, cvz]
    
  cartesian_coords=np.loadtxt(filename, skiprows=2, usecols=(1, 2, 3))
  lattice_matrix = np.array([av, bv, cv])
    
  # Compute the inverse of the lattice matrix
  inv_lattice_matrix = np.linalg.inv(lattice_matrix)
    
  # Multiply the Cartesian coordinates by the inverse lattice matrix
  fractional_coords = np.dot(cartesian_coords, inv_lattice_matrix)
    
  lattice = Lattice.from_parameters(a,b,c,alpha,beta,gamma)
  species = np.loadtxt(filename, dtype=str, skiprows=2, usecols=(0))
    
  struct = Structure(lattice,species,fractional_coords)
    
  CifWriter(struct).write_file(output)

def main():
  parser = argparse.ArgumentParser(description='Convert XYZ file to CIF file')
  parser.add_argument('filename', help='Input XYZ file name')
  parser.add_argument('a', type=float, help='Lattice parameter a')
  parser.add_argument('b', type=float, help='Lattice parameter b')
  parser.add_argument('c', type=float, help='Lattice parameter c')
  parser.add_argument('alpha', type=float, help='Angle alpha')
  parser.add_argument('beta', type=float, help='Angle beta')
  parser.add_argument('gamma', type=float, help='Angle gamma')
  parser.add_argument('output', help='Output CIF file name')
    
  args = parser.parse_args()
    
  xyz_to_cif(args.filename, args.a, args.b, args.c, args.alpha, args.beta, args.gamma, args.output)

if __name__ == '__main__':
  main()
