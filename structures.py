import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter

def xyz_to_cif(filename, a, b, c, alpha, beta, gamma,output):
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
