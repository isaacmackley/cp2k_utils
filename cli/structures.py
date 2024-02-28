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
  # Your existing function code here...

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
