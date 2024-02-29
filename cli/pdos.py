import argparse
from cp2k_utils.pdos import pdos_plot


def main():
  parser = argparse.ArgumentParser(description='Automatic plotting of DOS')
  parser.add_argument('elements', nargs='+', help='list of elements (in correct order)')
  parser.add_argument('-p', '--spin', type=bool, help='whether calculation is spin polarised')
  parser.add_argument('-s', '--sigma', type=float, help='sigma for smearing')
    
  args = parser.parse_args()
    
  pdos_plot(args.elements, args.spin, args.sigma)

if __name__ == '__main__':
  main()
