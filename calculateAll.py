import argparse, sys, os

parser = argparse.ArgumentParser(description='Calculate grids for ksi analysis for hydrosheds 500 m data')
parser.add_argument('-b','--base_directory', help='base directory containing input files.')
parser.add_argument('prefix', metavar='dir', type=str, 
                   help='prefix of hydrosheds dataset directory')
parser.add_argument('Ao', metavar='Ao', type=float, 
                   help='reference area at x = xo')
parser.add_argument('theta', metavar='theta', type=float, 
                   help='concavity')

args = parser.parse_args()

from demMethods import processAll
if args.base_directory:
    processAll(args.prefix, args.Ao, args.theta, args.base_directory)
else:
    processAll(args.prefix, args.Ao, args.theta)
