from math import sin
import argparse

parser = argparse.ArgumentParser(description='Generate time series file: data.txt')
parser.add_argument('N', help='data set size')
args = parser.parse_args()

print 'Generating', args.N, 'data points...'
with open("data.txt", "w") as ts:
  for i in range (0, int(args.N)):
    x = float(i)/100.0
    y = sin (x)
    ts.write("%g %g\n" % (x, y))
