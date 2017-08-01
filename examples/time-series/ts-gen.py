#generate unique sin-based time series
from datetime import datetime
from math import sin
from math import pi
import argparse
import random

parser = argparse.ArgumentParser(description='Generate sin-based time series with t=[0 .. 1.0], v=[0 .. 1.0]: ts-data-{1 .. 10}.txt')
parser.add_argument('N', help='data set size')
args = parser.parse_args()
fN = float(args.N)
random.seed(datetime.now())
print 'Generating 10 data sets with', args.N, 'data points...'
for j in range (1, 11):
  name = 'ts-data-%d.txt' % j
  ampl = random.uniform(0, 1)
  nper = random.uniform (1, 5)
  print 'File: %s' % name
  print 'Amplitude: %.2f' % ampl
  print 'Number of periods: %.2f' % nper
  with open(name, "w") as ts:
    for i in range (0, int(args.N)):
      t = float(i)/fN
      x = t*nper*2.0*pi
      y = ampl * sin (x)
      ts.write("%g %g\n" % (t, y))
