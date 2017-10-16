from collections import defaultdict
import matplotlib.pyplot as plt
from subprocess import call
from math import log
import os, sys
import ast

if len(sys.argv) < 2:
  print 'SYNOPSIS: python path/to/scaling-plots.py path/to/top/level/output/directory'
  print '          e.g. python scaling-plots.py out/array-of-cubes'
  sys.exit(0)

path = sys.argv[1]
timings_path = os.path.join(path,'TIMINGS')
runtimes_path = os.path.join(path,'RUNTIMES')

# remove TIMINGS if exist
if os.path.isfile(timings_path):
  os.remove(timings_path)

# generate fresh TIMINGS
for d0 in os.listdir(path):
   d1 = os.path.join(path,d0)
   if os.path.isdir(d1) and 'XDMF' not in d1:
     call (['solfec',d1])

# generate TIMINGS_{dataset}.png figures
timings = defaultdict(list)
f = open(timings_path, 'r')
lines = set([line for line in f])
for x in sorted(lines):
  y = x.split (' = ')
  key = y[0][:y[0].rfind('_')]
  num = int(y[0][y[0].rfind('_')+1:][1:])
  val = ast.literal_eval(y[1])
  timings[key].append((num,val))

for key in timings:
  sncpu = []
  lncpu = []
  total = []
  timint = []
  condet = []
  conupd = []
  locdyn = []
  consol = []
  parbal = []
  for pair in timings[key]:
    sncpu.append(str(pair[0]))
    lncpu.append(log(pair[0]))
    total.append(pair[1]['TOTAL'])
    timint.append(pair[1]['TIMINT'])
    condet.append(pair[1]['CONDET'])
    conupd.append(pair[1]['CONUPD'])
    locdyn.append(pair[1]['LOCDYN'])
    consol.append(pair[1]['CONSOL'])
    parbal.append(pair[1]['PARBAL'])

  plt.clf ()
  plt.plot (lncpu, total, label = 'total')
  plt.plot (lncpu, timint, label = 'timint')
  plt.plot (lncpu, condet, label = 'condet')
  plt.plot (lncpu, conupd, label = 'conupd')
  plt.plot (lncpu, locdyn, label = 'locdyn')
  plt.plot (lncpu, consol, label = 'consol')
  plt.plot (lncpu, parbal, label = 'parbal')
  plt.xticks (lncpu, sncpu)
  plt.legend (loc = 'upper right')
  plt.xlabel ('NCPU')
  plt.ylabel ('TIME (s)')
  plt.title (key)
  plt.savefig (timings_path + '_' + key + '.png')

# generate RUNTIMES_{dataset}.png figures
timings = defaultdict(list)
f = open(runtimes_path, 'r')
lines = set([line for line in f])
for x in sorted(lines):
  y = x.split (' = ')
  key = y[0][:y[0].rfind('_')]
  num = int(y[0][y[0].rfind('_')+1:][1:])
  val = float(y[1])
  timings[key].append((num,val))

plt.clf ()

for key in timings:
  sncpu = []
  lncpu = []
  total = []
  for pair in timings[key]:
    sncpu.append(str(pair[0]))
    lncpu.append(log(pair[0]))
    total.append(pair[1])
  plt.plot (lncpu, total, label = key)

plt.xticks (lncpu, sncpu)
plt.legend (loc = 'upper right')
plt.xlabel ('NCPU')
plt.ylabel ('TIME (s)')
plt.title ('Total runtimes')
plt.savefig (runtimes_path + '.png')
