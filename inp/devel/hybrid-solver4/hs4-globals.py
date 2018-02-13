###
import os, sys
def where(program):
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, program)):
      return path
  return None
path = where('parmec4')
if path == None:
  print 'ERROR: parmec4 not found in PATH!'
  print '       Download and compile parmec;'
  print '       Add parmec directory to PATH variable.'
  sys.exit(1)
sys.path.append(os.path.join (path, 'python'))
from acc_sweep import *
from mesh_hex import *
###

dO0 = 0.2
d01 = 0.2
d12 = 0.05
d47 = 0.01
d34 = (d12-d47)/2.
d45 = 0.02
dOz = 0.5
gap = 0.001
step = 1E-3
stop = 5.0
lofq_sweep = 2.0
hifq_sweep = 8.0
lohi_dwell = 3.0
lofq = lohi_dwell
hifq = lohi_dwell
amag = 5.0
