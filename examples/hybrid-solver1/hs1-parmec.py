M = 2 # must be same as hs1-solfec.py
N = 3 # must be same as hs1-solfec.py
gap = 0.002 # must be same as hs1-solfec.py
lofq = 1.0
hifq = 2.0
amag = 0.2
step = 1E-4
stop = 10

import os, sys

def where(program):
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, program)):
      return path
  return None

path = where('parmec4')

if path == None:
  print 'ERROR: parmec4 not found in PATH!'
  print '       Download and compile parmec; add parmec directory to PATH variable;'
  sys.exit(1)

print '(Found parmec4 at:', path + ')'

sys.path.append(os.path.join (path, 'python'))

from acc_sweep import *

(vt, vd, vv, va) = acc_sweep (step, stop, lofq, hifq, amag)

try:
  from scipy.interpolate import interp1d
except:
  print 'ERROR: SciPy interp1d failed to load -->'
  print '       perhaps SciPy needs to be installed'
  sys.exit(1)

vel = interp1d(vt, vv) # linear spline of velocity history ...
def linvel(t): return (vel(t), 0, 0) # ... based on the acceleration sweep function
def angvel(t): return (0, 0, 0) # zero angular velocity signal

matnum = MATERIAL (100, 1E6, 0.25)

def cube (x):
  nodes = [x+0.0, 0.0, 0.0,
	   x+0.1, 0.0, 0.0,
	   x+0.1, 0.1, 0.0,
	   x+0.0, 0.1, 0.0,
	   x+0.0, 0.0, 0.1,
	   x+0.1, 0.0, 0.1,
	   x+0.1, 0.1, 0.1,
	   x+0.0, 0.1, 0.1]
  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]
  colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]
  parnum = MESH (nodes, elements, matnum, colors)
  CONSTRAIN (parnum, [0, 1, 0, 0, 0, 1], [1, 0, 0, 0, 1, 0, 0, 0, 1])
  ANALYTICAL (particle=parnum)

for i in range (0,M): cube (i*(0.1+gap))
for i in range (0,M): cube ((M+N+i)*(0.1+gap))

PRESCRIBE (0, linvel, angvel) # first body
PRESCRIBE (2*M-1, linvel, angvel) # last body

spring_curve = [-1-gap, -1E3, -gap, 0, 1, 0]
damper_curve = [-1, 10, 0, 0, 1, 10]

for i in range (1,M):
  p1 = (i*(0.1+gap)-gap, 0.05, 0.05)
  p2 = (i*(0.1+gap), 0.05, 0.05)
  SPRING (i-1, p1, i, p2, spring_curve, damper_curve, (1, 0, 0))

for i in range (1,M):
  p1 = ((M+N+i)*(0.1+gap)-gap, 0.05, 0.05)
  p2 = ((M+N+i)*(0.1+gap), 0.05, 0.05)
  SPRING (M+i-1, p1, M+i, p2, spring_curve, damper_curve, (1, 0, 0))

print 'PARMEC estimated critical time step:', CRITICAL()

#DEM (5.0, step, 0.01)
