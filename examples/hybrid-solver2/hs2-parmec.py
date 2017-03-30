M = 5 # must be same as in hs2-solfec.py
N = 3 # must be same as in hs2-solfec.py
gap = 0.001 # must be same as in hs2-solfec.py
lofq = 1 # low excitation frequency
hifq = 1 # high excitation freqency
amag = 1 # acceleration magnitude
step = 1E-4 # time step
stop = 5 # must be >= stop in hs2-solfec.py

# find path to parmec source directory in order
# to load the acceleration sweep signal script
import os, sys
def where(program):
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, program)):
      return path
  return None
path = where('parmec4')
if path == None:
  print 'ERROR: parmec4 not found in PATH!'
  print '       Download and compile parmec;',
  print 'add parmec directory to PATH variable;'
  sys.exit(1)
print '(Found parmec4 at:', path + ')'
sys.path.append(os.path.join (path, 'python'))

# generate acceleration sweep signal;
# note that parmec will use the associated
# velocity signal, rather than the acceleration itself
from acc_sweep import *
(vt, vd, vv, va) = acc_sweep (step, stop, lofq, hifq, amag)
tsv = [None]*(len(vt)+len(vd))
tsv[::2] = vt
tsv[1::2] = vv
tsv = TSERIES (tsv) # velocity time series
ts0 = TSERIES (0.0) # constant, zero "time series"
linvel = (tsv, tsv, ts0) # linear motion excitation vector
angvel = (ts0, ts0, ts0) # angular motion excitation vector

# create material
matnum = MATERIAL (100, 1E6, 0.25)

# cube creation function
def cube (x, y):
  nodes = [x+0.0, y+0.0, 0.0,
	   x+0.1, y+0.0, 0.0,
	   x+0.1, y+0.1, 0.0,
	   x+0.0, y+0.1, 0.0,
	   x+0.0, y+0.0, 0.1,
	   x+0.1, y+0.0, 0.1,
	   x+0.1, y+0.1, 0.1,
	   x+0.0, y+0.1, 0.1]
  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]
  colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]
  parnum = MESH (nodes, elements, matnum, colors)
  RESTRAIN (parnum, [0, 0, 1], [1, 0, 0, 0, 1, 0, 0, 0, 1])
  ANALYTICAL (particle=parnum)
  return parnum

# generate cube pattern
ijmap = {}
for i in range (0,M+N+M):
  for j in range (0,M+N+M):
    # skip the inner NxN Solfec block
    if i >= M and j >= M and i < M+N and j < M+N: continue
    else: # create the outer MxM, NxM, MxN blocks
      num = cube (i*(0.1+gap), j*(0.1+gap))
      ijmap[(i,j)] = num # map body numbers to (i,j)-grid

# prescribe sine dwell motion 
# of the outer-most shell of bodies
for (i,j) in ijmap:
  outer = [0, M+N+M-1]
  if i in outer or j in outer:
    num = ijmap[(i,j)]
    PRESCRIBE (num, linvel, angvel)

# spring-damper curves
spring_curve = [-1-gap, -1E3, -gap, 0, 1, 0]
damper_curve = [-1, -7, 1, 7]

# insert contact springs
ijmax = M+N+M-1
for (i, j) in ijmap:
  if i < ijmax and not (i == M-1 and j in range(M,M+N)):
    p1 = (i*(0.1+gap)+0.1, j*(0.1+gap)+0.05, 0.05)
    p2 = (i*(0.1+gap)+0.1+gap, j*(0.1+gap)+0.05, 0.05)
    n1 = ijmap[(i,j)]
    n2 = ijmap[(i+1,j)]
    SPRING (n1, p1, n2, p2, spring_curve,
            damper_curve, (1, 0, 0))
  if j < ijmax and not (j == M-1 and i in range(M,M+N)):
    p1 = (i*(0.1+gap)+0.05, j*(0.1+gap)+0.1, 0.05)
    p2 = (i*(0.1+gap)+0.05, j*(0.1+gap)+0.1+gap, 0.05)
    n1 = ijmap[(i,j)]
    n2 = ijmap[(i,j+1)]
    SPRING (n1, p1, n2, p2, spring_curve,
            damper_curve, (0, 1, 0))

# FYI, print out critical time step information 
print 'PARMEC estimated critical time step:', CRITICAL()

#DEM (stop, step, 0.01)
