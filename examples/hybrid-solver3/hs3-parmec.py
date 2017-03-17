M = 5 # must be same as in hs2-solfec.py
N = 3 # must be same as in hs2-solfec.py
gap = 0.001 # must be same as in hs2-solfec.py
lofq = 1
hifq = 1
amag = 1
step = 1E-4
stop = 5 # must be >= stop in hs2-solfec.py

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
tsv = [None]*(len(vt)+len(vd))
tsv[::2] = vt
tsv[1::2] = vv
tsv = TSERIES (tsv)
ts0 = TSERIES (0.0)
linvel = (tsv, tsv, tsv)
angvel = (ts0, ts0, ts0)

matnum = MATERIAL (100, 1E6, 0.25)

def cube (x, y, z):
  nodes = [x+0.0, y+0.0, z+0.0,
	   x+0.1, y+0.0, z+0.0,
	   x+0.1, y+0.1, z+0.0,
	   x+0.0, y+0.1, z+0.0,
	   x+0.0, y+0.0, z+0.1,
	   x+0.1, y+0.0, z+0.1,
	   x+0.1, y+0.1, z+0.1,
	   x+0.0, y+0.1, z+0.1]
  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]
  colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]
  parnum = MESH (nodes, elements, matnum, colors)
  CONSTRAIN (parnum, angular=[1, 0, 0, 0, 1, 0, 0, 0, 1])
  ANALYTICAL (particle=parnum)
  return parnum

ijkmap = {}
for i in range (0,M+N+M):
  for j in range (0,M+N+M):
    for k in range (0,M+N+M):
      if i >= M and j >= M and i < M+N and j < M+N and k >= M and k < M+N: continue
      else:
	num = cube (i*(0.1+gap), j*(0.1+gap), k*(0.1+gap))
	ijkmap[(i,j,k)] = num

for (i,j,k) in ijkmap:
  outer = [0, M+N+M-1]
  if i in outer or j in outer or k in outer:
    num = ijkmap[(i,j,k)]
    PRESCRIBE (num, linvel, angvel) # outer most shell of bodies

spring_curve = [-1-gap, -1E3, -gap, 0, 1, 0]
damper_curve = [-1, -7, 1, 7]

ijkmax = M+N+M-1
for (i,j,k) in ijkmap:
  if i < ijkmax and not (i == M-1 and j in range(M,M+N) and k in range(M,M+N)):
    p1 = (i*(0.1+gap)+0.1, j*(0.1+gap)+0.05, k*(0.1+gap)+0.05)
    p2 = (i*(0.1+gap)+0.1+gap, j*(0.1+gap)+0.05, k*(0.1+gap)+0.05)
    n1 = ijkmap[(i,j,k)]
    n2 = ijkmap[(i+1,j,k)]
    SPRING (n1, p1, n2, p2, spring_curve, damper_curve, (1, 0, 0))
  if j < ijkmax and not (j == M-1 and i in range(M,M+N) and k in range(M,M+N)):
    p1 = (i*(0.1+gap)+0.05, j*(0.1+gap)+0.1, k*(0.1+gap)+0.05)
    p2 = (i*(0.1+gap)+0.05, j*(0.1+gap)+0.1+gap, k*(0.1+gap)+0.05)
    n1 = ijkmap[(i,j,k)]
    n2 = ijkmap[(i,j+1,k)]
    SPRING (n1, p1, n2, p2, spring_curve, damper_curve, (0, 1, 0))
  if k < ijkmax and not (k == M-1 and i in range(M,M+N) and j in range(M,M+N)):
    p1 = (i*(0.1+gap)+0.05, j*(0.1+gap)+0.05, k*(0.1+gap)+0.1)
    p2 = (i*(0.1+gap)+0.05, j*(0.1+gap)+0.05, k*(0.1+gap)+0.1+gap)
    n1 = ijkmap[(i,j,k)]
    n2 = ijkmap[(i,j,k+1)]
    SPRING (n1, p1, n2, p2, spring_curve, damper_curve, (0, 0, 1))

print 'PARMEC estimated critical time step:', CRITICAL()

#DEM (step, step, 0.01)
