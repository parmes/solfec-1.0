# FE tetrahedral parts "crushing"
import sys
sys.path.append ('inp/mesh')
from netgenread import *

step = 1E-3
stop = 1.5
ni = 10
nj = 10
nk = 10

solfec = SOLFEC ('DYNAMIC', step, 'out/tetcrush')

solver = GAUSS_SEIDEL_SOLVER (1E-2, 100, 1E-5)

material = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.35, restitution = 0.0)

tetmesh = NETGEN_READ ('inp/mesh/part0.mesh', 1, 1)

d = 3.0
for i in range (0, ni):
  for j in range (0, nj):
    for k in range (0, nk):
      shape = TRANSLATE (COPY (tetmesh), (d*i, d*j, 2*k))
      bod = BODY (solfec, 'FINITE_ELEMENT', shape, material, form = 'BC')
      bod.damping = 1E-3
      bod.scheme = 'DEF_LIM'


shape = HULL ([0, 0, 0,
               1, 0, 0,
	       1, 1, 0,
	       0, 1, 0,
	       0, 0, 1,
	       1, 0, 1,
	       1, 1, 1,
	       0, 1, 1], 1, 1)

shp = TRANSLATE (SCALE (COPY (shape), (2.0, d*nj, 2*nk)), (-4.0, -1.0, 0.25))
p1 = shp.vertex (5)
p2 = shp.vertex (7)
bod = BODY (solfec, 'RIGID', shp, material)
FIX_DIRECTION (bod, p1, (0, 0, 1))
FIX_DIRECTION (bod, p2, (0, 0, 1))
frc = TIME_SERIES ([0, 1E6*bod.volume, 0.5, 1E6*bod.volume, 0.5+step, 1E6*bod.volume, stop, 1E8*bod.volume])
FORCE (bod, 'SPATIAL', MASS_CENTER (bod), (1, 0, 0),frc)
tms = TIME_SERIES ([0, 0, 0.5, 0, 0.5+step, 1, stop, 1])
SET_VELOCITY (bod, p2, (0, 1, 0), tms)

shp = TRANSLATE (SCALE (COPY (shape), (2*d*ni, 2*d*nj, 2.0)), (-0.5*d*ni-1.25, -0.5*d*nj-1.0, -2.0))
BODY (solfec, 'OBSTACLE', shp, material)

shp = TRANSLATE (SCALE (COPY (shape), (2.0, d*nj, 2*nk)), (d*(ni-1)+2.0, -1.0, 0))
bod = BODY (solfec, 'OBSTACLE', shp, material)

GRAVITY (solfec, (0, 0, -10))

OUTPUT (solfec, 0.01)

RUN (solfec, solver, stop)

if not VIEWER() and solfec.mode == 'READ':

  timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL', 'PARBAL']
  dur = DURATION (solfec)
  th = HISTORY (solfec, timers, dur[0], dur[1])
  total = 0.0
  for i in range (0, len(timers)):
    sum = 0.0
    for tt in th [i+1]: sum += tt
    print timers [i], 'TIME:', sum
    if i < 6: total += sum
  print 'TOTAL TIME:', total

  try:
    import matplotlib.pyplot as plt
    th = HISTORY (solfec, [(solfec, 'KINETIC'), (solfec, 'INTERNAL'), (solfec, 'EXTERNAL'), (solfec, 'CONTACT'), (solfec, 'FRICTION')], 0, stop)
    plt.plot (th [0], th [1], label='KIN')
    plt.plot (th [0], th [2], label='INT')
    plt.plot (th [0], th [3], label='EXT')
    plt.plot (th [0], th [4], label='CON')
    plt.plot (th [0], th [5], label='FRI')
    tot = []
    #for i in range(0, len (th[0])): tot.append (th[1][i] + th[2][i] - th[3][i])
    #plt.plot (th [0], tot, label='TOT')
    plt.axis (xmin = 0, xmax = stop)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Energy [J]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/tetcrush/tetene.eps')
  except ImportError:
    pass # no reaction
