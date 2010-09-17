# FE tetrahedral parts "crushing"
import sys
sys.path.append ('inp/mesh')
from netgenread import *

step = 1E-3
stop = 1
frc = 1E6
ni = 3
nj = 3
nk = 3

solfec = SOLFEC ('DYNAMIC', step, 'out/tetcrush')

solver = GAUSS_SEIDEL_SOLVER (1E-4, 100)

material = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5, restitution = 0.25)

tetmesh = NETGEN_READ ('inp/mesh/part0.mesh', 1, 1)

for i in range (0, ni):
  for j in range (0, nj):
    for k in range (0, nk):
      shape = TRANSLATE (COPY (tetmesh), (4*i, 4*j, 4*k))
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

shp = TRANSLATE (SCALE (COPY (shape), (8*ni, 8*nj, 2.0)), (-2*ni-1.5, -2*nj-1.5, -2.0))
BODY (solfec, 'OBSTACLE', shp, material)

shp = TRANSLATE (SCALE (COPY (shape), (2.0, 4*nj, 4*nk)), (-4.0, -1.5, 0.25))
p1 = shp.vertex (0)
p2 = shp.vertex (1)
p3 = shp.vertex (6)
p4 = shp.vertex (3)
bod = BODY (solfec, 'RIGID', shp, material)
FIX_DIRECTION (bod, p1, (0, 1, 0))
FIX_DIRECTION (bod, p1, (0, 0, 1))
FIX_DIRECTION (bod, p2, (0, 0, 1))
FIX_DIRECTION (bod, p3, (0, 1, 0))
FIX_DIRECTION (bod, p3, (0, 0, 1))
FIX_DIRECTION (bod, p4, (0, 0, 1))
FORCE (bod, 'SPATIAL', MASS_CENTER (bod), (1, 0, 0), frc)

shp = TRANSLATE (SCALE (COPY (shape), (2.0, 4*nj, 4*nk)), (4*ni-1, -1.5, 0.25))
p1 = shp.vertex (5)
p2 = shp.vertex (2)
p3 = shp.vertex (7)
p4 = shp.vertex (4)
bod = BODY (solfec, 'RIGID', shp, material)
FIX_DIRECTION (bod, p1, (0, 1, 0))
FIX_DIRECTION (bod, p1, (0, 0, 1))
FIX_DIRECTION (bod, p2, (0, 0, 1))
FIX_DIRECTION (bod, p3, (0, 1, 0))
FIX_DIRECTION (bod, p3, (0, 0, 1))
FIX_DIRECTION (bod, p4, (0, 0, 1))
FORCE (bod, 'SPATIAL', MASS_CENTER (bod), (-1, 0, 0), frc)

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
