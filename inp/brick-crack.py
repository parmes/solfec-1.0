# prescribed fragmentation with background FE mesh

import sys
sys.path.append ('inp/cores/inc')
from simple_core_base import *

PI = 3.14159265358979323846 

def scene_base (material, solfec):

  vertices = [1, 0, 0,
              1, 1, 0,
              0, 1, 0,
              0, 0, 0,
              1, 0, 1,
              1, 1, 1,
              0, 1, 1,
              0, 0, 1]

  faces = [4, 0, 1, 5, 4, 3,
           4, 1, 2, 6, 5, 3,
           4, 2, 3, 7, 6, 3,
           4, 3, 0, 4, 7, 3,
           4, 0, 3, 2, 1, 3,
           4, 4, 5, 6, 7, 3]

  outd = 0.4598
  margin = 0.05
  thick = 0.1
  lx = outd + (margin + thick)
  ly = outd + (margin + thick)

  cvx = CONVEX (vertices, faces, 3)
  scl = (lx,  ly,  thick)
  vec = (-lx/2, -ly/2, -thick)
  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)

  BODY (solfec, 'OBSTACLE', cvx, material)

# main module
#import rpdb2; rpdb2.start_embedded_debugger('a')
step = 1E-3
stop = 1

solver = GAUSS_SEIDEL_SOLVER (1E-4, 1000)
solfec = SOLFEC ('DYNAMIC', step, 'out/brick-crack/')
bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E7, poisson = 0.2, density = 2E3)
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5, restitution = 0.0)
GRAVITY (solfec, (0, 0, -10))

h = 0.1
scene_base (bulkmat, solfec)
shp = gcore_brick (0, 0, h)
msh = PIPE ((0, 0, h), (0, 0, 0.45), 0.125, 0.13, 3, 8, 2, 0, [0, 0, 0, 0])
bod = BODY (solfec , 'FINITE_ELEMENT', shp, bulkmat, mesh = msh, form = 'BC')
#bod = BODY (solfec , 'FINITE_ELEMENT', msh, bulkmat, form = 'BC')
INITIAL_VELOCITY (bod, (0, 0, -1), (0, 0, 0))

SIMPLIFIED_CRACK (bod, TRANSLATE (bod.center, (0, 0.137, 0)), (1, 0, 0), 3, 'TENSILE', ft=1E3, topoadj = 'OFF')
SIMPLIFIED_CRACK (bod, TRANSLATE (bod.center, (0.137, 0, 0)), (0, 1, 0), 3, 'TENSILE', ft=1E3, topoadj = 'ON') #FIXME => interpenteration (splitting in 3?)

#bod.scheme = 'DEF_LIM'
#bod.damping = 1E-3

RUN (solfec, solver, stop)
