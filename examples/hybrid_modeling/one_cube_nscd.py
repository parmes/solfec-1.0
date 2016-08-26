# two-body hybrid modeling example 
import matplotlib.pyplot as plt
import time, sys, os
sys.path.append(os.path.dirname(__file__))
from acc_sweep import *

step = 1E-5  # time step
stop = 0.05   # duration of the simulation
damp = 1E-5  # amount of stiffness proportional damping
nele = 2     # number of elements per body (along x, y, z)
l = 0.1      # length of one body
w = 0.1      # widhth of one body
h = 0.1      # height of one body
gap = 0.002  # gap
ivel = 0.3   # impact velocity

# geometrical tolerance << gap
GEOMETRIC_EPSILON (1E-9)

# solfec and material
solfec = SOLFEC ('DYNAMIC', step, 'out/hybrid_modeling/one_cube_nscd')
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.05, restitution = 0.0)
bulk = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E9, poisson = 0.25, density = 1E3)

# create body
nodes = [0, 0, 0, w, 0, 0, w, l, 0, 0, l, 0, 0, 0, h, w, 0, h, w, l, h, 0, l, h]
shape = HEX (nodes, nele, nele, nele, 1, [1]*6)
body = BODY (solfec, 'FINITE_ELEMENT', shape, bulk, form = 'BC')
body.scheme = 'DEF_LIM'
body.damping = damp
INITIAL_VELOCITY (body, (0, -ivel, 0), (0, 0, 0))

# create obstacle
shape = HEX (nodes, 1, 1, 1, 1, [1]*6)
TRANSLATE (shape, (0, -gap-l, 0))
obstacle = BODY (solfec, 'OBSTACLE', shape, bulk)
FIX_DIRECTION (obstacle, (w,-(gap+l),0), (1, 0, 0))
FIX_DIRECTION (obstacle, (w,-(gap+l),0), (0, 1, 0))
FIX_DIRECTION (obstacle, (w,-(gap+l),0), (0, 0, 1))
FIX_DIRECTION (obstacle, (w,-gap,0), (1, 0, 0))
FIX_DIRECTION (obstacle, (w,-gap,h), (1, 0, 0))
FIX_DIRECTION (obstacle, (w,-gap,0), (0, 0, 1))

# create constraints solver
slv = NEWTON_SOLVER ()

# run simulation
t0 = time.time()
RUN (solfec, slv, stop)
if solfec.mode == 'WRITE':
  print "Analysis run time:", (time.time() - t0), "seconds"

# post-process results
if not VIEWER() and solfec.mode == 'READ':
  data = []
  data.append ((body, 'KINETIC'))
  data.append ((body, 'INTERNAL'))
  data.append ((body, body.center, 'VY'))
  th = HISTORY (solfec, data, 0, stop)

  plt.plot (th[0], th[1], lw = 2, label = 'kinetic')
  plt.plot (th[0], th[2], lw = 2, label = 'internal')
  plt.legend(loc = 'best')
  plt.xlabel ('Time $(s)$')
  plt.ylabel ('Energy $(J)$')
  plt.savefig ('out/hybrid_modeling/one_cube_nscd/ene.png')

  plt.clf ()
  plt.plot (th[0], th[3], lw = 2)
  plt.xlabel ('Time $(s)$')
  plt.ylabel ('Velocity $(m/s)$')
  plt.savefig ('out/hybrid_modeling/one_cube_nscd/vel.png')

  print 'Velocity restitution:', th[3][len(th[3])-1]/abs(ivel)
