# --------------------- #
# binary impact example #
# --------------------- #
import matplotlib.pyplot as plt
from math import sin, cos, pi

step = 1E-4  # time step
stop = 0.01  # duration of the simulation
damp = 1E-4  # amount of stiffness proportional damping
ivel = 0.2   # inital velocity
nele = 10    # number of elements per body (along y)
l = 0.1      # length of one body
w = 0.1      # widhth of one body
h = 0.1      # height of one body
gap = 0.001  # gap

GEOMETRIC_EPSILON (1E-9) # tiny geometrical tolerance (<< gap)

solfec = SOLFEC ('DYNAMIC', step, 'out/binimp')
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = 0.0)
bulk = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E9, poisson = 0.25, density = 1E3)

nodes = [0, 0, 0,
	 w, 0, 0,
	 w, l, 0,
	 0, l, 0,
	 0, 0, h,
	 w, 0, h,
	 w, l, h,
	 0, l, h]

shape = HEX (nodes, 1, nele, 1, 1, [1]*6)
bod1 = BODY (solfec, 'FINITE_ELEMENT', shape, bulk, form = 'BC')
INITIAL_VELOCITY (bod1, (0, ivel, 0), (0, 0, 0))
bod1.scheme = 'DEF_LIM' # implicit time integration
bod1.damping = damp

shape = HEX (nodes, 1, nele, 1, 1, [1]*6)
TRANSLATE (shape, (0, l+gap, 0))
bod2 = BODY (solfec, 'FINITE_ELEMENT', shape, bulk, form = 'BC')
INITIAL_VELOCITY (bod2, (0, -ivel, 0), (0, 0, 0))
bod2.scheme = 'DEF_LIM' # implicit time integration
bod2.damping = damp

# create constraints solver
slv = NEWTON_SOLVER ()

# run simulation
RUN (solfec, slv, stop)

# post-process results
if not VIEWER() and solfec.mode == 'READ':
  th = HISTORY (solfec, [(solfec, 'KINETIC'), (solfec, 'INTERNAL'),
       (bod1, bod1.center, 'VY'), (bod2, bod2.center, 'VY')], 0, stop)
  plt.clf ()
  plt.plot (th[0], th[1], lw = 2, label = 'kinetic')
  plt.plot (th[0], th[2], lw = 2, label = 'internal')
  plt.legend(loc = 'best')
  plt.xlabel ('Time $(s)$')
  plt.ylabel ('Energy $(J)$')
  plt.savefig ('out/binimp/ene.png')

  plt.clf ()
  plt.plot (th[0], th[3], lw = 2, label = 'vy1')
  plt.plot (th[0], th[4], lw = 2, label = 'vy2')
  plt.legend(loc = 'best')
  plt.xlabel ('Time $(s)$')
  plt.ylabel ('Velocity $(m/s)$')
  plt.savefig ('out/binimp/vy.png')

  ovel = 0.0
  for i in range (len(th[3])-10, len(th[3])):
    ovel += th[3][i]
  ovel /= 10.0;

  print 'Coefficient of restitution is', -ovel/ivel
