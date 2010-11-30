# double pendulum
from math import sin
from math import cos

PI = 3.14159265358979323846 

GEOMETRIC_EPSILON (1E-5)

def double_pendulum_create (material, solfec):
  x1 =  sin (PI/3.0)
  z1 = -cos (PI/3.0)
  x2 = x1 + sin (PI/5.0)
  z2 = z1 - cos (PI/5.0)

  a = (0, 0, 2.0)
  b = (x1, 0, 2.0 + z1)
  c = (x2, 0, 2.0 + z2)
  J = (1, 0, 0, 0, 1, 0, 0, 0, 1)

  sph = SPHERE (b, 0.05, 1, 1)
  bod1 = BODY (solfec, 'RIGID', sph, material, label = 'BALL1')
  BODY_CHARS (bod1, 1.0, 1.0, b, J)
  PUT_RIGID_LINK (bod1, None, b, a)

  sph = SPHERE (c, 0.05, 1, 1);
  bod2 = BODY (solfec, 'RIGID', sph, material, label = 'BALL2')
  BODY_CHARS (bod2, 1.0, 1.0, c, J)
  PUT_RIGID_LINK (bod2, bod1, c, b)

  nodes = [-0.10, -0.2, 0.0,
           -0.05, -0.2, 0.0,
           -0.05,  0.2, 0.0,
           -0.10,  0.2, 0.0,
           -0.10, -0.2, 2.0,
           -0.05, -0.2, 2.0,
           -0.05,  0.2, 2.0,
           -0.10,  0.2, 2.0]

  surfaces = [0, 0, 0, 0, 0, 0]

  msh = HEX (nodes, 1, 1, 1, 0, surfaces)
  BODY (solfec, 'OBSTACLE', msh, material)

def double_pendulum_run (step, stop, gravity, skip):

  solfec = SOLFEC ('DYNAMIC', step, 'out/tests/double-pendulum')
  if not VIEWER(): solfec.verbose = 'OFF'

  bulkmat = BULK_MATERIAL (solfec)

  SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', restitution = 0.1, friction = 0.0)

  GRAVITY (solfec, (0, 0, -gravity))

  double_pendulum_create (bulkmat, solfec)

  gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'EXIT')

  if not VIEWER():
    if solfec.mode == 'READ' and not skip: print '\nPrevious test results exist. Please "make del" and rerun tests'
    else: RUN (solfec, gs, stop)

    if not skip: return solfec.mode
    else: return solfec
  else:
    RUN (solfec, gs, stop)
    return (gs, solfec) # return solver and solfec to avoid their deallocation

# main module
#import rpdb2; rpdb2.start_embedded_debugger('a')

gravity = 9.81
step = 0.001
stop = 2.5

if VIEWER():
  data = double_pendulum_run (step, stop, gravity, 1) # store return data to avoid memory deallocation
else:
  mode = double_pendulum_run (step, stop, gravity, 0)
  if mode == 'WRITE':

    solfec = double_pendulum_run (step, stop, gravity, 1) # reopen in read mode

    if solfec.mode == 'WRITE': print 'ERROR'
    else:
      bod1 = BYLABEL (solfec, 'BODY', 'BALL1')
      bod2 = BYLABEL (solfec, 'BODY', 'BALL2')
      th = HISTORY (solfec, [(solfec, 'KINETIC'), (solfec, 'INTERNAL'), (solfec, 'EXTERNAL'), (solfec, 'CONTACT'), (solfec, 'FRICTION'),
			     (bod1, bod1.center, 'DX'), (bod2, bod2.center, 'DX'), (bod1, bod1.center, 'DZ'), (bod2, bod2.center, 'DZ')], 0, stop)

      for i in range (0, len (th[6])): th [6][i] += bod1.center [0]
      for i in range (0, len (th[7])): th [7][i] += bod2.center [0]

      tot = []
      for i in range (0, len (th[0])):
	tot.append (th[1][i]
		    + bod1.mass * gravity * max (bod1.center [2] + th[8][i], 0)
		    + bod2.mass * gravity * max (bod2.center [2] + th[9][i], 0))

      shift = tot [len (th[0]) - 1]
      for i in range (0, len (th[0])): tot [i] -= shift

      REF = [(0, 11.6757089834),
	     (600, 9.02830283563),
	     (1000, 0.891186821736),
	     (1500, 0.64127866263),
	     (2400, 5.16299181186e-09)]

      passed = 1

      for ref in REF:
	error = abs (tot[ref[0]] - ref[1]) / max (1, ref[1])
	if error > 1E-3:
	  passed = 0
	  print 'FAILED'
	  print '(', 'Computed total energy was %.3f' % tot[ref[0]], 'while the reference is %.3f' % ref[1]
	  break

      if passed: print 'PASSED'

      try:
	import matplotlib.pyplot as plt
        plt.clf ()
	plt.plot (th [0], th [6], label='$x_1$')
	plt.plot (th [0], th [7], label='$x_2$')
	plt.axis (xmin = 0, xmax = stop, ymin = 0, ymax = 1.5)
	plt.xlabel ('Time [s]')
	plt.ylabel ('Position [m]')
	plt.legend(loc = 'upper right')
	plt.savefig ('out/tests/double-pendulum/double-pendulum-dx.eps')
	plt.clf ()
	plt.plot (th [0], tot, label='total energy')
	plt.axis (xmin = 0, xmax = stop, ymin = 0)
	plt.xlabel ('Time [s]')
	plt.ylabel ('Energy [J]')
	plt.legend(loc = 'upper right')
	plt.savefig ('out/tests/double-pendulum/double-pendulum-ene.eps')
      except ImportError:
	pass # no reaction
