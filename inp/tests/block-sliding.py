# block sliding on a frictional table
from math import cos

T = [] # plots
X = []
V = []

# sliding force
def block_force (bod, q, u, time, step):
  f = (8 * cos (time), 0, 0, 0, 0, 0, 0, 0, 0)
  return f

# plotting callback
def callback_function (bod, solfec):
  T.append (solfec.time)
  X.append (bod.conf [9])
  V.append (bod.velo [3])
  return 1

# main module
base = [-0.5, -0.15, -0.1,
         0.5, -0.15, -0.1,
         0.5,  0.15, -0.1,
        -0.5,  0.15, -0.1,
        -0.5, -0.15,  0.0,
         0.5, -0.15,  0.0,
         0.5,  0.15,  0.0,
        -0.5,  0.15,  0.0]
cube = [-0.15, -0.15, 0.0,
         0.15, -0.15, 0.0,
         0.15,  0.15, 0.0,
        -0.15,  0.15, 0.0,
        -0.15, -0.15, 0.1,
         0.15, -0.15, 0.1,
         0.15,  0.15, 0.1,
        -0.15,  0.15, 0.1]
surfaces = [0, 0, 0, 0, 0, 0]
vector = (3, 0, 0)

step = 0.001
stop = 10.0
gravity = 9.81
solfec = SOLFEC ('DYNAMIC', step, 'out/tests/block-sliding')
if not VIEWER(): solfec.verbose = 'OFF'
material = BULK_MATERIAL (solfec, density = 111.11111111111111111) # so that 0.3*0.3*0.1*111.(1) = 1.0
SURFACE_MATERIAL (solfec, friction = 0.8)
gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000)
GRAVITY (solfec, (0, 0, -gravity))

msh = HEX (base, 1, 1, 1, 0, surfaces)
TRANSLATE (msh, vector)
BODY (solfec, 'OBSTACLE', msh, material)
msh = HEX (cube, 2, 2, 1, 0, surfaces)
TRANSLATE (msh, vector)
bod = BODY (solfec, 'RIGID', msh, material)
FORCE (bod, 'SPATIAL', (0,0,0), (0,0,0), block_force, bod);

if not VIEWER():
  if solfec.mode == 'READ': print '\nPrevious test results exist. Please "make del" and rerun tests'
  else:
    CALLBACK (solfec, step, (bod, solfec), callback_function)
    RUN (solfec, gs, stop)

    REF = [(0, 2000, 3.00435080898),
	   (0, 5000, 2.99129121231),
	   (0, 8000, 3.00441870359),
	   (1, 150, 0.0183699246991),
	   (1, 3150, -0.021223275136),
	   (1, 6150, 0.00288766940501)]

    passed = 1

    for ref in REF:
      if ref [0] == 0: val = X[ref[1]]
      else: error = val = V[ref[1]]
      error = abs (val - ref[2]) / max (1, ref[2])
      if error > 1E-3:
	passed = 0
	print 'FAILED'
	print '(', 'Computed value was %.3f' % val, 'while the reference is %.3f' % ref[2], ')'
	break

    if passed: print 'PASSED'

    try:
      import matplotlib.pyplot as plt
      plt.clf ()
      plt.plot (T, X, label='x')
      plt.axis (xmin = 0, xmax = stop, ymin = 2.99, ymax = 3.006)
      plt.xlabel ('Time [s]')
      plt.ylabel ('Position [m]')
      plt.legend (loc = 'upper right')
      plt.savefig ('out/tests/block-sliding/block-sliding-x.eps')
      plt.clf ()
      plt.plot (T, V, label='$v_x$')
      plt.axis (xmin = 0, xmax = stop, ymin = -0.04, ymax = 0.04)
      plt.xlabel ('Time [s]')
      plt.ylabel ('Velocity [m/s]')
      plt.legend (loc = 'upper right')
      plt.savefig ('out/tests/block-sliding/block-sliding-v.eps')
    except ImportError:
      pass # no reaction

else: RUN (solfec, gs, stop)
