# projectile test

T = []
X = []
Y = []
Z = []

# drag force
def projectile_force (bod, q, u, time, step):
  f = ( -bod.mass * u [3],
        -bod.mass * u [4],
	-bod.mass * u [5],
	 0, 0, 0, 0, 0, 0)
  return f

# plotting callback
def callback_function (bod, solfec):
  T.append (solfec.time)
  X.append (bod.conf [9] / 0.0254)
  Y.append (bod.conf [10] / 0.0254)
  Z.append (bod.conf [11] / 0.0254)
  return 1

# main module
#import rpdb2; rpdb2.start_embedded_debugger('a')

c = (0, 0, 0)
J = (1, 0, 0, 0, 1, 0, 0, 0, 1)
v0 = (2.54, 0, 12.7)
w0 = (0, 0, 0)
gravity = 9.81456;
stop = 1.976
step = stop / 1024.0

solfec = SOLFEC ('DYNAMIC', step, 'out/tests/projectile')
solfec.verbose = 'OFF'
material = BULK_MATERIAL (solfec)
GRAVITY (solfec, (0, 0, -1), gravity)

sph = SPHERE (c, 1, 0, 0)
bod = BODY (solfec, 'RIGID', sph, material);
BODY_CHARS (bod, 0.45359237, 1, c, J);
FORCE (bod, 'SPATIAL', (0,0,0), (0,0,0), projectile_force, bod);
INITIAL_VELOCITY (bod, v0, w0)

if not VIEWER():
  if solfec.mode == 'READ': print '\nPrevious test results exist. Please "make del" and rerun tests'
  else:
    CALLBACK (solfec, step, (bod, solfec), callback_function)
    RUN (solfec, PENALTY_SOLVER(), stop)
    value = bod.conf [9] / 0.0254 # conver to inches
    exact = 86.138
    error = abs (value - exact) / exact
    if (error < 0.005): print 'PASSED'
    else:
      print 'FAILED'
      print '(', 'Computed x distance was %.3f' % value, 'while the reference value is %.3f' % exact, ')'

    try:
      import matplotlib.pyplot as plt
      plt.clf ()
      plt.plot (T, X, label='x')
      plt.plot (T, Y, label='y')
      plt.plot (T, Z, label='z')
      plt.axis (xmin = 0, xmax = stop, ymin = 0, ymax = 180)
      plt.xlabel ('Time [s]')
      plt.ylabel ('Position [in]')
      plt.legend (loc = 'upper right')
      plt.savefig ('out/tests/projectile/projectile.eps')
    except ImportError:
      pass # no reaction

else: RUN (solfec, PENALTY_SOLVER(), stop)
