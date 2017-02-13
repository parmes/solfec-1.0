# drag force
def projectile_force (bod, q, u, time, step):
  f = ( -bod.mass * u [3],
        -bod.mass * u [4],
	-bod.mass * u [5],
	 0, 0, 0, 0, 0, 0)
  return f

T = []
X = []
Y = []
Z = []

# plotting callback
def callback_function (bod, solfec):
  T.append (solfec.time)
  X.append (bod.conf [9] / 0.0254)
  Y.append (bod.conf [10] / 0.0254)
  Z.append (bod.conf [11] / 0.0254)
  return 1

c = (0, 0, 0)
J = (1, 0, 0, 0, 1, 0, 0, 0, 1)
v0 = (2.54, 0, 12.7)
w0 = (0, 0, 0)
gravity = 9.81456;
stop = 1.976
step = stop / 1024.0

solfec = SOLFEC ('DYNAMIC', step, 'out/projectile')
material = BULK_MATERIAL (solfec)
GRAVITY (solfec, (0, 0, -gravity))

sph = SPHERE (c, 1, 0, 0)
bod = BODY (solfec, 'RIGID', sph, material);
BODY_CHARS (bod, 0.45359237, 1, c, J);
FORCE (bod, 'SPATIAL', (0,0,0), (0,0,0), projectile_force, bod);
INITIAL_VELOCITY (bod, v0, w0)

CALLBACK (solfec, step, (bod, solfec), callback_function)
RUN (solfec, PENALTY_SOLVER(), stop)
value = bod.conf [9] / 0.0254 # conver to inches
exact = 86.138
error = abs (value - exact) / exact
print 'Computed x distance: %.3f;' % value, 'Reference distance: %.3f;' % exact, 'Relative error: %.3f' % error

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
  plt.savefig ('validation/projectile/projectile.png')
except (ImportError, RuntimeError):
  import sys
  print "Unexpected error:", sys.exc_info()[1]
  print "Plotting has failed!"
  pass
