from math import sin, cos
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

def double_pendulum_run (step, stop, gravity):
  solfec = SOLFEC ('DYNAMIC', step, 'out/double-pendulum')
  bulkmat = BULK_MATERIAL (solfec)
  SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', restitution = 0.1, friction = 0.0)
  GRAVITY (solfec, (0, 0, -gravity))
  double_pendulum_create (bulkmat, solfec)
  gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'EXIT')
  RUN (solfec, gs, stop)
  return (gs, solfec) # return solver and solfec to avoid their deallocation

gravity = 9.81
step = 0.001
stop = 2.5

data = double_pendulum_run (step, stop, gravity)

if not VIEWER():
  if data[1].mode == 'WRITE':
    solfec = double_pendulum_run (step, stop, gravity) [1] # reopen in read mode
    if solfec.mode == 'WRITE': print 'ERROR --> delete results (e.g. make del) and rerun'
  else: solfec = data[1]

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

  try:
    import matplotlib.pyplot as plt
    plt.clf ()
    plt.plot (th [0], th [6], label='$x_1(t)$')
    plt.plot (th [0], th [7], label='$x_2(t)$')
    plt.axis (xmin = 0, xmax = stop, ymin = 0, ymax = 1.5)
    plt.yticks([0.0, 0.5, 1.0, 1.5])
    plt.xlabel ('Time [s]')
    plt.ylabel ('Position [m]')
    plt.legend(loc = 'upper right')
    plt.title('Solfec')
    plt.savefig ('validation/double-pendulum/double-solfec-x.png')
    plt.clf ()
    plt.plot (th [0], tot)
    plt.axis (xmin = 0, xmax = stop, ymin = 0, ymax = 12)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Energy [J]')
    plt.title('Solfec')
    plt.savefig ('validation/double-pendulum/double-solfec-energy.png')
  except (ImportError, RuntimeError):
    import sys
    print "Unexpected error:", sys.exc_info()[1]
    print "Plotting has failed!"
    pass
