# ellipsoids in a box

step = 1E-3
stop = 10

solfec = SOLFEC ('DYNAMIC', step, 'out/ellip')

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E6, poisson = 0.3, density = 1E3)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0, restitution = 1, spring = 1E6, dashpot = 0)

side = HULL ([0, 0, 0,
              0, 1, 0,
	      1, 1, 0,
	      1, 0, 0,
              0, 0, -0.1,
              0, 1, -0.1,
	      1, 1, -0.1,
	      1, 0, -0.1], 1, 1)
SCALE (side, (0.1, 0.1, 0.05))
BODY (solfec, 'OBSTACLE', COPY (side), bulkmat)
BODY (solfec, 'OBSTACLE', TRANSLATE (COPY (side), (0, 0, 0.1)), bulkmat)
BODY (solfec, 'OBSTACLE', ROTATE (COPY (side), (0, 0, 0), (1, 0, 0), 90), bulkmat)
BODY (solfec, 'OBSTACLE', TRANSLATE (ROTATE (COPY (side), (0, 0, 0), (1, 0, 0), 90), (0, 0.095, 0)), bulkmat)
BODY (solfec, 'OBSTACLE', ROTATE (COPY (side), (0, 0, 0), (0, 1, 0), -90), bulkmat)
BODY (solfec, 'OBSTACLE', TRANSLATE (ROTATE (COPY (side), (0, 0, 0), (0, 1, 0), -90), (0.095, 0, 0)), bulkmat)

n = 1
for i in range (0, n):
  for j in range (0, n):
    for k in range (0, n):
      s = 0.095/float(n+1)
      x = s + i*s
      y = s + j*s
      z = s + k*s
      shp = ELLIP ((x, y, z), (0.2*s, 0.3*s, 0.4*s), 3, 3)
      #shp = SPHERE ((x, y, z), 0.2*s, 3, 3)
      ROTATE (shp, (x, y, z), (i+1, j+1, k+1), 7.0*(i+1)*(j+1)*(k+1))
      bod = BODY (solfec, 'RIGID', shp, bulkmat)
      len = (x*x+y*y+z*z)**0.5;      
      INITIAL_VELOCITY (bod, (0.01*x/len, -0.005*y/len, 0.07*z/len), (0, 0, 0))

#sv = PENALTY_SOLVER()
#sv = GAUSS_SEIDEL_SOLVER (1E-4, 1000)
sv = NEWTON_SOLVER ()

OUTPUT (solfec, 0.01)
RUN (solfec, sv, stop)

if not VIEWER() and solfec.mode == 'READ':
  try:
    import matplotlib.pyplot as plt
    th = HISTORY (solfec, [(solfec, 'KINETIC')], 0, stop)
    plt.plot (th [0], th [1], label='kinetic')
    plt.axis (xmin = 0, xmax = stop)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Kinetic energy [J]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/ellip/energy.eps')
  except ImportError:
    pass # no reaction
