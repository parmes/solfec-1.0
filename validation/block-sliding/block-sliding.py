from math import cos

# sliding force
def block_force (bod, q, u, time, step):
  f = (8 * cos (time), 0, 0, 0, 0, 0, 0, 0, 0)
  return f

# analysis
def block_sliding_run (step, stop):
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

  solfec = SOLFEC ('DYNAMIC', step, 'out/block-sliding')
  material = BULK_MATERIAL (solfec, density = 111.11111111111111111) # so that 0.3*0.3*0.1*111.(1) = 1.0
  SURFACE_MATERIAL (solfec, friction = 0.8)
  gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000)
  GRAVITY (solfec, (0, 0, -9.81))

  msh = HEX (base, 1, 1, 1, 0, surfaces)
  TRANSLATE (msh, vector)
  BODY (solfec, 'OBSTACLE', msh, material)
  msh = HEX (cube, 2, 2, 1, 0, surfaces)
  TRANSLATE (msh, vector)
  bod = BODY (solfec, 'RIGID', msh, material)
  FORCE (bod, 'SPATIAL', (0,0,0), (0,0,0), block_force, bod);
  RUN (solfec, gs, stop)
  return (bod, solfec)

# main module
step = 0.001
stop = 10.0
(bod, solfec) = block_sliding_run (step, stop)

if not VIEWER():
  if solfec.mode == 'WRITE':
    (bod, solfec) = block_sliding_run (step, stop) # reopen in read mode
    if solfec.mode == 'WRITE': print 'ERROR --> delete results (e.g. make del) and rerun'
  th = HISTORY (solfec, [(bod, bod.center, 'CX'), (bod, bod.center, 'VX')], 0, stop)
  try:
    import matplotlib.pyplot as plt
    plt.clf ()
    plt.plot (th[0], th[1], label='x')
    plt.axis (xmin = 0, xmax = stop, ymin = 2.99, ymax = 3.006)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Position [m]')
    plt.legend (loc = 'upper right')
    plt.title ('Solfec')
    plt.savefig ('validation/block-sliding/block-sliding-x.png')
    plt.clf ()
    plt.plot (th[0], th[2], label='$v_x$')
    plt.axis (xmin = 0, xmax = stop, ymin = -0.04, ymax = 0.04)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Velocity [m/s]')
    plt.legend (loc = 'upper right')
    plt.title ('Solfec')
    plt.savefig ('validation/block-sliding/block-sliding-v.png')
  except (ImportError, RuntimeError):
    import sys
    print "Unexpected error:", sys.exc_info()[1]
    print "Plotting has failed!"
    pass
