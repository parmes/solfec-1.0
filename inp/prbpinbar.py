# pinned bar example

def pinned_bar_create (material, solfec, kind, scheme):

  nodes = [-0.05, -0.05, 0.0,
            0.05, -0.05, 0.0,
            0.05,  0.05, 0.0,
	   -0.05,  0.05, 0.0,
	   -0.05, -0.05, 1.0,
	    0.05, -0.05, 1.0,
	    0.05,  0.05, 1.0,
	   -0.05,  0.05, 1.0]

  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, 0]

  point = (0, 0, 0.75)

  vector = (0, 1, 0)

  fix1 = (0, -0.05, 0.75)

  fix2 = (0, 0.05, 0.75)

  msh = MESH (nodes, elements, 0)

  ROTATE (msh, point, vector, -30)

  bod = BODY (solfec, kind, msh, material)
  bod.scheme = scheme

  FIX_POINT (solfec, bod, fix1)
  FIX_DIRECTION (solfec, bod, fix2, (1, 0, 0))
  FIX_DIRECTION (solfec, bod, fix2, (0, 0, 1))

  msh = MESH (nodes, elements, 0)

  vector = (-0.15, 0, -0.2)

  TRANSLATE (msh, vector)

  BODY (solfec, 'OBSTACLE', msh, material)

  return bod

### solfec context ###

def create_simulation (output, step, stop, solver, kind, scheme):

  solfec = SOLFEC ('DYNAMIC', step, output)
  surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = 1.0)
  bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)
  GRAVITY (solfec, (0, 0, -1), 9.8)
  bod = pinned_bar_create (bulkmat, solfec, kind, scheme)
  RUN (solfec, solver, stop)

  return (solfec, bod)

### main module ###
#import rpdb2; rpdb2.start_embedded_debugger('a')

step = 0.001
stop = 1

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'CONTINUE', diagsolver = 'PROJECTED_GRADIENT')

s1 = create_simulation ('out/prbpinbar/rig', step, stop, gs, 'RIGID', 'RIG_NEG')
s2 = create_simulation ('out/prbpinbar/imp', step, stop, gs, 'PSEUDO_RIGID', 'DEF_IMP')
s3 = create_simulation ('out/prbpinbar/exp', step, stop, gs, 'PSEUDO_RIGID', 'DEF_EXP')

if not VIEWER() and s1[0].mode == 'READ':
  try:
    import matplotlib.pyplot as plt
    th1 = HISTORY (s1[0], [(s1[1], 'KINETIC'), (s1[1], 'INTERNAL'), (s1[1], 'EXTERNAL'), (s1[1], 'CONTACT'), (s1[1], 'FRICTION')], 0, stop)
    th2 = HISTORY (s2[0], [(s2[1], 'KINETIC'), (s2[1], 'INTERNAL'), (s2[1], 'EXTERNAL'), (s2[1], 'CONTACT'), (s2[1], 'FRICTION')], 0, stop)
    th3 = HISTORY (s3[0], [(s3[1], 'KINETIC'), (s3[1], 'INTERNAL'), (s3[1], 'EXTERNAL'), (s3[1], 'CONTACT'), (s3[1], 'FRICTION')], 0, stop)
    plt.plot (th3 [0], th3 [1], label='PRB-EXP')
    plt.plot (th2 [0], th2 [1], label='PRB-IMP')
    plt.plot (th1 [0], th1 [1], label='RIG')
    plt.axis (xmin = 0, xmax = stop)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Kinetic energy [J]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/prbpinbar/prbpinbar.eps')
  except ImportError:
    pass # no reaction
