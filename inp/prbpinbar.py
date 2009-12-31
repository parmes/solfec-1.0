# pinned bar example

def pinned_bar_create (material, solfec, scheme, base):

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

  point = TRANSLATE (point, base)

  vector = (0, 1, 0)

  fix1 = (0, -0.05, 0.75)

  fix1 = TRANSLATE (fix1, base)

  fix2 = (0, 0.05, 0.75)

  fix2 = TRANSLATE (fix2, base)

  msh = MESH (nodes, elements, 0)

  TRANSLATE (msh, base)

  ROTATE (msh, point, vector, -30)

  bod = BODY (solfec, 'PSEUDO_RIGID', msh, material)
  bod.scheme = scheme

  FIX_POINT (solfec, bod, fix1)
  FIX_DIRECTION (solfec, bod, fix2, (1, 0, 0))
  FIX_DIRECTION (solfec, bod, fix2, (0, 0, 1))

  msh = MESH (nodes, elements, 0)

  vector = (-0.15, 0, -0.2)

  vector = TRANSLATE (vector, base)

  TRANSLATE (msh, vector)

  BODY (solfec, 'OBSTACLE', msh, material)

  return bod

### solfec context ###

def create_simulation (output, step, stop, solver, scheme, base):

  solfec = SOLFEC ('DYNAMIC', step, output)
  surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = 1.0)
  bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)
  GRAVITY (solfec, (0, 0, -1), 9.8)
  bod = pinned_bar_create (bulkmat, solfec, scheme, base)
  RUN (solfec, solver, stop)

  return (solfec, bod)

### main module ###

step = 0.001

stop = 1.5

#import rpdb2; rpdb2.start_embedded_debugger('a')

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'CONTINUE', diagsolver = 'PROJECTED_GRADIENT')

s1 = create_simulation ('out/prbpinbar1', step, stop, gs, 'DEF_EXP', (0, 0, 0))

s2 = create_simulation ('out/prbpinbar2', step, stop, gs, 'DEF_IMP', (0, 1, 0))

if not VIEWER() and s1[0].mode == 'READ':
  import matplotlib.pyplot as plt
  th1 = HISTORY (s1[0], [(s1[1], 'KINETIC'), (s1[1], 'INTERNAL'), (s1[1], 'EXTERNAL'), (s1[1], 'CONTACT'), (s1[1], 'FRICTION')], 0.4, stop)
  th2 = HISTORY (s2[0], [(s2[1], 'KINETIC'), (s2[1], 'INTERNAL'), (s2[1], 'EXTERNAL'), (s2[1], 'CONTACT'), (s2[1], 'FRICTION')], 0.4, stop)
  plt.plot (th1 [0], th1 [1], label='kinE')
  plt.plot (th1 [0], th1 [2], label='intE')
  plt.plot (th2 [0], th2 [1], label='kinI')
  plt.plot (th2 [0], th2 [2], label='intI')
  plt.legend(loc=1)
  plt.show ()
