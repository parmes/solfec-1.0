# bar drop example

def bar_drop (step, stop, formul, scheme, damp):

  a = 0.1
  b = 0.1
  c = 1.0
  V0 = 1.0

  nodes = [-a, -b, 0,
	    a, -b, 0,
	    a,  b, 0,
	   -a,  b, 0,
	   -a, -b, c,
	    a, -b, c,
	    a,  b, c,
	   -a,  b, c]

  msh = HEX (nodes, 1, 1, 5, 0, [0, 0, 0, 0, 0, 0])

  sol = SOLFEC ('DYNAMIC', step, 'out/bar-drop/' + scheme + '-' + str (step))

  bulk = BULK_MATERIAL (sol,
			model = 'KIRCHHOFF',
			young = 15E9,
			poisson = 0.25,
			density = 1.8E3)

  bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk, form = formul)
  bod.scheme = scheme
  bod.damping = damp
  INITIAL_VELOCITY (bod, (0, 0, -V0), (0, 0, 0))

  a = 2 * a
  b = 2 * b
  c = 0.1
  nodes = [-a, -b, -c,
	    a, -b, -c,
	    a,  b, -c,
	   -a,  b, -c,
	   -a, -b, 0,
	    a, -b, 0,
	    a,  b, 0,
	   -a,  b, 0]

  msh = HEX (nodes, 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])
  BODY (sol, 'OBSTACLE', msh, bulk)

  SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = 0.0)

  gs = GAUSS_SEIDEL_SOLVER (1E-10, 1000, diagsolver = 'PROJECTED_GRADIENT')
  GRAVITY (sol, (0, 0, -10))
  OUTPUT (sol, step)
  RUN (sol, gs, stop)

  if not VIEWER() and sol.mode == 'READ':
    return HISTORY (sol, [(sol, 'KINETIC'), (sol, 'INTERNAL'), (sol, 'EXTERNAL'), (bod, (0, 0, 0), 'VZ'), (bod, (0, 0, 0), 'DZ')], 0, stop)
  else: return [sol, gs]

step = 1E-4
stop = 1
damp = 2E-5
formul = 'BC'
th1 = bar_drop (step, stop, formul, 'DEF_LIM', damp)
th2 = bar_drop (step, stop, formul, 'DEF_LIM2', damp)
th3 = bar_drop (step, stop, formul, 'DEF_IMP', damp)

if not VIEWER() and isinstance (th1, tuple):
  try:
    import matplotlib.pyplot as plt
    plt.title ('ENE')
    plt.plot (th1 [0], th1 [1], label='LIM')
    plt.plot (th2 [0], th2 [1], label='LIM2')
    plt.plot (th3 [0], th3 [1], label='IMP')
    plt.axis (xmin = 0, xmax = stop)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Energy [J]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/bar-drop/bar-drop-ene.eps')
    plt.clf ()
    plt.title ('VZ')
    plt.plot (th1 [0], th1 [4], label='LIM')
    plt.plot (th2 [0], th2 [4], label='LIM2')
    plt.plot (th3 [0], th3 [4], label='IMP')
    plt.ylabel ('Velocity [m/s]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/bar-drop/bar-drop-vz.eps')
    plt.clf ()
    plt.title ('DZ')
    plt.plot (th1 [0], th1 [5], label='LIM')
    plt.plot (th2 [0], th2 [5], label='LIM2')
    plt.plot (th3 [0], th3 [5], label='IMP')
    plt.ylabel ('Displacement [m]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/bar-drop/bar-drop-dz.eps')
 
  except ImportError:
    pass # no reaction
