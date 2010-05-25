# bar drop example

def bar_drop (step, stop, scheme):

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

  bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)
  bod.scheme = scheme
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

stop = 1
th1 = bar_drop (1E-3, stop, 'DEF_IMP')
th2 = bar_drop (1E-4, stop, 'DEF_IMP')
th3 = bar_drop (1E-5, stop, 'DEF_EXP')

if not VIEWER() and isinstance (th1, tuple):
  try:
    import matplotlib.pyplot as plt
    plt.plot (th1 [0], th1 [1], label='KIN-IMP(1E-3)')
    plt.plot (th2 [0], th2 [1], label='KIN-IMP(1E-4)')
    plt.plot (th3 [0], th3 [1], label='KIN-EXP(1E-5)')
    plt.axis (xmin = 0, xmax = stop)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Energy [J]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/bar-drop/bar-drop-ene.eps')
    plt.clf ()
    plt.plot (th1 [0], th1 [4], label='VZ-IMP(1E-3)')
    plt.plot (th2 [0], th2 [4], label='VZ-IMP(1E-4)')
    plt.plot (th3 [0], th3 [4], label='VZ-EXP(1E-5)')
    plt.ylabel ('Velocity [m/s]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/bar-drop/bar-drop-vz.eps')
    plt.clf ()
    plt.plot (th1 [0], th1 [5], label='DZ-IMP(1E-3)')
    plt.plot (th2 [0], th2 [5], label='DZ-IMP(1E-4)')
    plt.plot (th3 [0], th3 [5], label='DZ-EXP(1E-5)')
    plt.ylabel ('Displacement [m]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/bar-drop/bar-drop-dz.eps')
 
  except ImportError:
    pass # no reaction
