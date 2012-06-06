# two-bar impact example from:
# F. Cirak and M. West, Decomposition Contact Response (DCR) for Explicit Finite Element Dynamics,
# International Journal for Numerical Methods in Engineering, 64(8):1078-1110, 2005. 

def bar_drop (step, stop, formul, scheme, damp, tet):

  a = 1.0
  b = 1.0
  c = 10.0
  V0 = 0.1

  nodes = [-a, -b, 0,
	    a, -b, 0,
	    a,  b, 0,
	   -a,  b, 0,
	   -a, -b, c,
	    a, -b, c,
	    a,  b, c,
	   -a,  b, c]

 
  if tet: end = '_tet'
  else: end = '_hex'
  sol = SOLFEC ('DYNAMIC', step, 'out/twobarimp/' + scheme + '-' + str (step) + end)

  SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = 0.0)

  bulk = BULK_MATERIAL (sol,
			model = 'KIRCHHOFF',
			young = 1.0,
			poisson = 0.0,
			density = 1.0)

  msh = HEX (nodes, 1, 1, 100, 0, [0, 0, 0, 0, 0, 0])
  if tet: msh = TETRAHEDRALIZE (msh, 'out/twobarimp/t1')
  b1 = BODY (sol, 'FINITE_ELEMENT', msh, bulk, form = formul)
  b1.scheme = scheme
  b1.damping = damp
  INITIAL_VELOCITY (b1, (0, 0, -V0), (0, 0, 0))

  msh = HEX (nodes, 1, 1, 100, 0, [0, 0, 0, 0, 0, 0])
  TRANSLATE (msh, (0, 0, -10))
  if tet: msh = TETRAHEDRALIZE (msh, 'out/twobarimp/t2')
  b2 = BODY (sol, 'FINITE_ELEMENT', msh, bulk, form = formul)
  b2.scheme = scheme
  b2.damping = damp
  INITIAL_VELOCITY (b2, (0, 0, V0), (0, 0, 0))

  gs = GAUSS_SEIDEL_SOLVER (1E-10, 1000)
  OUTPUT (sol, step)
  RUN (sol, gs, stop)

  if not VIEWER() and sol.mode == 'READ':
    return HISTORY (sol, [(sol, 'KINETIC'), (sol, 'INTERNAL'), (sol, 'EXTERNAL'), (b1, (0, 0, 0), 'VZ'), (b1, (0, 0, 0), 'DZ'), (b2, (0, 0, 0), 'VZ'), (b2, (0, 0, 0), 'DZ')], 0, stop)
  else: return [sol, gs]

stop = 50
damp = 0.0
formul = 'BC'
th1 = bar_drop (0.05, stop, formul, 'DEF_LIM', damp, 0)
th2 = bar_drop (0.05, stop, formul, 'DEF_LIM', damp, 1)

if not VIEWER() and isinstance (th1, tuple):
  try:
    import matplotlib.pyplot as plt
    plt.title ('Two-bar impact energy')
    plt.plot (th1 [0], th1 [1], label='KIN-hex')
    plt.plot (th1 [0], th1 [2], label='INT-hex')
    plt.plot (th2 [0], th2 [1], label='KIN-tet')
    plt.plot (th2 [0], th2 [2], label='INT-tet')
    plt.axis (xmin = 0, xmax = stop)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Energy [J]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/twobarimp/twobarimp_ene.eps')
    plt.clf ()
    plt.title ('Two-bar impact velocity')
    plt.plot (th1 [0], th1 [4], label='bar1-hex')
    plt.plot (th1 [0], th1 [6], label='bar2-hex')
    plt.plot (th2 [0], th2 [4], label='bar1-tet')
    plt.plot (th2 [0], th2 [6], label='bar2-tet')
    plt.ylabel ('Velocity [m/s]')
    plt.legend(loc = 'upper left')
    plt.savefig ('out/twobarimp/twobarimp_vz.eps')
    plt.clf ()
    plt.title ('Two-bar impact displacement')
    plt.plot (th1 [0], th1 [5], label='bar1-hex')
    plt.plot (th1 [0], th1 [7], label='bar2-hex')
    plt.plot (th2 [0], th2 [5], label='bar1-tet')
    plt.plot (th2 [0], th2 [7], label='bar1-tet')
    plt.ylabel ('Displacement [m]')
    plt.legend(loc = 'upper left')
    plt.savefig ('out/twobarimp/twobarimp_dz.eps')
 
  except ImportError:
    pass # no reaction
