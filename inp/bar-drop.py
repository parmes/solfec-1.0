# bar drop example

a = 0.1
b = 0.1
c = 1.0
step = 0.001
stop = 5
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

sol = SOLFEC ('DYNAMIC', step, 'out/bar-drop')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 15E9,
		      poisson = 0.25,
		      density = 1.8E3)

bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)
bod.scheme = 'DEF_LIM'
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
  try:
    import matplotlib.pyplot as plt
    th = HISTORY (sol, [(sol, 'KINETIC'), (sol, 'INTERNAL'), (sol, 'EXTERNAL')], 0, stop)
    plt.plot (th [0], th [1], label='KIN')
    plt.plot (th [0], th [2], label='INT')
    plt.plot (th [0], th [3], label='EXT')
    tot = []
    for i in range(0, len (th[0])): tot.append (th[1][i] + th[2][i] - th[3][i])
    plt.plot (th [0], tot, label='TOT')
    plt.axis (xmin = 0, xmax = stop)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Energy [J]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/bar-drop/bar-drop.eps')
  except ImportError:
    pass # no reaction
