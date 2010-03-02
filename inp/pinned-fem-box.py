# pinned fem box example 

a = 1.0
b = 1.0
c = 2.0
step = 0.01 # let the critical step rule
stop = 10

nodes = [-a, -b, 0,
          a, -b, 0,
          a,  b, 0,
         -a,  b, 0,
         -a, -b, c,
          a, -b, c,
          a,  b, c,
         -a,  b, c]

msh = HEX (nodes, 2, 2, 2, 0, [0, 1, 2, 3, 4, 5])

sol = SOLFEC ('DYNAMIC', step, 'out/pinned-fem-box')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 15E9,
		      poisson = 0.25,
		      density = 1.8E3)

bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)
bod.scheme = 'DEF_LIM'
FIX_POINT (bod, (-a, -b, 0))
FIX_POINT (bod, (-a, b, 0))

gs = GAUSS_SEIDEL_SOLVER (1E-5, 1000)
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
    plt.savefig ('out/pinned-fem-box/pinned-fem-box.eps')
  except ImportError:
    pass # no reaction
