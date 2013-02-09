# bent mesh example
step = 0.002
stop = 1

sol = SOLFEC ('DYNAMIC', step, 'out/bend')
bulk = BULK_MATERIAL (sol, model = 'KIRCHHOFF', young = 1E7, poisson = 0.3, density = 1E3)
SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

nodes = [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1]
msh = HEX (nodes, 15, 10, 2, 0, [0, 0, 0, 0, 0, 0])
SCALE (msh, (15, 10, 1))
BEND (msh, (0, 0, -3), (-1, 0, 0), 270)
BEND (msh, (5, 7, 0), (0, 0, 1), 90)
bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)

shp = HEX (nodes, 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])
SCALE (shp, (15, 15, 1))
TRANSLATE (shp, (-1, -2, -10))
bod = BODY (sol, 'OBSTACLE', shp, bulk)

sv = NEWTON_SOLVER ()
GRAVITY (sol, (0, 0, -10))
OUTPUT (sol, step)
RUN (sol, sv, stop)

if not VIEWER() and sol.mode == 'READ':
  try:
    import matplotlib.pyplot as plt
    dur = DURATION (sol)
    th = HISTORY (sol, [(sol, 'KINETIC'), (sol, 'INTERNAL'), (sol, 'EXTERNAL'), (sol, 'CONTACT'), (sol, 'FRICTION')], dur [0], dur [1])
    plt.plot (th [0], th [1], label='KIN')
    plt.plot (th [0], th [2], label='INT')
    plt.plot (th [0], th [3], label='EXT')
    plt.plot (th [0], th [4], label='CON')
    plt.plot (th [0], th [5], label='FRI')
    plt.xlabel ('Time [s]')
    plt.ylabel ('Energy [J]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/bend/bendene.eps')
  except ImportError:
    pass # no reaction
