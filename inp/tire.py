# bent mesh example

step = 0.001
stop = 1.0

sol = SOLFEC ('QUASI_STATIC', step, 'out/tire')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E7,
		      poisson = 0.45,
		      density = 1E3)

SURFACE_MATERIAL (sol, 0, 1, model = 'SIGNORINI_COULOMB', label = 'A', friction = 1, cohesion = 1E6)
SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', label = 'B', friction = 0.3)

nodes = [0, 0, 0,
         1, 0, 0,
         1, 1, 0,
         0, 1, 0,
         0, 0, 1,
         1, 0, 1,
         1, 1, 1,
         0, 1, 1]

msh = PIPE ((0, 0, 0), (0, 0, -1.88495559215), 0.3, 0.1, 10, 10, 1, 1, [1, 1, 2, 3])
BEND (msh, (0, 1, 0), (1, 0, 0), 180)
bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)
PRESSURE (bod, 2, -1E3)

shp = HEX (nodes, 1, 1, 1, 1, [0, 0, 0, 0, 0, 0])
SCALE (shp, (2, 4, 1))
TRANSLATE (shp, (-1, -1, 0))
bod = BODY (sol, 'RIGID', shp, bulk)
ts0 = TIME_SERIES ([0, 1, 0.5, 1, 0.5+step, 0, 1.0, 0])
ts1 = TIME_SERIES ([0, 0, 0.5, 0, 0.5+step, 1, 1.0, 1])
SET_VELOCITY (bod, (-1, -1, 1), (0, 0, -1), ts0)
SET_VELOCITY (bod, (1, -1, 1), (0, 0, -1), ts0)
SET_VELOCITY (bod, (-1, 3, 1), (0, 0, -1), ts0)
SET_VELOCITY (bod, (1, 3, 1), (0, 0, -1), ts0)
SET_VELOCITY (bod, bod.center, (0, 1, 0), ts1)

shp = HEX (nodes, 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])
SCALE (shp, (2, 4, 0.1))
TRANSLATE (shp, (-1, -1, -1.5))
bod = BODY (sol, 'OBSTACLE', shp, bulk)

#sv = GAUSS_SEIDEL_SOLVER (1E-6, 1000)
sv = NEWTON_SOLVER (locdyn = 'OFF')
OUTPUT (sol, step)
RUN (sol, sv, stop)

if not VIEWER() and sol.mode == 'READ':

  timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL', 'PARBAL']
  dur = DURATION (sol)
  th = HISTORY (sol, timers, dur[0], dur[1])
  total = 0.0

  for i in range (0, len(timers)):
    sum = 0.0
    for tt in th [i+1]: sum += tt
    print timers [i], 'TIME:', sum
    if i < 6: total += sum

  print 'TOTAL TIME:', total

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
    plt.savefig ('out/tire/ene.eps')
  except ImportError:
    pass # no reaction
