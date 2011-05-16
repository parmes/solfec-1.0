# bent mesh example

step = 0.002
stop = 1.0

sol = SOLFEC ('DYNAMIC', step, 'out/cohesion')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E7,
		      poisson = 0.45,
		      density = 1E3)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 1, cohesion = 5E5)

nodes = [0, 0, 0,
         1, 0, 0,
         1, 1, 0,
         0, 1, 0,
         0, 0, 1,
         1, 0, 1,
         1, 1, 1,
         0, 1, 1]

shp = HEX (nodes, 4, 8, 8, 1, [0, 1, 2, 3, 4, 5])
SCALE (shp, (2, 4, 8))
bod = BODY (sol, 'FINITE_ELEMENT', shp, bulk)
bod.nodecontact = 'ON'
bod.damping = 0.001
PRESSURE (bod, 4, -90000)

shp = HEX (nodes, 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])
SCALE (shp, (2.2, 4.2, 0.1))
TRANSLATE (shp, (-0.1, -0.1, -0.1))
bod = BODY (sol, 'OBSTACLE', shp, bulk)
bod.nodecontact = 'ON'

sv = NEWTON_SOLVER (locdyn = 'OFF')
GRAVITY (sol, (0, 0, -10))
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
    plt.savefig ('out/cohesion/ene.eps')
  except ImportError:
    pass # no reaction
