# bent mesh example

step = 0.002
stop = 0.2

sol = SOLFEC ('DYNAMIC', step, 'out/tire')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E7,
		      poisson = 0.45,
		      density = 1E3)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

nodes = [0, 0, 0,
         1, 0, 0,
         1, 1, 0,
         0, 1, 0,
         0, 0, 1,
         1, 0, 1,
         1, 1, 1,
         0, 1, 1]

msh = PIPE ((0, 0, 0), (0, 1.5, 0), 0.3, 0.1, 10, 20, 2, 1, [1, 2, 3, 4])
BEND (msh, (0, 0, 1), (1, 0, 0), 180)
bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)
bod.nodecontact = 'ON'
PRESSURE (bod, 3, -10000)
PARTITION (bod, NCPU (sol)) #FIXME: in parallel this high presure separates meshes!

shp = HEX (nodes, 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])
SCALE (shp, (2, 4, 0.1))
TRANSLATE (shp, (-1, -1, -0.5))
bod = BODY (sol, 'OBSTACLE', shp, bulk)

gs = GAUSS_SEIDEL_SOLVER (1E-6 1000)
GRAVITY (sol, (0, 0, -10))
OUTPUT (sol, step)
RUN (sol, gs, stop)

if not VIEWER() and sol.mode == 'READ':

  timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL', 'PARBAL',
            'GSINIT', 'GSRUN', 'GSCOM', 'GSMCOM']
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
    plt.savefig ('out/bend/bendene.eps')
  except ImportError:
    pass # no reaction
