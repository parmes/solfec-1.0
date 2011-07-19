# strange silo
from random import randint
from random import random
from random import seed
from math import sin

# random convex particle
def particle (x, y, z, r, material, solfec):
  m = randint (8, 32)
  points = []
  r = 0.3 * r + 0.7 * r * random ()
  for n in range (m):
    points.append (x + r * (1.0 - random()) * 2.0)
    points.append (y + r * (1.0 - random()) * 2.0)
    points.append (z + r * (1.0 - random()) * 2.0)
  hull = HULL (points, 2, 2)
  bod = BODY (solfec, 'RIGID', hull, material)


# main module
step = 0.002
stop = 2

seed (1)

sol = SOLFEC ('DYNAMIC', step, 'out/silo')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E7,
		      poisson = 0.3,
		      density = 1E3)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

GRAVITY (sol, (0, 0, -10))

nodes = [0, 0, 0,
         1, 0, 0,
         1, 1, 0,
         0, 1, 0,
         0, 0, 1,
         1, 0, 1,
         1, 1, 1,
         0, 1, 1]

for i in range (0, 5):
  for j in range (0, 5):
    for k in range (0, 5):
	msh = HEX (nodes, 15, 10, 3, 0, [0, 0, 0, 0, 0, 0])
	SCALE (msh, (15, 10, 1))
	BEND (msh, (0, 0, -3), (-1, 0, 0), 270)
	BEND (msh, (5, 7, 0), (0, 0, 1), 90)
	TRANSLATE (msh, (i*13, j*12-2, k*9))
	bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)

for i in range (1, 12):
  for j in range (1, 11):
    for k in range (1, 12):
      particle (i*5, j*5, 35 + k*5, 2.5, bulk, sol)

shp = HEX (nodes, 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])
SCALE (shp, (40, 65, 1))
ROTATE (shp, (0, 0, 0), (0, 1, 0), 45)
TRANSLATE (shp, (-5, -5, -10))
shp2 = COPY (shp)
bod = BODY (sol, 'OBSTACLE', shp, bulk)

ROTATE (shp2, (23.9914, -5, -37.5772), (0, 0, 1), 180)
TRANSLATE (shp2, (15, 65, 0))
bod = BODY (sol, 'OBSTACLE', shp2, bulk)

shp3 = HULL ([67.2757, 60, -9.29289,
              38.9914, 60, -37.5772,
              23.9914, 60, -37.5772,
              -4.29289, 60, -9.29289,
              67.2757, 61, -9.29289,
              38.9914, 61, -37.5772,
              23.9914, 61, -37.5772,
              -4.29289, 61, -9.29289], 1, 1)
shp4 = TRANSLATE (COPY (shp3), (0, -65, 0))
bod = BODY (sol, 'OBSTACLE', shp3, bulk)
bod = BODY (sol, 'OBSTACLE', shp4, bulk)

shp5 = HULL ([67.2757, -5, -50,
            -4.29289, -5, -50,
            -4.29289, 60, -50,
             67.2757, 60, -50,
             67.2757, -5, -51,
            -4.29289, -5, -51,
            -4.29289, 60, -51,
             67.2757, 60, -51], 1, 1)
bod = BODY (sol, 'OBSTACLE', shp5, bulk)

#sv = GAUSS_SEIDEL_SOLVER (1E-4, 1000)
sv = NEWTON_SOLVER (delta = 1E-7)

OUTPUT (sol, 0.01)
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
    plt.savefig ('out/slo/energy.eps')
  except ImportError:
    pass # no reaction
