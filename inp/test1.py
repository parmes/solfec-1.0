# test of Solfec functionality with geometry defined by convex polyhedra
import math
from sys import getrefcount


a = 1.0
b = 1.0
c = 1.0
step = 1E-3
time = 10 * step

shape = CONVEX (
        [0, 0, 0,
	 a, 0, 0,
	 a, b, 0,
	 0, b, 0,
	 0, 0, c,
	 a, 0, c,
	 a, b, c,
	 0, b, c],
	[4, 0, 3, 2, 1, 2,
	 4, 1, 2, 6, 5, 2,
	 4, 2, 3, 7, 6, 2,
	 4, 3, 0, 4, 7, 2,
	 4, 0, 1, 5, 4, 2,
	 4, 4, 5, 6, 7, 2], 2)

#shape = HULL (
#       [0, 0, 0,
#	 a, 0, 0,
#	 a, b, 0,
#	 0, b, 0,
#	 0, 0, c,
#	 a, 0, c,
#	 a, b, c,
#	 0, b, c], 1, 1)

sol = SOLFEC ('DYNAMIC', step, 'out/test1a')

sol2 = SOLFEC ('DYNAMIC', step, 'out/test1b')

sur = SURFACE_MATERIAL (sol,
                        model = 'SIGNORINI_COULOMB',
                        friction = 0.5,
			spring = 1E3,
			dashpot = 1E2)

sur2 = SURFACE_MATERIAL (sol2,
                        model = 'SIGNORINI_COULOMB',
                        friction = 0.3)

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E4,
		      poisson = 0.25,
		      density = 1E2)

bulk2 = BULK_MATERIAL (sol2,
                      model = 'KIRCHHOFF',
		      young = 5E6,
		      poisson = 0.25,
		      density = 1E3)


point = (a/2, b/2, c/2)
normal = (1, 0, 0)

copy = TRANSLATE (COPY (shape), (0, 0, 0))
(one, two) = SPLIT (copy, point, normal)
(b1, b2) = SPLIT (COPY (one), point, (0, 1, 0))
(b3, b4) = SPLIT (COPY (two), point, (0, 1, 0))

copy = TRANSLATE ([b1, b2, b3, b4], (0, 0, 1.5))
BODY (sol, 'OBSTACLE', shape, bulk)
BODY (sol, 'PSEUDO_RIGID', copy, bulk)

#BODY (sol, 'RIGID', [one, two], bulk)

#BODY (sol, 'RIGID', shape, bulk)

BODY (sol2, 'RIGID', one, bulk2)
BODY (sol2, 'RIGID', two, bulk2)

def gscallback (gs):
  print gs.error
  return 0

gs = GAUSS_SEIDEL_SOLVER (1E-3, 10000, failure = 'CALLBACK', callback=gscallback)
gs.history = 'ON'

ex = EXPLICIT_SOLVER ()

GRAVITY (sol, (0, 0, -1), 10)

OUTPUT (sol, 10 * step)

def runcallback (gs):
  if gs.error != 'OK': print gs.error
  else: print 'gs iters = ', gs.iters
  return 1

print '@@@ gs refcnt 1 =', getrefcount (gs)

CALLBACK (sol, step, (gs), runcallback)

print '@@@ gs refcnt 2 =', getrefcount (gs)

RUN (sol, gs, time)
