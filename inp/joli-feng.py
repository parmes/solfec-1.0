# example from:
# Joli, P. and Feng, Z.-Q.,
# Uzawa and Newton algorithms to solve frictional contact problems within the bi-potential framework,
# International Journal for Numerical Methods in Engineering (2008), pp. 317-330
#
# It seem like it's going to be hard to get this to work, without prescribing the motion by reducing
# the size of the configuration space (handling all constraints for this tiny, stiff cube seems too
# hard for the solvers)

c = 1E-3
d = 0.2 * c
step = 1E-3
stop = 10 * step

#sv = NEWTON_SOLVER (1E-10, 100)
sv = GAUSS_SEIDEL_SOLVER (1E-5, 1000, 1E-5)

GEOMETRIC_EPSILON (1E-8)

solfec = SOLFEC ('QUASI_STATIC', step, 'out/joli-feng')

table = HEX ([-d, -d,  0,
	       c+d, -d,  0,
	       c+d,  c+d,  0,
	      -d,  c+d,  0,
              -d, -d, -d,
	       c+d, -d, -d,
	       c+d,  c+d, -d,
	      -d,  c+d, -d], 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 210E9, poisson = 0.3, density = 10E3)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3)

BODY (solfec, 'OBSTACLE', table, bulkmat)

box = HEX ([0, 0, 0,
	    c, 0, 0,
	    c, c, 0,
	    0, c, 0,
            0, 0, c,
	    c, 0, c,
	    c, c, c,
	    0, c, c], 2, 2, 2, 2, [2, 2, 2, 2, 2, 2])

point = []
for i in range (18, 27): point.append (box.node (i))

q = []
for i in range (0, 9): q.append (box.node (i))

bod = BODY (solfec, 'FINITE_ELEMENT', box, bulkmat, form = 'BC')
for i in range (0, 9): SET_VELOCITY (bod, point [i], (0, 0, -1), 1E-3)
for i in range (0, 9): FIX_DIRECTION (bod, point [i], (1, 0, 0))
for i in range (0, 9): FIX_DIRECTION (bod, point [i], (0, 1, 0))

RUN (solfec, sv, stop)
