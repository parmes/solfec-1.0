# example from:
# Joli, P. and Feng, Z.-Q.,
# Uzawa and Newton algorithms to solve frictional contact problems within the bi-potential framework,
# International Journal for Numerical Methods in Engineering (2008), pp. 317-330

c = 1E-3
d = 0.2 * c

step = 1E-3

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
pp = box.node (4)

bod = BODY (solfec, 'FINITE_ELEMENT', box, bulkmat)
for i in range (0, 9): SET_VELOCITY (bod, point [i], (0, 0, -1), 1E-3)

sv = NEWTON_SOLVER (1E-10, 100, 'NONSMOOTH_HYBRID')
sv.linminiter = 20
sv.resdec = 0.001
#sv = GAUSS_SEIDEL_SOLVER (1E-5, 100, 1E-5)

RUN (solfec, sv, 100 * step)
