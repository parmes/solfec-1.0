c = 0.1
d = 0.2 * c
step = 1E-3
stop = 10 * step

sv = GAUSS_SEIDEL_SOLVER (1E-5, 1000, 1E-5)
#sv = NEWTON_SOLVER (maxiter = 200, theta = 0.15, locdyn = 'OFF')

solfec = SOLFEC ('QUASI_STATIC', step, 'out/joli-feng')

table = HEX ([-d, -d,  0,
	       c+d, -d,  0,
	       c+d,  c+d,  0,
	      -d,  c+d,  0,
              -d, -d, -d,
	       c+d, -d, -d,
	       c+d,  c+d, -d,
	      -d,  c+d, -d], 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E6, poisson = 0.3, density = 500)
surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3)

BODY (solfec, 'OBSTACLE', table, bulkmat)

box = HEX ([0, 0, 0,
	    c, 0, 0,
	    c, c, 0,
	    0, c, 0,
            0, 0, c,
	    c, 0, c,
	    c, c, c,
	    0, c, c], 8, 8, 8, 2, [2, 2, 2, 2, 2, 2])

bod = BODY (solfec, 'FINITE_ELEMENT', box, bulkmat)

GRAVITY (solfec, (0, 0, -1000))
RUN (solfec, sv, stop)
