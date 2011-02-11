# finite elements test

step = 1E-3

solfec = SOLFEC ('DYNAMIC', step, 'out/domino')
solver = NEWTON_SOLVER ()
bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E3, poisson = 0.3, density = 10)

tet = MESH ([0, 0, 0,
             1, 0, 0,
	     0, 1, 0,
	     0, 0, 1],
	     [4, 0, 1, 2, 3, 1], 1)

pyr = MESH ([0, 0, 0,
             1, 0, 0,
	     1, 1, 0,
	     0, 1, 0,
             0.5, 0.5, 1],
	     [5, 0, 1, 2, 3, 4, 2], 2)

TRANSLATE (pyr, (2, 0, 0))

wed = MESH ([0, 0, 0,
             1, 0, 0,
	     0, 1, 0,
             0, 0, 1,
             1, 0, 1,
	     0, 1, 1],
	     [6, 0, 1, 2, 3, 4, 5, 3], 3)

TRANSLATE (wed, (0, 2, 0))

hex = MESH ([0, 0, 0,
             1, 0, 0,
	     1, 1, 0,
	     0, 1, 0,
             0, 0, 1,
             1, 0, 1,
	     1, 1, 1.1,
	     0, 1, 1.1],
	     [8, 0, 1, 2, 3, 4, 5, 6, 7, 4], 4)

#TRANSLATE (hex, (2, 2, 0))

#BODY (solfec, 'FINITE_ELEMENT', tet, bulkmat)
#BODY (solfec, 'FINITE_ELEMENT', pyr, bulkmat)
#BODY (solfec, 'FINITE_ELEMENT', wed, bulkmat)


shp = HULL ([0, 0, 0,
             1, 0, 0,
	     0, 1, 0,
	     0, 0, 1], 0, 0)

BODY (solfec, 'FINITE_ELEMENT', shp, bulkmat, mesh = hex)

table = HULL ([0, 0, 0,
               0, 1, 0,
	       1, 1, 0,
	       1, 0, 0,
               0, 0, -0.1,
               0, 1, -0.1,
	       1, 1, -0.1,
	       1, 0, -0.1], 1, 1)

SCALE (table, (5, 5, 1))
TRANSLATE (table, (-1, -1, 0))

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5)

BODY (solfec, 'OBSTACLE', table, bulkmat)

GRAVITY (solfec, (0, 0, -10))

RUN (solfec, solver, 1.0)
