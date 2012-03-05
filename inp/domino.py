# a domino example

step = 1E-3

solfec = SOLFEC ('DYNAMIC', step, 'out/domino')

GRAVITY (solfec, (0, 0, -9.81))

table = HULL ([0, 0, 0,
               0, 1, 0,
	       1, 1, 0,
	       1, 0, 0,
               0, 0, -0.1,
               0, 1, -0.1,
	       1, 1, -0.1,
	       1, 0, -0.1], 1, 1)

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5, restitution = 0.25)

BODY (solfec, 'OBSTACLE', table, bulkmat)

hex = HEX ([0, 0, 0,
	    1, 0, 0,
	    1, 1, 0,
	    0, 1, 0,
	    0, 0, 1,
	    1, 0, 1,
	    1, 1, 1,
	    0, 1, 1], 2, 2, 2, 2, [2, 2, 2, 2, 2, 2])

SCALE (hex, (0.2, 0.05, 0.4))

TRANSLATE (hex, (0.4, 0, 0))

for i in range (0, 4):
  shp = COPY (hex)
  TRANSLATE (shp, (0, i * 0.2, 0))
  BODY (solfec, 'RIGID', shp, bulkmat)

shp = SPHERE ((0.5, -0.5, 0.3), 0.1, 3, 3)

ball = BODY (solfec, 'RIGID', shp, bulkmat)

INITIAL_VELOCITY (ball, (0, 3, 0), (0, 0, 0))

gs = GAUSS_SEIDEL_SOLVER (1E-3, 100)

OUTPUT (solfec, step)

RUN (solfec, gs, 1.0)
