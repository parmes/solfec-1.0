# balls

step = 1E-3

solfec = SOLFEC ('DYNAMIC', step, 'out/domino')

GRAVITY (solfec, (0, 0, -10))

table = HULL ([0, 0, 0,
               0, 1, 0,
	       1, 1, 0,
	       1, 0, 0,
               0, 0, -0.1,
               0, 1, -0.1,
	       1, 1, -0.1,
	       1, 0, -0.1], 1, 1)

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5, restitution = 0)

BODY (solfec, 'OBSTACLE', table, bulkmat)

for i in range (0, 2):
  shp = SPHERE ((0.5, 0.5 + i * 0.01, 0.1 + i*0.2), 0.1, 3, 3)
  bod = BODY (solfec, 'RIGID', shp, bulkmat)

gs = GAUSS_SEIDEL_SOLVER (1E1, 1000, 1E-6)

OUTPUT (solfec, step)

RUN (solfec, gs, 0.497)
