# balls

step = 1E-2
dura = 3.0

solfec = SOLFEC ('DYNAMIC', step, 'out/balls')

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

pipe = PIPE ((0.5, 0.5, 0), (0, 0, 2), 0.5, 0.1, 1, 8, 1, 1, [1, 2, 3, 4])
BODY (solfec, 'OBSTACLE', pipe, bulkmat)

for i in range (0, 100):
  shp = SPHERE ((0.5 - i % 2 * 0.01, 0.5 + i % 2 * 0.01, 0.1 + i*0.2), 0.1, 3, 3)
  bod = BODY (solfec, 'RIGID', shp, bulkmat)

gs = GAUSS_SEIDEL_SOLVER (1E1, 100, 1E-6)

OUTPUT (solfec, step)

IMBALANCE_TOLERANCE (solfec, 1.3, updatefreq = 10)

RUN (solfec, gs, dura)
