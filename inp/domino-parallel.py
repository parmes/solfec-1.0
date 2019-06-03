# a domino example -- body generation on rank 0 process

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

# For parallel runs it may sometimes be more memory-effective to setup all bodies and boundary
# conditions on rank 0 process (although the unmodified script -- domino.py -- also work fine);
# This approach may be appropriate when the Python part of model genration consumes much resources,
# which are no longer needed during the solution phase
if RANK() == 0:

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

  vb = [0]
  for i in range (0, 4):
    shp = COPY (hex)
    TRANSLATE (shp, (0, i * 0.2, 0))
    b = BODY (solfec, 'RIGID', shp, bulkmat)
    DISPLAY_POINT (b, b.center, ' Label ' + str(i))
    vb.append(b)

  CONTACT_EXCLUDE_BODIES (vb[2], vb[3])

  shp = SPHERE ((0.5, -0.5, 0.3), 0.1, 3, 3)

  ball = BODY (solfec, 'RIGID', shp, bulkmat)

  INITIAL_VELOCITY (ball, (0, 3, 0), (0, 0, 0))

gs = GAUSS_SEIDEL_SOLVER (1E-3, 100)

OUTPUT (solfec, step)

RUN (solfec, gs, 1.0)
