# RIGID_TO_FEM test example 2
# Here we initialise FEM bodies with rigid motion from example 1
# ----------
# 19/02/2016

step = 1E-4

solfec = SOLFEC ('DYNAMIC', step, 'out/devel/rigid-to-fem-2')

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
  b = BODY (solfec, 'FINITE_ELEMENT', shp, bulkmat)
  b.scheme = 'DEF_LIM'
  b.damping = step

shp = SPHERE ((0.5, -0.5, 0.3), 0.1, 3, 3)

ball = BODY (solfec, 'RIGID', shp, bulkmat)

INITIAL_VELOCITY (ball, (0, 3, 0), (0, 0, 0))

gs = GAUSS_SEIDEL_SOLVER (1E-3, 100)

OUTPUT (solfec, step)

# initialise from the end of the rigid simulation
RIGID_TO_FEM ('out/devel/rigid-to-fem-1', 0.5, solfec)

RUN (solfec, gs, 0.5)
