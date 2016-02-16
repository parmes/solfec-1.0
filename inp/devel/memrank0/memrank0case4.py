# This is a modified domino example, to exemplify possible sources of issues
# with excessive use of memory on rank 0 during MPI runs with FEM meshes
# -----------
# 16 Feb 2016
# -----------
# case 4: many meshes created and assigned on all ranks;
#         along this auxiliary larger meshe created on all ranks, within the loop scope, and saved
#         into a list stored outside of the loop scope; this is done only on rank 0
#
# result: uniform use of memory on all ranks (approx 300 GB) but rank 0 (approx 1 GB);
#         suggestion to Steve: perhaps there is a similar construct in your input decks?

step = 1E-3

solfec = SOLFEC ('DYNAMIC', step, 'out/devel/memrank0case4')

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

SCALE (table, (1, 200, 1))
BODY (solfec, 'OBSTACLE', table, bulkmat)

hex = HEX ([0, 0, 0,
	    1, 0, 0,
	    1, 1, 0,
	    0, 1, 0,
	    0, 0, 1,
	    1, 0, 1,
	    1, 1, 1,
	    0, 1, 1], 4, 4, 4, 2, [2, 2, 2, 2, 2, 2])

SCALE (hex, (0.2, 0.05, 0.4))
TRANSLATE (hex, (0.4, 0, 0))

temp = []
for i in range (0, 1000):
  shp = COPY (hex)

  shp2 = HEX ([0, 0, 0,
	      1, 0, 0,
	      1, 1, 0,
	      0, 1, 0,
	      0, 0, 1,
	      1, 0, 1,
	      1, 1, 1,
	      0, 1, 1], 16, 16, 16, 2, [2, 2, 2, 2, 2, 2])

  if RANK() == 0:
    temp.append(shp2)

  TRANSLATE (shp, (0, i * 0.2, 0))
  b = BODY (solfec, 'FINITE_ELEMENT', shp, bulkmat)
  b.scheme = 'DEF_LIM'

shp = SPHERE ((0.5, -0.5, 0.3), 0.1, 3, 3)

ball = BODY (solfec, 'RIGID', shp, bulkmat)

INITIAL_VELOCITY (ball, (0, 3, 0), (0, 0, 0))

gs = GAUSS_SEIDEL_SOLVER (1E-3, 100)

OUTPUT (solfec, step)

RUN (solfec, gs, 10*step)
