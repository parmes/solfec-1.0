# a finite-element domino example

step = 1E-3

solfec = SOLFEC ('DYNAMIC', step, 'out/fedomino')

GRAVITY (solfec, (0, 0, -9.81))

sep = 0.06
num = 10
b = 0.0705
w = 0.01926
h = 0.2551 

table = HULL ([0, 0, 0,
               0, 1, 0,
	       1, 1, 0,
	       1, 0, 0,
               0, 0, -0.02,
               0, 1, -0.02,
	       1, 1, -0.02,
	       1, 0, -0.02], 1, 1)

SCALE (table, (b, (w+sep)*num+ 2*h, h))
TRANSLATE (table, (0, -h, 0))

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.2, restitution = 0)

BODY (solfec, 'OBSTACLE', table, bulkmat)

hex = HEX ([0, 0, 0,
	    1, 0, 0,
	    1, 1, 0,
	    0, 1, 0,
	    0, 0, 1,
	    1, 0, 1,
	    1, 1, 1,
	    0, 1, 1], 4, 3, 10, 2, [2, 2, 2, 2, 2, 2])

SCALE (hex, (b, w, h))

for i in range (0, num):
  shp = COPY (hex)
  TRANSLATE (shp, (0, i * sep, 0))
  bod = BODY (solfec, 'FINITE_ELEMENT', shp, bulkmat, form = 'BC')
  bod.scheme = 'DEF_LIM'
  bod.damping = step
  if i == 0:
    con = SET_VELOCITY (bod, (b/2.,0,h), (0, 1, 0), 1)

sv = GAUSS_SEIDEL_SOLVER  (1E-3, 200)

RUN (solfec, sv, 0.025)

DELETE (solfec, con)

RUN (solfec, sv, 1.0)
