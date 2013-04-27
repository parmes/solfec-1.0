# splitting tests
import sys

sv = GAUSS_SEIDEL_SOLVER (1E-4, 1000)

sol = SOLFEC ('DYNAMIC', 0.001, 'out/splits')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E6,
		      poisson = 0.3,
		      density = 1E3)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

GRAVITY (sol, (0, 0, -10))

def split_pipe (base, point, normal, topoadj):
  shp = PIPE  (base, (0, 0, 2), 0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])
  point = TRANSLATE (point, base)
  (a, b) = SPLIT (COPY (shp), point, normal, 0, topoadj)
  if topoadj == 'ON':
    if a == None: shp = b
    else:
      if b != None:
	print 'Error in topologically adjacent splitting!'
	sys.exit (1)
      shp = a
    BODY (sol, 'FINITE_ELEMENT', shp, bulk)
  else:
    if a == None or b == None:
      print 'Error in splitting!'
      sys.exit (1)
    else:
      BODY (sol, 'FINITE_ELEMENT', a, bulk)
      BODY (sol, 'FINITE_ELEMENT', b, bulk)

  shp = PIPE  (TRANSLATE (base, (0, 0, -1)), (0, 0, -1), 0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])
  BODY (sol, 'OBSTACLE', shp, bulk)

split_pipe ((0, 0, 0), (1, 0, 0), (0,  1, 0), 'ON')
split_pipe ((4, 0, 0), (0, 1, 0), (1,  0, 0), 'ON')
split_pipe ((8, 0, 0), (1, 1, 0), (1, -1, 0), 'ON')
split_pipe ((0, 4, 0), (1, 0, 0), (0,  1, 0), 'OFF')
split_pipe ((4, 4, 0), (0, 1, 0), (1,  0, 0), 'OFF')
split_pipe ((8, 4, 0), (1, 1, 0), (1, -1, 0), 'OFF')

RUN (sol, sv, 1.0)
