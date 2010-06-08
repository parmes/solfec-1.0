# stack of cubes

KINEM = 'PSEUDO_RIGID'
M = 2
N = 1

def cube (x, y, z, a, b, c, sur, vol):

  nodes = [0, 0, 0,
	   a, 0, 0,
	   a, b, 0,
	   0, b, 0,
	   0, 0, c,
	   a, 0, c,
	   a, b, c,
	   0, b, c]

  shp = HEX (nodes, 1, 1, 1, vol, [sur, sur, sur, sur, sur, sur])

  TRANSLATE (shp, (x, y, z))

  return shp

def stack_of_cubes_create (material, solfec):

  # create an obstacle base
  shp = cube (0, 0, -1, M, M, 1, 1, 1)
  BODY (solfec, 'OBSTACLE', shp, material)

  # create the remaining bricks
  for x in range (M):
    for y in range (M):
      for z in range (N):
	#shp = cube (x, y, z, 1, 1, 1, 2, 2)
	shp = SPHERE ((x+0.5, y+0.5, z+0.5), 0.5, 2, 2)
        BODY (solfec, KINEM, shp, material)

### main module ###
#import rpdb2; rpdb2.start_embedded_debugger('a')

step = 0.001

solfec = SOLFEC ('DYNAMIC', step, 'out/tests/locdyn')

CONTACT_SPARSIFY (solfec, 0.005, 0.001)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3)

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E5, poisson = 0.25, density = 1E1)

GRAVITY (solfec, (0, 0, -10))

stack_of_cubes_create (bulkmat, solfec)

sv = GAUSS_SEIDEL_SOLVER (1E-3, 10)

IMBALANCE_TOLERANCE (solfec, 1.1, 'ON', 2.0)

OUTPUT (solfec, step)

RUN (solfec, sv, step)

LOCDYN_DUMP (solfec, 'out/tests/locdyn-' + str (NCPU (solfec)))
