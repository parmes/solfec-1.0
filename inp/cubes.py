# stack of cubes example (CONVEX and PSEUDO_RIGID)

from math import pow

# test kind flag:
# * 1 => parallel growing along one dimension (Gauss-Seidel with empty middle nodes sets)
# * 2 => parallel growing along two dimensions
# * 3 => parallel growing along three dimensions
# * a number > 3 => fixed size model

TEST = 8
KINEM = 'PSEUDO_RIGID'
VARIANT = 'FULL'

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

  if TEST == 1:
    N = 10 * NCPU (solfec);
    M = 10
  elif TEST == 2:
    N = 10
    M = int (pow (100 * NCPU (solfec), .5) + 1.)
  elif TEST == 3:
    N = M = int (pow (1000 * NCPU (solfec), 1./3.) + 1.)
  else:
    N = M = TEST

  # create an obstacle base
  shp = cube (0, 0, -1, M, M, 1, 1, 1)
  BODY (solfec, 'OBSTACLE', shp, material)

  # create the remaining bricks
  for x in range (M):
    for y in range (M):
      for z in range (N):
	shp = cube (x, y, z, 1, 1, 1, 2, 2)
        BODY (solfec, KINEM, shp, material)

### main module ###
#import rpdb2; rpdb2.start_embedded_debugger('a')

step = 0.001

solfec = SOLFEC ('DYNAMIC', step, 'out/cubes_' + str (TEST) + '_' + KINEM + '_' + VARIANT)

CONTACT_SPARSIFY (solfec, 0.005)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3)

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E5, poisson = 0.25, density = 1E1)

GRAVITY (solfec, (0, 0, -9.81))

stack_of_cubes_create (bulkmat, solfec)

#gs = GAUSS_SEIDEL_SOLVER (1E-3, 10000, failure = 'CONTINUE', diagsolver = 'PROJECTED_GRADIENT')
#gs.variant = VARIANT

gs = NEWTON_SOLVER ('SMOOTHED_VARIATIONAL', 1E1, 10, 1E-8)
gs.nonmonlength = 5

IMBALANCE_TOLERANCE (solfec, 1.1, 'ON', 2.0)

OUTPUT (solfec, 50 * step, 'FASTLZ')

RUN (solfec, gs, 1000 * step)

if not VIEWER() and solfec.mode == 'READ':

  timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL', 'PARBAL', 'GSINIT', 'GSRUN', 'GSCOM', 'GSMCOM']
  dur = DURATION (solfec)
  th = HISTORY (solfec, timers, dur[0], dur[1])
  total = 0.0

  for i in range (0, len(timers)):
    sum = 0.0
    for tt in th [i+1]: sum += tt
    print timers [i], 'TIME:', sum
    if i < 6: total += sum

  print 'TOTAL TIME:', total
