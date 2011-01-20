# stack of cubes example (CONVEX and PSEUDO_RIGID)

from math import pow

# test kind flag:
# * 1 => parallel growing along one dimension (Gauss-Seidel with empty middle nodes sets)
# * 2 => parallel growing along two dimensions
# * 3 => parallel growing along three dimensions
# * a number > 3 => fixed size model

TEST = 8
KINEM = 'RIGID'
SOLVER = 'nt'
SAREA = 0.05

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
        b = BODY (solfec, KINEM, shp, material)
	if KINEM != 'RIGID': b.scheme = 'DEF_LIM2'

### main module ###
#import rpdb2; rpdb2.start_embedded_debugger('a')

step = 0.001

outdir = 'out/cubes_' + str (TEST) + '_' + KINEM

solfec = SOLFEC ('DYNAMIC', step, outdir)

CONTACT_SPARSIFY (solfec, 0.005, SAREA)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3)

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E9, poisson = 0.25, density = 1E3)

GRAVITY (solfec, (0, 0, -10))

stack_of_cubes_create (bulkmat, solfec)

if SOLVER == 'gs':
  sv = GAUSS_SEIDEL_SOLVER (1E-2, 1000, 1E-7)
else:
  sv = NEWTON_SOLVER (1E-7, 1000)

IMBALANCE_TOLERANCE (solfec, 1.1, 'ON', 2.0)

OUTPUT (solfec, 1 * step, 'ON')

MERIT = []

def callback (sv):
  if solfec.time > 0: MERIT.append (sv.merhist)
  return 1

if not VIEWER() and NCPU(solfec) == 1: CALLBACK (solfec, step, sv, callback)

RUN (solfec, sv, 1 * step)

# Post-processing
def write_data (t, v, path):
  out = open (path, 'w')
  n = min (len(t), len(v))
  for i in range (n):
    out.write ('%g\t%g\n' % (t[i], v[i]))
  out.close ()

if not VIEWER() and NCPU(solfec) == 1 and solfec.mode == 'WRITE':
  try:
    import matplotlib.pyplot as plt

    i = 0
    for M in MERIT:
      t = list (range (0, len(M)))
      plt.plot (t, M)
      if i == 0:
	if SAREA > 0.0: kind = 'welld'
	else: kind = 'overd'
	path = 'out/cubes_' + SOLVER + '_' + kind + '.dat'
	write_data (t, M, path)
      i = i + 1

    plt.semilogy (10)
    plt.xlabel ('Iteration')
    plt.ylabel ('Merit function f')
    plt.savefig (outdir + '/cubes.eps')
 
  except ImportError:
    pass # no reaction

if not VIEWER() and solfec.mode == 'READ':

  timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL', 'PARBAL']
  dur = DURATION (solfec)
  th = HISTORY (solfec, timers, dur[0], dur[1])
  total = 0.0

  for i in range (0, len(timers)):
    sum = 0.0
    for tt in th [i+1]: sum += tt
    print timers [i], 'TIME:', sum
    if i < 6: total += sum

  print 'TOTAL TIME:', total
