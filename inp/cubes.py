# stack of cubes example (CONVEX and PSEUDO_RIGID)

from math import pow
import sys

# test kind flag:
# * 1 => parallel growing along one dimension (Gauss-Seidel with empty middle nodes sets)
# * 2 => parallel growing along two dimensions
# * 3 => parallel growing along three dimensions
# * a number > 3 => fixed size model

TEST = 8
KINEM = 'FINITE_ELEMENT'
SOLVER = 'ns'
SAREA = 0.05
step = 0.001
duration =  100 * step
MAKE_TESTS = 0 # make convergence tests

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
	if KINEM != 'RIGID':
	  b.scheme = 'DEF_LIM'
	  b.damping = 1E-4

def create_solver (solver, kinem, sarea, meritval):
  if solver == 'gs':
    sv = GAUSS_SEIDEL_SOLVER (1E-2, 1000, meritval)
  elif solver == 'ts':
    if sarea > 0.0:
	sv = TEST_SOLVER (meritval, 1000, 1000, 100, 0.5, 0.5E-7, 1E-11)
    else:
      if kinem == 'RIGID':
	sv = TEST_SOLVER (meritval, 1000, 1000, 10, 0.4, 0.25E-7, 1E-11)
      elif kinem == 'PSEUDO_RIGID':
	sv = TEST_SOLVER (meritval, 1000, 1000, 10, 0.7, 1E-7, 1E-11)
      else:
        sv = TEST_SOLVER (meritval, 1000, 1000, 10, 0.4, 0.8E-7, 1E-11)
  elif solver == 'ns':
    if sarea > 0.0: sv = NEWTON_SOLVER (meritval, 1000)
    else: sv = NEWTON_SOLVER (meritval, 1000, delta = 1E-7)
  elif solver == 'ps':
    sv = PENALTY_SOLVER ()
  else:
    print 'Unkown solver! => using Gauss-Seidel (1E-4, 1000)'
    sv = GAUSS_SEIDEL_SOLVER (1E-4, 1000)

  return sv

def merhist_callback (sname, sv, mh):
  v = sv.merhist
  if sname != 'ts': t = list (range (0, len (v)))
  else: t = sv.mvhist
  mh.append ((t, v))
  return 1

def create_simulation (nsteps, ncubes, kinem, solver, sarea, frict, meritval):

  outdir = 'out/cubes_' + str (ncubes) + '_' + kinem + '_' + solver + '_' + str (frict)
  solfec = SOLFEC ('DYNAMIC', step, outdir)
  CONTACT_SPARSIFY (solfec, 0.005, sarea)
  SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = frict)
  bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E10, poisson = 0.25, density = 2E3)
  GRAVITY (solfec, (0, 0, -10))

  if kinem == 'rig': kin = 'RIGID'
  elif kinem == 'prb': kin = 'PSEUDO_RIGID'
  else: kin = 'FINITE_ELEMENT'

  # create an obstacle base
  shp = cube (0, 0, -1, ncubes, ncubes, 1, 1, 1)
  BODY (solfec, 'OBSTACLE', shp, bulkmat)

  # create the remaining bricks
  for x in range (ncubes):
    for y in range (ncubes):
      for z in range (ncubes):
	shp = cube (x, y, z, 1, 1, 1, 2, 2)
        b = BODY (solfec, kin, shp, bulkmat)
	if kinem != 'rig':
	  b.scheme = 'DEF_LIM'
	  b.damping = 1E-4

  sv = create_solver (solver, kin, sarea, meritval)
  mh = []
  CALLBACK (solfec, step, (solver, sv, mh), merhist_callback)
  RUN (solfec, sv, step * nsteps)
  return mh



### main module ###
#import rpdb2; rpdb2.start_embedded_debugger('a')

# conv. tests

if not VIEWER() and MAKE_TESTS == 1:

  meritval = 1e-15
  ncubes = 8
  mu = 0.3
  tw1 = create_simulation (1, ncubes, 'rig', 'gs', 0.01, mu, meritval)
  tw2 = create_simulation (1, ncubes, 'rig', 'ts', 0.01, mu, meritval)
  tw3 = create_simulation (1, ncubes, 'prb', 'gs', 0.01, mu, meritval)
  tw4 = create_simulation (1, ncubes, 'prb', 'ts', 0.01, mu, meritval)
  tw5 = create_simulation (1, ncubes, 'fem', 'gs', 0.01, mu, meritval)
  tw6 = create_simulation (1, ncubes, 'fem', 'ts', 0.01, mu, meritval)

  ti1 = create_simulation (1, ncubes, 'rig', 'gs', 0.0, mu, meritval)
  ti2 = create_simulation (1, ncubes, 'rig', 'ts', 0.0, mu, meritval)
  ti3 = create_simulation (1, ncubes, 'prb', 'gs', 0.0, mu, meritval)
  ti4 = create_simulation (1, ncubes, 'prb', 'ts', 0.0, mu, meritval)
  ti5 = create_simulation (1, ncubes, 'fem', 'gs', 0.0, mu, meritval)
  ti6 = create_simulation (1, ncubes, 'fem', 'ts', 0.0, mu, meritval)

  try:
    import matplotlib.pyplot as plt
 
    plt.clf ()
    for th in tw1: plt.plot (th [0], th[1], color = 'b', label = 'GS')
    for th in tw2: plt.plot (th [0], th[1], color = 'r', label = 'PQN1')
    plt.semilogy (10)
    plt.legend(loc = 'upper right')
    plt.title ('Rigid model, well-conditioned, $\mu=0.3$')
    plt.xlabel ('Matric-vector products')
    plt.ylabel ('Merit function g')
    plt.savefig ('out/cubes_8_0.3_well_rig_GS_PQN1.eps')

    plt.clf ()
    for th in tw3: plt.plot (th [0], th[1], color = 'b', label = 'GS')
    for th in tw4: plt.plot (th [0], th[1], color = 'r', label = 'PQN1')
    plt.semilogy (10)
    plt.legend(loc = 'upper right')
    plt.title ('Pseudo-rigid model, well-conditioned, $\mu=0.3$')
    plt.xlabel ('Matric-vector products')
    plt.ylabel ('Merit function g')
    plt.savefig ('out/cubes_8_0.3_well_prb_GS_PQN1.eps')

    plt.clf ()
    for th in tw5: plt.plot (th [0], th[1], color = 'b', label = 'GS')
    for th in tw6: plt.plot (th [0], th[1], color = 'r', label = 'PQN1')
    plt.semilogy (10)
    plt.legend(loc = 'upper right')
    plt.title ('Finite-element model, well-conditioned, $\mu=0.3$')
    plt.xlabel ('Matric-vector products')
    plt.ylabel ('Merit function g')
    plt.savefig ('out/cubes_8_0.3_well_fem_GS_PQN1.eps')

    plt.clf ()
    for th in ti1: plt.plot (th [0], th[1], color = 'b', label = 'GS')
    for th in ti2: plt.plot (th [0], th[1], color = 'r', label = 'PQN1')
    plt.semilogy (10)
    plt.legend(loc = 'upper right')
    plt.title ('Rigid model, ill-conditioned, $\mu=0.3$')
    plt.xlabel ('Matric-vector products')
    plt.ylabel ('Merit function g')
    plt.savefig ('out/cubes_8_0.3_ill_rig_GS_PQN1.eps')

    plt.clf ()
    for th in ti3: plt.plot (th [0], th[1], color = 'b', label = 'GS')
    for th in ti4: plt.plot (th [0], th[1], color = 'r', label = 'PQN1')
    plt.semilogy (10)
    plt.legend(loc = 'upper right')
    plt.title ('Pseudo-rigid model, ill-conditioned, $\mu=0.3$')
    plt.xlabel ('Matric-vector products')
    plt.ylabel ('Merit function g')
    plt.savefig ('out/cubes_8_0.3_ill_prb_GS_PQN1.eps')

    plt.clf ()
    for th in ti5: plt.plot (th [0], th[1], color = 'b', label = 'GS')
    for th in ti6: plt.plot (th [0], th[1], color = 'r', label = 'PQN1')
    plt.semilogy (10)
    plt.legend(loc = 'upper right')
    plt.title ('Finite-element model, ill-conditioned, $\mu=0.3$')
    plt.xlabel ('Matric-vector products')
    plt.ylabel ('Merit function g')
    plt.savefig ('out/cubes_8_0.3_ill_fem_GS_PQN1.eps')

  except ImportError:
    pass # no reaction
 
  sys.exit ()

elif not VIEWER() and MAKE_TESTS == 2:

  meritval = 1e-15
  ncubes = 8
  mu = 0.3
  tw1 = create_simulation (1, ncubes, 'rig', 'gs', 0.01, mu, meritval)
  tw2 = create_simulation (1, ncubes, 'rig', 'ns', 0.01, mu, meritval)
  tw3 = create_simulation (1, ncubes, 'prb', 'gs', 0.01, mu, meritval)
  tw4 = create_simulation (1, ncubes, 'prb', 'ns', 0.01, mu, meritval)
  tw5 = create_simulation (1, ncubes, 'fem', 'gs', 0.01, mu, meritval)
  tw6 = create_simulation (1, ncubes, 'fem', 'ns', 0.01, mu, meritval)

  ti1 = create_simulation (1, ncubes, 'rig', 'gs', 0.0, mu, meritval)
  ti2 = create_simulation (1, ncubes, 'rig', 'ns', 0.0, mu, meritval)
  ti3 = create_simulation (1, ncubes, 'prb', 'gs', 0.0, mu, meritval)
  ti4 = create_simulation (1, ncubes, 'prb', 'ns', 0.0, mu, meritval)
  ti5 = create_simulation (1, ncubes, 'fem', 'gs', 0.0, mu, meritval)
  ti6 = create_simulation (1, ncubes, 'fem', 'ns', 0.0, mu, meritval)

  try:
    import matplotlib.pyplot as plt
 
    plt.clf ()
    for th in tw1: plt.plot (th [0], th[1], color = 'b', label = 'GS')
    for th in tw2: plt.plot (th [0], th[1], color = 'r', label = 'PQN2')
    plt.semilogy (10)
    plt.legend(loc = 'upper right')
    plt.title ('Rigid model, well-conditioned, $\mu=0.3$')
    plt.xlabel ('Matric-vector products')
    plt.ylabel ('Merit function g')
    plt.savefig ('out/cubes_8_0.3_well_rig_GS_PQN2.eps')

    plt.clf ()
    for th in tw3: plt.plot (th [0], th[1], color = 'b', label = 'GS')
    for th in tw4: plt.plot (th [0], th[1], color = 'r', label = 'PQN2')
    plt.semilogy (10)
    plt.legend(loc = 'upper right')
    plt.title ('Pseudo-rigid model, well-conditioned, $\mu=0.3$')
    plt.xlabel ('Matric-vector products')
    plt.ylabel ('Merit function g')
    plt.savefig ('out/cubes_8_0.3_well_prb_GS_PQN2.eps')

    plt.clf ()
    for th in tw5: plt.plot (th [0], th[1], color = 'b', label = 'GS')
    for th in tw6: plt.plot (th [0], th[1], color = 'r', label = 'PQN2')
    plt.semilogy (10)
    plt.legend(loc = 'upper right')
    plt.title ('Finite-element model, well-conditioned, $\mu=0.3$')
    plt.xlabel ('Matric-vector products')
    plt.ylabel ('Merit function g')
    plt.savefig ('out/cubes_8_0.3_well_fem_GS_PQN2.eps')

    plt.clf ()
    for th in ti1: plt.plot (th [0], th[1], color = 'b', label = 'GS')
    for th in ti2: plt.plot (th [0], th[1], color = 'r', label = 'PQN2')
    plt.semilogy (10)
    plt.legend(loc = 'upper right')
    plt.title ('Rigid model, ill-conditioned, $\mu=0.3$')
    plt.xlabel ('Matric-vector products')
    plt.ylabel ('Merit function g')
    plt.savefig ('out/cubes_8_0.3_ill_rig_GS_PQN2.eps')

    plt.clf ()
    for th in ti3: plt.plot (th [0], th[1], color = 'b', label = 'GS')
    for th in ti4: plt.plot (th [0], th[1], color = 'r', label = 'PQN2')
    plt.semilogy (10)
    plt.legend(loc = 'upper right')
    plt.title ('Pseudo-rigid model, ill-conditioned, $\mu=0.3$')
    plt.xlabel ('Matric-vector products')
    plt.ylabel ('Merit function g')
    plt.savefig ('out/cubes_8_0.3_ill_prb_GS_PQN2.eps')

    plt.clf ()
    for th in ti5: plt.plot (th [0], th[1], color = 'b', label = 'GS')
    for th in ti6: plt.plot (th [0], th[1], color = 'r', label = 'PQN2')
    plt.semilogy (10)
    plt.legend(loc = 'upper right')
    plt.title ('Finite-element model, ill-conditioned, $\mu=0.3$')
    plt.xlabel ('Matric-vector products')
    plt.ylabel ('Merit function g')
    plt.savefig ('out/cubes_8_0.3_ill_fem_GS_PQN2.eps')

  except ImportError:
    pass # no reaction

  sys.exit ()


# regular code

outdir = 'out/cubes_' + str (TEST) + '_' + KINEM

solfec = SOLFEC ('DYNAMIC', step, outdir)

CONTACT_SPARSIFY (solfec, 0.005, SAREA)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3, spring = 1E7, dashpot = 1E3)

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E10, poisson = 0.25, density = 2E3)

GRAVITY (solfec, (0, 0, -10))

stack_of_cubes_create (bulkmat, solfec)

sv = create_solver (SOLVER, KINEM, SAREA, 1E-8)

OUTPUT (solfec, 1 * step, 'ON')

MERIT = []

def callback (sv):
  if solfec.time > 0: MERIT.append (sv.merhist)
  return 1

if not VIEWER() and NCPU(solfec) == 1: CALLBACK (solfec, step, sv, callback)

IMBALANCE_TOLERANCE (solfec, 1.1)

RUN (solfec, sv, duration)

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
