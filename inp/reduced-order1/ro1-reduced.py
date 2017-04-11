step = 1E-3
stop = 0.1

sol = SOLFEC ('DYNAMIC', step, 'out/reduced-order1/ro1-reduced')
if 'percentage' in locals(): sol.verbose = '%'
GRAVITY (sol, (0, 0, -10))
mat = BULK_MATERIAL (sol, model = 'KIRCHHOFF',
       young = 1E6, poisson = 0.25, density = 1E3)
SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.1)

nodes = [0, 0, 0,   1, 0, 0,   1, 1, 0,   0, 1, 0,
	 0, 0, 1,   1, 0, 1,   1, 1, 1,   0, 1, 1]
msh = HEX (nodes, 1, 1, 1, 0, [0]*6)
SCALE (msh, (0.01, 0.1, 0.01))
BODY (sol, 'OBSTACLE', msh, mat)

try:
  import pickle
  mod = pickle.load(open('out/reduced-order1/mod.pickle', 'rb'))
  val = pickle.load(open('out/reduced-order1/val.pickle', 'rb'))
  base = [x for vec in mod for x in vec]
  reduced = (val[0:len(mod)], base)
except:
  print 'Any of'
  print 'out/reduced-order1/rig.pickle'
  print 'out/reduced-order1/mod.pickle'
  print 'out/reduced-order1/val.pickle'
  print 'files has not been found'
  print '--> first run solfec examples/reduced-order1/ro1-fem.py'
  print '--> then run python examples/reduced-order1/ro1-modred.py'
  print 'to generate these files!'
  import sys
  sys.exit(0)

msh = PIPE ((0.005, 0.05, 0), (0, 0, 0.1),
            0.01, 0.005, 36, 36, 4, 1, [1]*4)
ROTATE (msh, (0.005, 0.05, 0.05), (0, 1, 0), 90)
TRANSLATE (msh, (-0.025, 0, 0))
bod = BODY (sol, 'FINITE_ELEMENT', msh, mat, form = 'BC-RO', base = reduced)
bod.damping = step

ns = NEWTON_SOLVER ()
OUTPUT (sol, 0.0025)

import time
t0 = time.time()
RUN (sol, ns, stop)
t1 = time.time()
print 'Total runtime: %.3f seconds' % (t1-t0)
