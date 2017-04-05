step = 1E-3
stop = 0.1

sol = SOLFEC ('DYNAMIC', step, 'out/reduced-order0/ro0-reduced')
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
  rig = pickle.load(open('out/reduced-order0/rig.pickle', 'rb'))
  mod = pickle.load(open('out/reduced-order0/mod.pickle', 'rb'))
  val = pickle.load(open('out/reduced-order0/val.pickle', 'rb'))
  base = [x for vec in rig for x in vec] + mod[0] # FIXME --> values in base (1.0s?)
  reduced = ([0,0,0,0,0,0,val[0]], base)
except:
  print 'Any of'
  print 'out/reduced-order0/rig.pickle'
  print 'out/reduced-order0/mod.pickle'
  print 'out/reduced-order0/val.pickle'
  print 'files has not been found'
  print '--> first run solfec examples/reduced-order0/ro0-fem.py'
  print '--> then run python examples/reduced-order0/ro0-modred.py'
  print 'to generate these files!'
  import sys
  sys.exit(0)

msh = PIPE ((0.005, 0.05, 0), (0, 0, 0.1),
            0.01, 0.005, 10, 36, 4, 1, [1]*4)
ROTATE (msh, (0.005, 0.05, 0.05), (0, 1, 0), 90)
bod = BODY (sol, 'FINITE_ELEMENT', msh, mat, form = 'BC-RO', base = reduced)
bod.damping = step

ns = NEWTON_SOLVER ()
OUTPUT (sol, 0.005)

import time
t0 = time.time()
RUN (sol, ns, stop)
t1 = time.time()
print 'Total runtime:', (t1-t0), 'seconds'
