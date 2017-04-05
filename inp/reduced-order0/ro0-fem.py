step = 1E-3
stop = 0.1

sol = SOLFEC ('DYNAMIC', step, 'out/reduced-order0/ro0-fem')
GRAVITY (sol, (0, 0, -10))
mat = BULK_MATERIAL (sol, model = 'KIRCHHOFF',
       young = 1E6, poisson = 0.25, density = 1E3)
SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.1)

nodes = [0, 0, 0,   1, 0, 0,   1, 1, 0,   0, 1, 0,
	 0, 0, 1,   1, 0, 1,   1, 1, 1,   0, 1, 1]
msh = HEX (nodes, 1, 1, 1, 0, [0]*6)
SCALE (msh, (0.01, 0.1, 0.01))
BODY (sol, 'OBSTACLE', msh, mat)

msh = PIPE ((0.005, 0.05, 0), (0, 0, 0.1),
            0.01, 0.005, 10, 36, 4, 1, [1]*4)
ROTATE (msh, (0.005, 0.05, 0.05), (0, 1, 0), 90)
bod = BODY (sol, 'FINITE_ELEMENT', msh, mat, form = 'BC')
bod.scheme = 'DEF_LIM'
bod.damping = step

ns = NEWTON_SOLVER ()
OUTPUT (sol, 0.005)
if sol.mode == 'WRITE' and not VIEWER():
  dsp = COROTATED_DISPLACEMENTS (sol, bod)
  rig = RIGID_DISPLACEMENTS (bod)

import time
t0 = time.time()
RUN (sol, ns, stop)
t1 = time.time()
print 'Total runtime:', (t1-t0), 'seconds'

if sol.mode == 'WRITE' and not VIEWER():
  print 'Saving displacement snapshots ...'
  import pickle
  pickle.dump(dsp, open('out/reduced-order0/dsp.pickle', 'wb')) 
  pickle.dump(rig, open('out/reduced-order0/rig.pickle', 'wb')) 
  print 'Calculating modal decomposition ...'
  MODAL_ANALYSIS (bod, 20, 'out/reduced-order0/modal')
