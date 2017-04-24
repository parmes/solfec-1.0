step = 1E-3
stop = 0.1

sol = SOLFEC ('DYNAMIC', step, 'out/reduced-order1/ro1-fem-tl')
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

msh = PIPE ((0.005, 0.05, 0), (0, 0, 0.1),
            0.01, 0.005, 36, 36, 4, 1, [1]*4)
ROTATE (msh, (0.005, 0.05, 0.05), (0, 1, 0), 90)
TRANSLATE (msh, (-0.025, 0, 0))
bod = BODY (sol, 'FINITE_ELEMENT', msh, mat, form = 'TL')
bod.scheme = 'DEF_LIM'
bod.damping = step

ns = NEWTON_SOLVER ()
OUTPUT (sol, 0.0025)
if sol.mode == 'WRITE' and not VIEWER():
  dsp = COROTATED_DISPLACEMENTS (sol, bod)
  rig = RIGID_DISPLACEMENTS (bod)

if sol.mode == 'WRITE':
  import time
  t0 = time.time()
  RUN (sol, ns, stop)
  t1 = time.time()
  print '\bTotal runtime: %.3f seconds' % (t1-t0)
 
  if not VIEWER():
    print 'Saving displacement snapshots ...'
    import pickle
    pickle.dump(dsp, open('out/reduced-order1/dsp.pickle', 'wb')) 
    pickle.dump(rig, open('out/reduced-order1/rig.pickle', 'wb')) 
    print 'Calculating modal decomposition ...'
    MODAL_ANALYSIS (bod, 20, 'out/reduced-order1/modal')

if sol.mode == 'READ' and not VIEWER():
  import pickle
  dur = DURATION (sol)
  th = HISTORY (sol, [(bod, 'KINETIC'), (bod, 'INTERNAL')], dur[0], dur[1])
  pickle.dump(th[0], open('out/reduced-order1/times.pickle', 'wb')) 
  pickle.dump(th[1], open('out/reduced-order1/kin-fem-tl.pickle', 'wb')) 
  pickle.dump(th[2], open('out/reduced-order1/int-fem-tl.pickle', 'wb')) 
