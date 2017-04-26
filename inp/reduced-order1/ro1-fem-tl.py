import os, sys
dirpath = os.path.dirname(os.path.realpath(__file__))
execfile (dirpath + '/ro1-lib.py') # import library

# set up simulation
(sol, bod, stop) = ro1_simulation ('fem-tl')

# collect displacement samples
if sol.mode == 'WRITE' and not VIEWER():
  dsp = COROTATED_DISPLACEMENTS (sol, bod)
  rig = RIGID_DISPLACEMENTS (bod)

# run simulation and save displacemement samples
if sol.mode == 'WRITE':
  import time
  t0 = time.time()
  RUN (sol, NEWTON_SOLVER(), stop)
  t1 = time.time()
  print '\bFEM-TL runtime: %.3f seconds' % (t1-t0)
  if not VIEWER():
    print 'Saving displacement snapshots ...'
    import pickle
    pickle.dump(dsp, open('out/reduced-order1/dsp.pickle', 'wb')) 
    pickle.dump(rig, open('out/reduced-order1/rig.pickle', 'wb')) 

# read and save time histories
if sol.mode == 'READ' and not VIEWER():
  ro1_read_histories (sol, bod, 'fem-tl')
