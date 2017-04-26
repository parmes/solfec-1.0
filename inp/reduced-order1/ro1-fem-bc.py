import os, sys
dirpath = os.path.dirname(os.path.realpath(__file__))
execfile (dirpath + '/ro1-lib.py') # import library

# set up simulation
(sol, bod, stop) = ro1_simulation ('fem-bc')

# run simulation
if sol.mode == 'WRITE':
  import time
  t0 = time.time()
  RUN (sol, NEWTON_SOLVER(), stop)
  t1 = time.time()
  print '\bFEM-BC runtime: %.3f seconds' % (t1-t0)

# read and save time histories
if sol.mode == 'READ' and not VIEWER():
  ro1_read_histories (sol, bod, 'fem-bc')
