import os, sys
dirpath = os.path.dirname(os.path.realpath(__file__))
execfile (dirpath + '/ro1-lib.py') # import library

# read simulations
(sol_tl, bod_tl, stop_tl) = ro1_simulation ('fem-tl')
if sol_tl.mode == 'WRITE': print 'Results all TL simulation are not present'

(sol_bc, bod_bc, stop_bc) = ro1_simulation ('fem-bc')
if sol_bc.mode == 'WRITE': print 'Results all BC simulation are not present'

sol_ro_mode = 'WRITE'
if sol_tl.mode == 'READ':
  (sol_ro, bod_ro, stop_ro) = ro1_simulation ('reduced')
  sol_ro_mode = sol_ro.mode
  if sol_ro.mode == 'WRITE': print 'Results all BC-RO simulation are not present'

if sol_tl.mode == 'WRITE' or sol_bc.mode == 'WRITE' or sol_ro_mode == 'WRITE':
  print 'Run first --> solfec exaples/reduced-order1/ro1-run-all.py'
  print 'Run next --> solfec -v exaples/reduced-order1/ro1-view.py'
  sys.exit(0)
