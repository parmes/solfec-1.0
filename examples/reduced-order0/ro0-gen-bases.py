# import model library routines
import os
dirpath = os.path.dirname(os.path.realpath(__file__))
execfile (dirpath + '/ro0-lib.py')

# Initialise POD and modal bases
try:
  import pickle
  pod_base = pickle.load(open('out/reduced-order0/pod_base.pickle', 'rb'))
  modal_base = pickle.load(open('out/reduced-order0/modal_base.pickle', 'rb'))
except:
  ro0_path_extension ('-init')
  print '=============='
  print 'Initialization'
  print '=============='
  sol = ro0_model (1.0/256, 0.0, 'TL', verbose='%\n', overwrite=True)
  print 'Samping displacements ...',
  rig = RIGID_DISPLACEMENTS (sol.bodies[0])
  defo = COROTATED_DISPLACEMENTS (sol, sol.bodies[0])
  RUN (sol, NEWTON_SOLVER(), 1.0)
  print '\bCalculating POD base ...'
  pod_base = ro0_POD_base (rig, defo, verbose=True)
  print 'Calculating modal base ...'
  modal_base = ro0_modal_base ()

  # Save POD and modal bases
  import pickle
  pickle.dump(pod_base, open('out/reduced-order0/pod_base.pickle', 'wb')) 
  pickle.dump(modal_base, open('out/reduced-order0/modal_base.pickle', 'wb')) 
