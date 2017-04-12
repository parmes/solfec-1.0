# import model preparation routines
import os, sys
dirpath = os.path.dirname(os.path.realpath(__file__))
execfile (dirpath + '/ro0-prep.py')

# find path to executable
def where(program):
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, program)):
      return path
  return None

if VIEWER():
  # create 'TL','BC','BC-RO', 'BC-MODAL'
  # simulations for viewing purposes
  ro0_model (1.0/256, 0.0, 'TL')
  ro0_model (1.0/256, 0.0, 'BC')

  # locate sofec executable
  path = where('solfec')
  if path == None:
    print 'ERROR: solfec executable not found in PATH!'
    print '       Download and compile solfec; add solfec directory to PATH variable;'
    sys.exit(1)

  # run solfec externally to generate POD and modal bases
  from subprocess import call
  call (['solfec', dirpath + '/ro0-gen-bases.py'])

  # load POD and modal bases
  import pickle
  pod_base = pickle.load(open('out/reduced-order0/pod_base.pickle', 'rb'))
  modal_base = pickle.load(open('out/reduced-order0/modal_base.pickle', 'rb'))

  # create reduced models
  ro0_model (1.0/256, 0.0, 'BC-RO', pod_base)
  ro0_model (1.0/256, 0.0, 'BC-MODAL', modal_base)

else:
  # Initialise POD and modal bases
  print '=============='
  print 'Initialization'
  print '=============='
  sol = ro0_model (1.0/256, 0.0, 'TL', verbose='%', overwrite=True)
  print 'Samping displacements ...',
  rig = RIGID_DISPLACEMENTS (sol.bodies[0])
  defo = COROTATED_DISPLACEMENTS (sol, sol.bodies[0])
  RUN (sol, NEWTON_SOLVER(), 1.0)
  print '\bCalculating POD base ...'
  pod_base = ro0_POD_base (rig, defo, verbose=True)
  print 'Calculating modal base ...'
  modal_base = ro0_modal_base ()

  # Plot bar elongation history for:
  # - all schemes
  # - steps 1/64s and 1/256s
  # - damping 0.0
  print '================'
  print 'Elongation plots'
  print '================'
  print 'Solving for step 1/64 ...'
  print 'FEM-TL ...',
  sol = ro0_model (1.0/64, 0.0, 'TL', verbose='%',
                   overwrite=True, runduration=1.0)
  print '\bFEM-BC ...',
  sol = ro0_model (1.0/64, 0.0, 'BC', verbose='%',
                   overwrite=True, runduration=1.0)
  print '\bFEM-BC-RO ...',
  sol = ro0_model (1.0/64, 0.0, 'BC-RO', pod_base,
      verbose='%', overwrite=True, runduration=1.0)
  print '\bFEM-BC-MODAL ...',
  sol = ro0_model (1.0/64, 0.0, 'BC-MODAL', modal_base,
          verbose='%', overwrite=True, runduration=1.0)

  # open in 'READ' mode and retrieve time histories
  times = ro0_times (ro0_model (1.0/64, 0.0, 'TL'))
  l_tl = ro0_elongation (ro0_model (1.0/64, 0.0, 'TL'))
  l_bc = ro0_elongation (ro0_model (1.0/64, 0.0, 'BC'))
  l_bc_ro = ro0_elongation (ro0_model (1.0/64, 0.0, 'BC-RO', pod_base))
  l_bc_modal = ro0_elongation (ro0_model (1.0/64, 0.0, 'BC-MODAL', modal_base))

  print '\bPlotting elongation history h=1/64, eta=0.0 ...'
  import matplotlib.pyplot as plt
  plt.clf ()
  plt.title ('Rotating bar: elongation $h = 1/64, \eta = 0.0$')
  plt.plot (times, l_tl, label='TL')
  plt.plot (times, l_bc, label='BC', linestyle='--', marker='.', markevery=1) # 'None' line style is possible
  plt.plot (times, l_bc_ro, label='BC-RO', linestyle='-', marker='x', markevery=1)
  plt.plot (times, l_bc_modal, label='BC-MODAL', linestyle='-.', marker='s', markevery=1)
  plt.xlabel ('Time [s]')
  plt.ylabel ('Elongation [m]')
  plt.legend(loc = 'upper right')
  plt.gcf().subplots_adjust(left=0.15)
  plt.savefig ('out/reduced-order0/ro0_elongation_h64_d0.png')
