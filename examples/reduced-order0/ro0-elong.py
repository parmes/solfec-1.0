# import model library routines
import os, sys
dirpath = os.path.dirname(os.path.realpath(__file__))
execfile (dirpath + '/ro0-lib.py')

if VIEWER():
  execfile (dirpath + '/ro0-view.py')
else:
  execfile (dirpath + '/ro0-gen-bases.py')
  import matplotlib.pyplot as plt
  ro0_path_extension ('-elong')
  print '================'
  print 'Elongation plots'
  print '================'
  # Plot bar elongation history for:
  # - all schemes
  # - steps 1/64s and 1/256s
  # - damping 0.0
  step_denom = [64, 256]
  for div in step_denom:
    print 'Solving for step 1/%d and damping 0.0 ...' % div
    print 'FEM-TL ...',
    sol = ro0_model (1.0/div, 0.0, 'TL', verbose='%',
		     overwrite=True, runduration=1.0)
    print '\bFEM-BC ...',
    sol = ro0_model (1.0/div, 0.0, 'BC', verbose='%',
		     overwrite=True, runduration=1.0)
    print '\bFEM-BC-RO ...',
    sol = ro0_model (1.0/div, 0.0, 'BC-RO', pod_base,
	verbose='%', overwrite=True, runduration=1.0)
    print '\bFEM-BC-MODAL ...',
    sol = ro0_model (1.0/div, 0.0, 'BC-MODAL', modal_base,
	    verbose='%', overwrite=True, runduration=1.0)

    # open in 'READ' mode and retrieve time histories
    times = ro0_times (ro0_model (1.0/div, 0.0, 'TL'))
    l_tl = ro0_elongation (ro0_model (1.0/div, 0.0, 'TL'))
    l_bc = ro0_elongation (ro0_model (1.0/div, 0.0, 'BC'))
    l_bc_ro = ro0_elongation (ro0_model (1.0/div, 0.0, 'BC-RO', pod_base))
    l_bc_modal = ro0_elongation (ro0_model (1.0/div, 0.0, 'BC-MODAL', modal_base))

    print '\bPlotting elongation history h=1/%d, eta=0.0 ...' % div
    plt.clf ()
    plt.title ('Rotating bar: elongation ($h = 1/%d, \eta = 0.0$)' % div)
    plt.plot (times, l_tl, label='TL')
    plt.plot (times, l_bc, label='BC', linestyle='--', marker='.', markevery=1) # 'None' line style is possible
    plt.plot (times, l_bc_ro, label='BC-RO', linestyle='-', marker='x', markevery=1)
    plt.plot (times, l_bc_modal, label='BC-MODAL', linestyle='-.', marker='s', markevery=1)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Elongation [m]')
    plt.legend(loc = 'upper right')
    plt.gcf().subplots_adjust(left=0.15)
    plt.savefig ('out/reduced-order0/ro0_elongation_h%d_d0.png' % div)

  # Plot bar elongation history for:
  # - all schemes
  # - steps 1/64s
  # - damping 0.01 and 0.05
  div = 64
  damping = [0.01, 0.05]
  for damp in damping:
    print 'Solving for step 1/%d and damping %g ...' % (div, damp)
    print 'FEM-TL ...',
    sol = ro0_model (1.0/div, damp, 'TL', verbose='%',
		     overwrite=True, runduration=1.0)
    print '\bFEM-BC ...',
    sol = ro0_model (1.0/div, damp, 'BC', verbose='%',
		     overwrite=True, runduration=1.0)
    print '\bFEM-BC-RO ...',
    sol = ro0_model (1.0/div, damp, 'BC-RO', pod_base,
	verbose='%', overwrite=True, runduration=1.0)
    print '\bFEM-BC-MODAL ...',
    sol = ro0_model (1.0/div, damp, 'BC-MODAL', modal_base,
	    verbose='%', overwrite=True, runduration=1.0)

    # open in 'READ' mode and retrieve time histories
    times = ro0_times (ro0_model (1.0/div, damp, 'TL'))
    l_tl = ro0_elongation (ro0_model (1.0/div, damp, 'TL'))
    l_bc = ro0_elongation (ro0_model (1.0/div, damp, 'BC'))
    l_bc_ro = ro0_elongation (ro0_model (1.0/div, damp, 'BC-RO', pod_base))
    l_bc_modal = ro0_elongation (ro0_model (1.0/div, damp, 'BC-MODAL', modal_base))

    print '\bPlotting elongation history h=1/%d, eta=%g ...' % (div, damp)
    plt.clf ()
    plt.title ('Rotating bar: elongation ($h = 1/%d, \eta = %g$)' % (div, damp))
    plt.plot (times, l_tl, label='TL')
    plt.plot (times, l_bc, label='BC', linestyle='--', marker='.', markevery=1) # 'None' line style is possible
    plt.plot (times, l_bc_ro, label='BC-RO', linestyle='-', marker='x', markevery=1)
    plt.plot (times, l_bc_modal, label='BC-MODAL', linestyle='-.', marker='s', markevery=1)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Elongation [m]')
    plt.legend(loc = 'lower right')
    plt.gcf().subplots_adjust(left=0.15)
    plt.savefig ('out/reduced-order0/ro0_elongation_h%d_d%g.png' % (div, damp))
