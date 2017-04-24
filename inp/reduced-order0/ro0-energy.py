# import model library routines
import os, sys
dirpath = os.path.dirname(os.path.realpath(__file__))
execfile (dirpath + '/ro0-lib.py')

if VIEWER():
  execfile (dirpath + '/ro0-view.py')
else:
  execfile (dirpath + '/ro0-gen-bases.py')
  import matplotlib.pyplot as plt
  ro0_path_extension ('-energy')
  print '============'
  print 'Energy plots'
  print '============'
  # Plot bar total energy history for:
  # - all schemes
  # - steps 1/64s and 1/256s
  # - damping 0.0
  # - duration 100s
  damp = 0.0
  dura = 100.0
  step_denom = [64, 256]
  for div in step_denom:
    print 'Solving for step 1/%d and damping %g ...' % (div, damp)
    print 'FEM-TL ...',
    sol = ro0_model (1.0/div, damp, 'TL', verbose='%',
		     overwrite=True, runduration=dura)
    print '\bFEM-BC ...',
    sol = ro0_model (1.0/div, damp, 'BC', verbose='%',
		     overwrite=True, runduration=dura)
    print '\bFEM-BC-RO ...',
    sol = ro0_model (1.0/div, damp, 'BC-RO', pod_base,
	verbose='%', overwrite=True, runduration=dura)
    print '\bFEM-BC-MODAL ...',
    sol = ro0_model (1.0/div, damp, 'BC-MODAL', modal_base,
	    verbose='%', overwrite=True, runduration=dura)

    # open in 'READ' mode and retrieve time histories
    print '\bTIME: ', 
    times = ro0_times (ro0_model (1.0/div, damp, 'TL'), progress='ON')
    print '\bFEM-TL: ',
    e_tl = ro0_energy (ro0_model (1.0/div, damp, 'TL'))
    print '\bFEM-BC: ',
    e_bc = ro0_energy (ro0_model (1.0/div, damp, 'BC'))
    print '\bFEM-BC-RO: ',
    e_bc_ro = ro0_energy (ro0_model (1.0/div, damp, 'BC-RO', pod_base))
    print '\bFEM-BC-MODAL: ',
    e_bc_modal = ro0_energy (ro0_model (1.0/div, damp, 'BC-MODAL', modal_base))

    print '\bPlotting energy history h=1/%d, eta=%g over %gs ...' % (div, damp, dura)
    plt.clf ()
    plt.title ('Rotating bar: total energy ($h = 1/%d, \eta = %g$)' % (div, damp))
    plt.plot (times, e_tl, label='TL')
    plt.plot (times, e_bc, label='BC', linestyle='--') # 'None' line style is possible
    plt.plot (times, e_bc_ro, label='BC-RO', linestyle='-')
    plt.plot (times, e_bc_modal, label='BC-MODAL', linestyle='-.')
    plt.xlabel ('Time [s]')
    plt.ylabel ('Energy [J]')
    plt.legend(loc = 'upper right')
    plt.gcf().subplots_adjust(left=0.15)
    plt.savefig ('out/reduced-order0/ro0_energy_h%d_d%g_%gs.png' % (div, damp, dura))
