# import model library routines
import os, sys
dirpath = os.path.dirname(os.path.realpath(__file__))
execfile (dirpath + '/ro0-lib.py')

if VIEWER():
  execfile (dirpath + '/ro0-view.py')
else:
  execfile (dirpath + '/ro0-gen-bases.py')
  import matplotlib.pyplot as plt
  ro0_path_extension ('-convtest')
  print '================='
  print 'Converfence plots'
  print '================='
  TL0 = ro0_convtest ('TL', 200E4, 0.000, 16, 7, 14, 4)
  BC0 = ro0_convtest ('BC', 200E4, 0.000, 16, 7, 14, 4)
  RO0 = ro0_convtest ('BC-RO', 200E4, 0.000, 16, 7, 14, 4, pod_base)
  MO0 = ro0_convtest ('BC-MODAL', 200E4, 0.000, 16, 7, 14, 4, modal_base)

  TL1 = ro0_convtest ('TL', 200E9, 1E-6, 16, 7, 14, 4)
  BC1 = ro0_convtest ('BC', 200E9, 1E-6, 16, 7, 14, 4)
  RO1 = ro0_convtest ('BC-RO', 200E9, 1E-6, 16, 7, 14, 4, pod_base)
  MO1 = ro0_convtest ('BC-MODAL', 200E9, 1E-6, 16, 7, 14, 4, modal_base)

  plt.clf ()
  plt.title ('Rotating bar: convergence rate (E = 200E4, $\eta=0$)')
  plt.loglog (TL0[0], TL0[1], label='TL', ls = '-.', lw=3)
  plt.loglog (BC0[0], BC0[1], label='BC')
  plt.loglog (RO0[0], RO0[1], label='BC-RO', ls = '-.', marker = 's')
  plt.loglog (MO0[0], MO0[1], label='BC-MODAL', ls = '--', marker = 'o')
  plt.xlabel ('Time step $h$ [s]')
  plt.ylabel ('Solution distance $\|q_{ref} - q_h\|$ [m]')
  plt.legend(loc = 'lower right')
  plt.savefig ('out/reduced-order0/ro0_convtest_undamped_E200E4.png')

  plt.clf ()
  plt.title ('Rotating bar: convergence rate (E = 200E9, $\eta=1/10^6$)')
  plt.loglog (TL1[0], TL1[1], label='TL', ls = '-.', lw=3)
  plt.loglog (BC1[0], BC1[1], label='BC')
  plt.loglog (RO1[0], RO1[1], label='BC-RO', ls = '-.', marker = 's')
  plt.loglog (MO1[0], MO1[1], label='BC-MODAL', ls = '--', marker = 'o')
  plt.xlabel ('Time step $h$ [s]')
  plt.ylabel ('Solution distance $\|q_{ref} - q_h\|$ [m]')
  plt.legend(loc = 'lower right')
  plt.savefig ('out/reduced-order0/ro0_convtest_damped_E200E9.png')
