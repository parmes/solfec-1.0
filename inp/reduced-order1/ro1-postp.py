import matplotlib.pyplot as plt
import pickle

times = pickle.load(open('out/reduced-order1/times.pickle', 'rb'))
kin_tl = pickle.load(open('out/reduced-order1/kin-fem-tl.pickle', 'rb'))
kin_bc = pickle.load(open('out/reduced-order1/kin-fem-bc.pickle', 'rb'))
kin_bc_ro = pickle.load(open('out/reduced-order1/kin-reduced.pickle', 'rb'))
kin_bc_modal = pickle.load(open('out/reduced-order1/kin-modal.pickle', 'rb'))
int_tl = pickle.load(open('out/reduced-order1/int-fem-tl.pickle', 'rb'))
int_bc = pickle.load(open('out/reduced-order1/int-fem-bc.pickle', 'rb'))
int_bc_ro = pickle.load(open('out/reduced-order1/int-reduced.pickle', 'rb'))
int_bc_modal = pickle.load(open('out/reduced-order1/int-modal.pickle', 'rb'))

opath = 'out/reduced-order1/ro1_kin_energy.png'
print '\bPlotting ', opath, '...'
plt.clf ()
plt.title ('Pipe impact: kinetic energy')
plt.plot (times, kin_tl, label='TL')
plt.plot (times, kin_bc, label='BC', linestyle='--') # 'None' line style is possible
plt.plot (times, kin_bc_ro, label='BC-RO', linestyle='-')
plt.plot (times, kin_bc_modal, label='BC-MODAL', linestyle='-.')
plt.xlabel ('Time [s]')
plt.ylabel ('Energy [J]')
plt.legend(loc = 'lower right')
plt.gcf().subplots_adjust(left=0.15)
plt.savefig (opath)

opath = 'out/reduced-order1/ro1_int_energy.png'
print '\bPlotting ', opath, '...'
plt.clf ()
plt.title ('Pipe impact: internal energy')
plt.plot (times, int_tl, label='TL')
plt.plot (times, int_bc, label='BC', linestyle='--') # 'None' line style is possible
plt.plot (times, int_bc_ro, label='BC-RO', linestyle='-')
plt.plot (times, int_bc_modal, label='BC-MODAL', linestyle='-.')
plt.xlabel ('Time [s]')
plt.ylabel ('Energy [J]')
plt.legend(loc = 'upper left')
plt.gcf().subplots_adjust(left=0.15)
plt.savefig (opath)
