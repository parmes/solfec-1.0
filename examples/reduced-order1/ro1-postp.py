import matplotlib.pyplot as plt
import pickle

# read time histories
try:
  times = pickle.load(open('out/reduced-order1/times.pickle', 'rb'))
  kin_tl = pickle.load(open('out/reduced-order1/kin-fem-tl.pickle', 'rb'))
  kin_bc = pickle.load(open('out/reduced-order1/kin-fem-bc.pickle', 'rb'))
  kin_bc_ro = pickle.load(open('out/reduced-order1/kin-reduced.pickle', 'rb'))
  int_tl = pickle.load(open('out/reduced-order1/int-fem-tl.pickle', 'rb'))
  int_bc = pickle.load(open('out/reduced-order1/int-fem-bc.pickle', 'rb'))
  int_bc_ro = pickle.load(open('out/reduced-order1/int-reduced.pickle', 'rb'))
  odm_tl = pickle.load(open('out/reduced-order1/odm-fem-tl.pickle', 'rb'))
  odm_bc = pickle.load(open('out/reduced-order1/odm-fem-bc.pickle', 'rb'))
  odm_bc_ro = pickle.load(open('out/reduced-order1/odm-reduced.pickle', 'rb'))
  idm_tl = pickle.load(open('out/reduced-order1/idm-fem-tl.pickle', 'rb'))
  idm_bc = pickle.load(open('out/reduced-order1/idm-fem-bc.pickle', 'rb'))
  idm_bc_ro = pickle.load(open('out/reduced-order1/idm-reduced.pickle', 'rb'))
  rot_tl = pickle.load(open('out/reduced-order1/rot-fem-tl.pickle', 'rb'))
  rot_bc = pickle.load(open('out/reduced-order1/rot-fem-bc.pickle', 'rb'))
  rot_bc_ro = pickle.load(open('out/reduced-order1/rot-reduced.pickle', 'rb'))
  p4z_tl = pickle.load(open('out/reduced-order1/p4z-fem-tl.pickle', 'rb'))
  p4z_bc = pickle.load(open('out/reduced-order1/p4z-fem-bc.pickle', 'rb'))
  p4z_bc_ro = pickle.load(open('out/reduced-order1/p4z-reduced.pickle', 'rb'))
  p5z_tl = pickle.load(open('out/reduced-order1/p5z-fem-tl.pickle', 'rb'))
  p5z_bc = pickle.load(open('out/reduced-order1/p5z-fem-bc.pickle', 'rb'))
  p5z_bc_ro = pickle.load(open('out/reduced-order1/p5z-reduced.pickle', 'rb'))
  p0mises_tl = pickle.load(open('out/reduced-order1/p0mises-fem-tl.pickle', 'rb'))
  p0mises_bc = pickle.load(open('out/reduced-order1/p0mises-fem-bc.pickle', 'rb'))
  p0mises_bc_ro = pickle.load(open('out/reduced-order1/p0mises-reduced.pickle', 'rb'))
  p2mises_tl = pickle.load(open('out/reduced-order1/p2mises-fem-tl.pickle', 'rb'))
  p2mises_bc = pickle.load(open('out/reduced-order1/p2mises-fem-bc.pickle', 'rb'))
  p2mises_bc_ro = pickle.load(open('out/reduced-order1/p2mises-reduced.pickle', 'rb'))
except:
  print 'A out/reduced-order1/[some file].pickle is missing',
  print '--> run ro0-run-all.py input file first!'
  import sys
  sys.exit(0)

# plots time window
t0 = 0.05
t1 = 0.1

# plot figures
opath = 'out/reduced-order1/ro1_kin_energy.png'
print '\bPlotting ', opath, '...'
plt.clf ()
plt.title ('Pipe impact: kinetic energy')
plt.plot (times, kin_tl, label='TL')
plt.plot (times, kin_bc, label='BC', linestyle='--') # 'None' line style is possible
plt.plot (times, kin_bc_ro, label='BC-RO', linestyle='-')
plt.xlim ((t0, t1))
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
plt.xlim ((t0, t1))
plt.xlabel ('Time [s]')
plt.ylabel ('Energy [J]')
plt.legend(loc = 'upper left')
plt.gcf().subplots_adjust(left=0.15)
plt.savefig (opath)

opath = 'out/reduced-order1/ro1_outer_diameter.png'
print '\bPlotting ', opath, '...'
plt.clf ()
plt.title ('Pipe impact: outer diameter in the contact area')
plt.plot (times, odm_tl, label='TL')
plt.plot (times, odm_bc, label='BC', linestyle='--') # 'None' line style is possible
plt.plot (times, odm_bc_ro, label='BC-RO', linestyle='-')
plt.xlim ((t0, t1))
plt.xlabel ('Time [s]')
plt.ylabel ('Length [m]')
plt.legend(loc = 'upper left')
plt.gcf().subplots_adjust(left=0.15)
plt.savefig (opath)

opath = 'out/reduced-order1/ro1_inner_diameter.png'
print '\bPlotting ', opath, '...'
plt.clf ()
plt.title ('Pipe impact: inner diameter in the contact area')
plt.plot (times, idm_tl, label='TL')
plt.plot (times, idm_bc, label='BC', linestyle='--') # 'None' line style is possible
plt.plot (times, idm_bc_ro, label='BC-RO', linestyle='-')
plt.xlim ((t0, t1))
plt.xlabel ('Time [s]')
plt.ylabel ('Length [m]')
plt.legend(loc = 'upper left')
plt.gcf().subplots_adjust(left=0.15)
plt.savefig (opath)

opath = 'out/reduced-order1/ro1_rotation.png'
print '\bPlotting ', opath, '...'
plt.clf ()
plt.title ('Pipe impact: rotation angle')
plt.plot (times, rot_tl, label='TL')
plt.plot (times, rot_bc, label='BC', linestyle='--') # 'None' line style is possible
plt.plot (times, rot_bc_ro, label='BC-RO', linestyle='-')
plt.xlim ((t0, t1))
plt.xlabel ('Time [s]')
plt.ylabel ('Angle [rad]')
plt.legend(loc = 'upper left')
plt.gcf().subplots_adjust(left=0.15)
plt.savefig (opath)

opath = 'out/reduced-order1/ro1_p4z_p5z.png'
print '\bPlotting ', opath, '...'
plt.clf ()
plt.title ('Pipe impact: p4z and p5z')
plt.plot (times, p4z_tl, label='p4z-TL')
plt.plot (times, p4z_bc, label='p4z-BC', linestyle='--') # 'None' line style is possible
plt.plot (times, p4z_bc_ro, label='p4z-BC-RO', linestyle='-')
plt.plot (times, p5z_tl, label='p5z-TL')
plt.plot (times, p5z_bc, label='p5z-BC', linestyle='--') # 'None' line style is possible
plt.plot (times, p5z_bc_ro, label='p5z-BC-RO', linestyle='-')
plt.xlim ((t0, t1))
plt.xlabel ('Time [s]')
plt.ylabel ('Position [m]')
plt.legend(loc = 'lower left')
plt.gcf().subplots_adjust(left=0.15)
plt.savefig (opath)

opath = 'out/reduced-order1/ro1_p0mises.png'
print '\bPlotting ', opath, '...'
plt.clf ()
plt.title ('Pipe impact: p0 MISES stress')
plt.plot (times, p0mises_tl, label='TL')
plt.plot (times, p0mises_bc, label='BC', linestyle='--') # 'None' line style is possible
plt.plot (times, p0mises_bc_ro, label='BC-RO', linestyle='-')
plt.xlim ((t0, t1))
plt.xlabel ('Time [s]')
plt.ylabel ('Stress [MPa]')
plt.legend(loc = 'upper left')
plt.gcf().subplots_adjust(left=0.15)
plt.savefig (opath)

opath = 'out/reduced-order1/ro1_p2mises.png'
print '\bPlotting ', opath, '...'
plt.clf ()
plt.title ('Pipe impact: p2 MISES stress')
plt.plot (times, p2mises_tl, label='TL')
plt.plot (times, p2mises_bc, label='BC', linestyle='--') # 'None' line style is possible
plt.plot (times, p2mises_bc_ro, label='BC-RO', linestyle='-')
plt.xlim ((t0, t1))
plt.xlabel ('Time [s]')
plt.ylabel ('Stress [MPa]')
plt.legend(loc = 'upper left')
plt.gcf().subplots_adjust(left=0.15)
plt.savefig (opath)
