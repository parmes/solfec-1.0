import matplotlib.pyplot as plt
import pickle, os
d0 = os.path.dirname(os.path.realpath(__file__))
execfile (d0 + '/hs4-globals.py')

# read time histories
try:
  n0t= pickle.load(open('out/hs4-solfec/node0t.pickle', 'rb'))
  n0dy_base = pickle.load(open('out/hs4-solfec/node0dy.pickle', 'rb'))
  n0dy_lwf1 = pickle.load(open('out/hs4-hybrid-lwf-1/node0dy.pickle', 'rb'))
  n0dy_lwf01 = pickle.load(open('out/hs4-hybrid-lwf-0.1/node0dy.pickle', 'rb'))
  n0dy_lwf001 = pickle.load(open('out/hs4-hybrid-lwf-0.01/node0dy.pickle', 'rb'))
  n0dy_base_swp = pickle.load(open('out/hs4-solfec-sweep/node0dy.pickle', 'rb'))
  n0dy_lwf1_swp = pickle.load(open('out/hs4-hybrid-lwf-1-sweep/node0dy.pickle', 'rb'))
  n0dy_lwf01_swp = pickle.load(open('out/hs4-hybrid-lwf-0.1-sweep/node0dy.pickle', 'rb'))
  n0dy_lwf001_swp = pickle.load(open('out/hs4-hybrid-lwf-0.01-sweep/node0dy.pickle', 'rb'))
except:
  print 'A out/hs4-*/[some file].pickle is missing',
  print '--> run hs4-run-all.py file first!'
  import sys
  sys.exit(0)

# plot figures
opath = 'out/hs4_dwell_node0_dy.png'
print '\bPlotting ', opath, '...'
plt.clf ()
plt.title ('ACC dwell at %gHz and mag. %gm/s2\nSolfec-only and hybrid results (lwf=leeway factor)' % (lohi_dwell, amag))
plt.plot (n0t, n0dy_base, label='Solfec-only')
plt.plot (n0t, n0dy_lwf1, label='lfw=1.0', linestyle='--')
plt.plot (n0t, n0dy_lwf01, label='lfw=0.1', linestyle='-.')
plt.plot (n0t, n0dy_lwf001, label='lfw=0.01', linestyle=':')
plt.xlabel ('Time [s]')
plt.ylabel ('DY of node 0 [m]')
plt.legend(loc = 'lower left')
plt.gcf().subplots_adjust(left=0.15)
plt.savefig (opath)

opath = 'out/hs4_sweep_node0_dy.png'
print '\bPlotting ', opath, '...'
plt.clf ()
plt.title ('ACC sweep from %gHz to %gHz and mag. %gm/s2\nSolfec-only and hybrid results (lwf=leeway factor)' % (lofq_sweep, hifq_sweep, amag))
plt.plot (n0t, n0dy_base_swp, label='Solfec-only')
plt.plot (n0t, n0dy_lwf1_swp, label='lfw=1.0', linestyle='--')
plt.plot (n0t, n0dy_lwf01_swp, label='lfw=0.1', linestyle='-.')
plt.plot (n0t, n0dy_lwf001_swp, label='lfw=0.01', linestyle=':')
plt.xlabel ('Time [s]')
plt.ylabel ('DY of node 0 [m]')
plt.legend(loc = 'lower left')
plt.gcf().subplots_adjust(left=0.15)
plt.savefig (opath)
