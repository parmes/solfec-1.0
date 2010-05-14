# simple core model
import sys
sys.path.append ('inp/cores/inc')
from simple_core_base import *

# main module
#import rpdb2; rpdb2.start_embedded_debugger('a')

step = 1E-3
stop = 100 * step
scheme = 'IMP'
solver = 'NEWTON'
plotconv = 0

solfec = SOLFEC ('DYNAMIC', step, 'out/cores/simple-small-prb-gs')
GRAVITY (solfec, (0, 0, -10))

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.7)
bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3)

if solver == 'GAUSS_SEIDEL_SOLVER':
  sv = GAUSS_SEIDEL_SOLVER (1E0, 100, 1E-6, failure = 'CONTINUE')
  sv.reverse = 'ON'
else:
  sv = NEWTON_SOLVER ('SMOOTHED_VARIATIONAL', 1E0, 5, 1E-5)
  sv.nonmonlength = 5
  sv.linmaxiter = 5

simple_core_create (0.0003, 0.0002, bulkmat, solfec, 'PSEUDO_RIGID', 'DEF_' + scheme, 'PSEUDO_RIGID', 'DEF_' + scheme, 4, 4, 4)

MERIT = []

def callback (sv):
  MERIT.append (sv.merhist)
  return 1

#OUTPUT (solfec, 0.03)

if not VIEWER() and plotconv == 1: CALLBACK (solfec, step, sv, callback)

RUN (solfec, sv, stop)

if not VIEWER() and solfec.mode == 'WRITE' and plotconv == 1:
  try:
    import matplotlib.pyplot as plt

    for M in MERIT:
      plt.plot (list (range (0, len(M))), M)

    plt.semilogy (10)
    plt.title (solver + ': PSEUDO-RIGID MODEL - ' + scheme)
    plt.xlabel ('Iteration')
    plt.ylabel ('Merit function f(R)')
    plt.savefig ('out/cores/simple-small-prb-gs/simple-small-prb-' + solver + '_' + scheme + '.eps')
 
  except ImportError:
    pass # no reaction

if not VIEWER() and solfec.mode == 'READ':

  timers = ['CONSOL', 'MERIT']
  dur = DURATION (solfec)
  th = HISTORY (solfec, timers, dur[0], dur[1])
  total = 0.0

  try:
    import matplotlib.pyplot as plt

    plt.plot (th[0], th[2])
    plt.semilogy (10)
    plt.title (solver + ': RIGID MODEL')
    plt.xlabel ('Time')
    plt.ylabel ('Merit function f(R)')
    plt.savefig ('out/cores/simple-small-prb-gs/simple-small-prb-' + solver + '-merit.eps')
 
  except ImportError:
    pass # no reaction
