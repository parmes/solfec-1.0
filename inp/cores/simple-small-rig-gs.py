# simple core model
import sys
sys.path.append ('inp/cores/inc')
from simple_core_base import *

# main module
#import rpdb2; rpdb2.start_embedded_debugger('a')

step = 1E-3
stop = 5
solver = 'GAUSS_SEIDEL_SOLVER'
plotconv = 0

solfec = SOLFEC ('DYNAMIC', step, 'out/cores/simple-small-rig-gs')
GRAVITY (solfec, (0, 0, -10))

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.7)
bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3)

if solver == 'GAUSS_SEIDEL_SOLVER':
  sv = GAUSS_SEIDEL_SOLVER (1E0, 8, 1E-6, failure = 'CONTINUE')
  sv.reverse = 'ON'
else:
  #sv = NEWTON_SOLVER ('SMOOTHED_VARIATIONAL', 1E0, 50, 1E-6)
  #sv.nonmonlength = 5
  sv = HYBRID_SOLVER ()

simple_core_create (0.0003, 0.0002, bulkmat, solfec, 'RIGID', 'DEFAULT', 'RIGID', 'DEFAULT', 4, 4, 4)

MERIT = []

def callback (sv):
  MERIT.append (sv.merhist)
  return 1

OUTPUT (solfec, 0.01)

if not VIEWER() and plotconv == 1: CALLBACK (solfec, step, sv, callback)

RUN (solfec, sv, stop)

if not VIEWER() and solfec.mode == 'WRITE' and plotconv == 1:
  try:
    import matplotlib.pyplot as plt

    for M in MERIT:
      plt.plot (list (range (0, len(M))), M)

    plt.semilogy (10)
    plt.title (solver + ': RIGID MODEL')
    plt.xlabel ('Iteration')
    plt.ylabel ('Merit function f(R)')
    plt.savefig ('out/cores/simple-small-rig-gs/simple-small-rig-' + solver + '.eps')
 
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
    plt.savefig ('out/cores/simple-small-rig-gs/simple-small-rig-' + solver + '-merit.eps')
 
  except ImportError:
    pass # no reaction
