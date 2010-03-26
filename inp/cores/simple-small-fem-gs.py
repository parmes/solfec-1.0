# simple core model
import sys
sys.path.append ('inp/cores/inc')
from simple_core_base import *

# main module
#import rpdb2; rpdb2.start_embedded_debugger('a')

step = 1E-5
stop = 100 * step
scheme = 'EXP'
solver = 'GAUSS_SEIDEL_SOLVER'
plotconv = 0

solfec = SOLFEC ('DYNAMIC', step, 'out/cores/simple-small-fem-gs')
GRAVITY (solfec, (0, 0, -10))

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.7)
bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3)

if solver == 'GAUSS_SEIDEL_SOLVER':
  sv = GAUSS_SEIDEL_SOLVER (1E0, 50, 1E-6, failure = 'CONTINUE')
else:
  sv = NEWTON_SOLVER ('NONSMOOTH_HSW', 1E0, 20, 1E-4)
  sv.nonmonlength = 5
  sv.linmaxiter = 100

simple_core_create (0.0003, 0.0002, bulkmat, solfec, 'FINITE_ELEMENT', 'DEF_' + scheme, 'FINITE_ELEMENT', 'DEF_' + scheme, 2, 2, 2)

MERIT = []

def callback (sv):
  MERIT.append (sv.merhist)
  return 1

OUTPUT (solfec, 0.03)

if not VIEWER() and plotconv == 1: CALLBACK (solfec, step, sv, callback)

RUN (solfec, sv, stop)

if not VIEWER() and solfec.mode == 'WRITE' and plotconv == 1:
  try:
    import matplotlib.pyplot as plt

    for M in MERIT:
      plt.plot (list (range (0, len(M))), M)

    plt.semilogy (10)
    plt.title (solver + ': FEM MODEL - ' + scheme)
    plt.xlabel ('Iteration')
    plt.ylabel ('Merit function f(R)')
    plt.savefig ('out/cores/simple-small-fem-gs/simple-small-fem-' + solver + '_' + scheme + '.eps')
 
  except ImportError:
    pass # no reaction

if not VIEWER() and solfec.mode == 'READ':

  timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL', 'PARBAL', 'GSINIT', 'GSRUN', 'GSCOM', 'GSMCOM']
  dur = DURATION (solfec)
  th = HISTORY (solfec, timers, dur[0], dur[1])
  total = 0.0

  for i in range (0, len(timers)):
    sum = 0.0
    for tt in th [i+1]: sum += tt
    print timers [i], 'TIME:', sum
    if i < 6: total += sum

  print 'TOTAL TIME:', total
