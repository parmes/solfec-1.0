# simple core model
import sys
sys.path.append ('inp/cores/inc')
from simple_core_base import *

# main module

step = 0.0002
stop = 10
outfrq = 0.03
kinem = 'PSEUDO_RIGID'
solver = 'PENALTY'
scheme = 'DEFAULT'
shake = 'TRUE'
plotconv = 0
GEOMETRIC_EPSILON (1E-6)
WARNINGS ('OFF')

if kinem == 'PSEUDO_RIGID' or kinem == 'FINITE_ELEMENT': scheme = 'DEF_IMP'

if kinem == 'PSEUDO_RIGID': kinstr = '-prb-'
elif kinem == 'RIGID': kinstr = '-rig-'
elif kinem == 'FINITE_ELEMENT': kinstr = '-fem-'
else:
  print 'Uknown kinematics'
  sys.exit (1)

if solver == 'GAUSS_SEIDEL': solstr = 'gs'
elif solver == 'NEWTON': solstr = 'nt'
elif solver == 'PENALTY': solstr = 'pn'
elif solver == 'TEST': solstr = 'ts'
else:
  print 'Uknown solver'
  sys.exit (1)

outdir =  'out/cores/simple-small' + kinstr + solstr
solfec = SOLFEC ('DYNAMIC', step, outdir)

if shake == 'TRUE':
  GRAVITY (solfec, (0, 0, -10))
else:
  acc = TIME_SERIES ('inp/cores/inc/acc-0.dat')
  GRAVITY (solfec, (acc, 0, -10))

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.7, spring = 1E6, dashpot = 1E3)
bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3)

if solver == 'GAUSS_SEIDEL':
  sv = GAUSS_SEIDEL_SOLVER (1E-3, 260, 1E-6, diagsolver = 'SEMISMOOTH_NEWTON')
  sv.reverse = 'ON'
elif solver == 'NEWTON':
  sv = NEWTON_SOLVER (1E-6, 250, delta = 1E-5)
elif solver == 'PENALTY':
  sv = PENALTY_SOLVER ('IMPLICIT')
elif solver == 'TEST':
  sv = TEST_SOLVER (1E-6, 250)

MERIT = []
def callback (sv):
  MERIT.append (sv.merhist)
  return 1

simple_core_create (0.0003, 0.0002, bulkmat, solfec, kinem, scheme, shake, 4, 4, 4)

UNPHYSICAL_PENETRATION (solfec, 0.02)
IMBALANCE_TOLERANCE (solfec, 1.1)
OUTPUT (solfec, outfrq)
if not VIEWER() and plotconv == 1: CALLBACK (solfec, step, sv, callback)
RUN (solfec, sv, stop)

if not VIEWER() and solfec.mode == 'WRITE' and plotconv == 1:
  try:
    import matplotlib.pyplot as plt

    for M in MERIT:
      plt.plot (list (range (0, len(M))), M)

    plt.semilogy (10)
    plt.title (solver + ': ' + str (int(stop/step)) + ' steps, ' + kinem + ' model')
    plt.xlabel ('Iteration')
    plt.ylabel ('Merit function f(R)')
    plt.savefig (outdir + '/iter-history' + kinstr + solstr + '.eps')
 
  except ImportError:
    pass # no reaction

if not VIEWER() and solfec.mode == 'READ':

  data = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL', 'PARBAL',
	  'MERIT']

  dur = DURATION (solfec)
  th = HISTORY (solfec, data, dur[0], dur[1], progress = 'ON')

  total = 0.0
  for i in range (0, 6):
    sum = 0.0
    for tt in th [i+1]: sum += tt
    print data [i], 'TIME:', sum
    if i < 6: total += sum

  print 'TOTAL TIME:', total

  try:
    import matplotlib.pyplot as plt

    plt.plot (th[0], th[7])
    plt.semilogy (10)
    plt.title (solver + ': ' + str (int(stop/step)) + ' steps, ' + kinem + ' model')
    plt.xlabel ('Time')
    plt.ylabel ('Merit function f(R)')
    plt.savefig (outdir + '/merit-history' + kinstr + solstr + '.eps')
 
  except ImportError:
    pass # no reaction
