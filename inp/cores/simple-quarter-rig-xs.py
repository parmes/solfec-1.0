# simple core model
import sys
sys.path.append ('inp/cores/inc')
from simple_core_base import *

# main module
#import rpdb2; rpdb2.start_embedded_debugger('a')

step = 1E-3
stop = 11

solfec = SOLFEC ('DYNAMIC', step, 'out/cores/simple-quarter-rig-xs')
GRAVITY (solfec, (0, 0, -1), 10)

= SURFACE_MATERIAL (solfec, model = 'SPRING_DASHPOT', spring = 1E6, dashpot=1E3, friction = 0.7)
bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3)

xs = EXPLICIT_SOLVER ()

simple_core_create (0.0003, 0.0002, bulkmat, solfec, 'RIGID', 'DEFAULT', 'RIGID', 'DEFAULT', 10, 10, 12)

OUTPUT (solfec, 0.03)
RUN (solfec, xs, stop)

if not VIEWER() and solfec.mode == 'READ':

  timers = ['TIMINT', 'CONDET', 'LOCDYN', 'CONSOL', 'PARBAL']
  dur = DURATION (solfec)
  th = HISTORY (solfec, timers, dur[0], dur[1])
  total = 0.0

  for i in range (0, 5):
    sum = 0.0
    for tt in th [i+1]: sum += tt
    print timers [i], 'TIME:', sum
    total += sum

  print 'TOTAL TIME:', total
