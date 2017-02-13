PI = 3.14159265358979323846 
previous = 0.0 # previous velocity component used by the termination callback
howmany = 0 # counter used by the termination callback
gravity = PI * PI
period = 2.36
numper = 10
step = 0.001
stop = numper*period+0.5*period

def pendulum_create (material, solfec):
  c = (1, 0, 1)
  J = (1, 0, 0, 0, 1, 0, 0, 0, 1)
  p = (1, 0, 1)
  q = (0, 0, 1)
  sph = SPHERE (p, 0.05, 0, 0);
  bod = BODY (solfec, 'RIGID', sph, material)
  BODY_CHARS (bod, 1.0, 1.0, c, J)
  PUT_RIGID_LINK (bod, None, p, q)
  return bod

solfec = SOLFEC ('DYNAMIC', step, 'out/pendulum')
bulkmat = BULK_MATERIAL (solfec)
GRAVITY (solfec, (0, 0, -gravity))
bod = pendulum_create (bulkmat, solfec)
gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'EXIT')
RUN (solfec, gs, stop)
if solfec.mode == 'WRITE': print 'Run once more now to produce figures...'

if not VIEWER() and solfec.mode == 'READ':
  print 'Producing figures...'
  his = HISTORY (solfec,[(bod, 'KINETIC'), (bod, bod.center, 'CZ'), (bod, bod.center, 'VZ')], 0, stop)
  T   = his[0]
  KIN = his[1]
  POT = []
  TOT = []
  signchange = 0
  for i in range(0,len(T)):
    kin = KIN[i]
    pot = bod.mass * gravity * max (his[2][i], 0)
    POT.append(pot)
    TOT.append(kin+pot)
    if i < len(T)-1 and his[3][i]*his[3][i+1] <= 0:
      if signchange == 4:
	exact = period
	error = abs (T[i] - exact) / exact
	print '1 swing -->',
	print 'Calculated pendulum period: %.3f;' % his[0][i], 'Feference value: %.3f;' % exact, 'Relative error: %.3f' % error
	print '           Total energy:', TOT[-1]
      if signchange == 4*numper:
	exact = period*numper
	error = abs (T[i] - exact) / exact
	print '%d swings -->' % numper,
	print 'Calculated pendulum period: %.3f;' % his[0][i], 'Feference value: %.3f;' % exact, 'Relative error: %.3f' % error
	print '           Total energy:', TOT[-1]
      signchange = signchange + 1
  try:
    import matplotlib.pyplot as plt
    plt.clf ()
    plt.plot (T, KIN, label='Kinetic')
    plt.plot (T, POT, label='Potential')
    plt.plot (T, TOT, label='Total')
    plt.axis (xmin = 0, xmax = period, ymin = 0, ymax = 10)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Energy [J]')
    plt.legend (loc = 'upper right')
    plt.savefig ('validation/pendulum/pendulum.png')
  except (ImportError, RuntimeError):
    import sys
    print "Unexpected error:", sys.exc_info()[1]
    print "Plotting has failed!"
    pass
