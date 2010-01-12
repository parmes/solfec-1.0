# mathematical pendulum

PI = 3.14159265358979323846 
previous = 0.0 # previous velocity component used by the termination callback
howmany = 0 # counter used by the termination callback

def pendulum_create (material, solfec):
  c = (1, 0, 1)
  J = (1, 0, 0, 0, 1, 0, 0, 0, 1)
  p = (1, 0, 1)
  q = (0, 0, 1)

  sph = SPHERE (p, 0.05, 0, 0);
  bod = BODY (solfec, 'RIGID', sph, material)
  BODY_CHARS (bod, 1.0, 1.0, c, J)
  PUT_RIGID_LINK (solfec, bod, None, p, q)

  return bod

# run is terminated after one swing
def termination (sol, bod):
  global previous, howmany

  if bod.velo [5] * previous <= 0:
    if howmany > 3:
      exact = 2.36
      error = abs (sol.time - exact) / exact
      if (error < 0.01): print 'PASSED'
      else:
	print 'FAILED'
        print '(', 'Computed pendulum period was %.3f' % sol.time, 'while the reference value is %.3f' % exact, ')'
      return 0
    else: howmany = howmany + 1

  previous = bod.velo [5]
  return 1

# main module
#import rpdb2; rpdb2.start_embedded_debugger('a')

step = 0.001
stop = 10

solfec = SOLFEC ('DYNAMIC', step, 'out/tests/math-pendulum')
solfec.verbose = 'OFF'

bulkmat = BULK_MATERIAL (solfec)

GRAVITY (solfec, (0, 0, -1), PI * PI)

bod = pendulum_create (bulkmat, solfec)

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'EXIT')

if not VIEWER(): CALLBACK (solfec, step, (solfec, bod), termination)

if not VIEWER () and solfec.mode == 'READ': print '\nPrevious test results exist. Please "make del" and rerun tests'
else: RUN (solfec, gs, stop)
