# mathematical pendulum

PI = 3.14159265358979323846 

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

### main module ###

step = 0.0001

solfec = SOLFEC ('DYNAMIC', step, 'out/pendulum')

bulkmat = BULK_MATERIAL (solfec)

GRAVITY (solfec, (0, 0, -1), PI * PI)

#import rpdb2; rpdb2.start_embedded_debugger('a')

bod = pendulum_create (bulkmat, solfec)

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'EXIT')

previous = 0.0 # previous velocity component used by the termination callback
howmany = 0 # counter used by the termination callback

def termination (sol, bod):
  global previous, howmany

  if bod.velo [5] * previous <= 0:
    if howmany > 3:
      print '-----------------------------'
      print 'Pendulum period =', sol.time
      print '-----------------------------'
      return 0
    else: howmany = howmany + 1

  previous = bod.velo [5]
  return 1

if not VIEWER(): CALLBACK (solfec, step, (solfec, bod), termination)

RUN (solfec, gs, 10)
