# 492 boxes falling into a container
import time

# global data

a = 0.13
b = 0.13
c = 0.010
r = 0.008

stop = 1.0

# simulation routine

def boxdrop (step, solver):

  solfec = SOLFEC ('DYNAMIC', step, 'out/sphdrop_' + str(step) + '_stop_' + str(solver) + '_solver')

  bulk = BULK_MATERIAL (solfec, model='KIRCHHOFF', young=1E9, poisson=0.25, density=1E3)

  surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.1,
	    spring = (4./3.)*(1E9/(1-0.25**2.0))*(r**0.5), dashpot = -1, hpow = 1.5)

  base = HEX ([ 0 , 0 , 0 ,
		a , 0 , 0 ,
		a , b , 0 ,
		0 , b , 0 ,
		0 , 0 , -c ,
		a , 0 , -c ,
		a , b , -c ,
		0 , b , -c ] , 1, 1, 1, 1, [1]*6)

  base  = COPY (base)
  wall1 = ROTATE (COPY (base), (a , 0 , -c), (a , 0 , 0), 90)
  wall2 = ROTATE (COPY (base), (a , 0 , -c), (0 , b , 0), 90)
  wall3 = ROTATE (COPY (base), (a , b , -c), (a , 0 , 0), -90)
  wall4 = ROTATE (COPY (base), (0 , b , -c), (0 , b , 0), -90)

  container = [base, wall1, wall2, wall3, wall4]

  BODY (solfec, 'OBSTACLE',container, bulk)

  for s in range (1,5):
    for t in range (1,5):
      for i in range (0,25):
	if s != t:
          p = (a*s/5 - i % 2 * 0.01, b*t/5 + i % 2 * 0.01, a + i*0.03)
          shp = HEX ([p[0]-r, p[1]-r, p[2]+r,
                      p[0]+r, p[1]-r, p[2]+r,
                      p[0]+r, p[1]+r, p[2]+r,
                      p[0]-r, p[1]+r, p[2]+r,
                      p[0]-r, p[1]-r, p[2]-r,
                      p[0]+r, p[1]-r, p[2]-r,
                      p[0]+r, p[1]+r, p[2]-r,
                      p[0]-r, p[1]+r, p[2]-r], 1, 1, 1, 1, [1]*6)
	  if solver == 'penalty': BODY (solfec, 'RIGID', shp, bulk)
	  else:
	    bod = BODY (solfec, 'FINITE_ELEMENT', shp, bulk, form = 'BC')
	    bod.scheme = 'DEF_LIM'
	    bod.damping = 1E-5

  GRAVITY (solfec, (0, 0, -10))
  if solver == 'penalty':
    slv = PENALTY_SOLVER ('EXPLICIT')
  else:
    slv = NEWTON_SOLVER (1E-6, delta = 1E-6)

  OUTPUT (solfec, 0.01)

  t0 = time.time()

  RUN (solfec, slv, stop)

  t1 = time.time() - t0

  return (solfec, t1)

# main module

s1, t1 = boxdrop (3E-5, 'penalty')
s2, t2 = boxdrop (3E-4, 'newton')

if not VIEWER() and s1.mode == 'WRITE' and s2.mode == 'WRITE':

  print 'DEM, step 3E-5, runtime = ', t1
  print 'CD, step 3E-4, runtime = ', t2
  print 'Run again without arguments for more refined timings!'

if not VIEWER() and s1.mode == 'READ' and s2.mode == 'READ':

  timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL']
  th1 = HISTORY (s1, timers, 0, 1)
  th2 = HISTORY (s2, timers, 0, 1)

  tt1 = [0, 0, 0]
  tt2 = [0, 0, 0]
  for i in range (0, 5):
    su = 0.0
    for tt in th1 [i+1]: su += tt
    if i == 0: tt1 [0] = su
    if i in (1, 2): tt1 [1] += su
    if i in (3, 4): tt1 [2] += su
    su = 0.0
    for tt in th2 [i+1]: su += tt
    if i == 0: tt2 [0] = su
    if i in (1, 2): tt2 [1] += su
    if i in (3, 4): tt2 [2] += su

  tot1 = tt1[0]+tt1[1]+tt1[2]
  tot2 = tt2[0]+tt2[1]+tt2[2]

  print 'DEM total time:', tot1
  print 'DEM time integration:', 100 * tt1[0] / tot1, '%'
  print 'DEM contact detection:', 100 * tt1[1] / tot1, '%'
  print 'DEM contact solution:', 100 * tt1[2] / tot1, '%'

  print 'CD (3E-4) total time:', tot2
  print 'CD (3E-4) time integration:', 100 * tt2[0] / tot2, '%'
  print 'CD (3E-4) contact detection:', 100 * tt2[1] / tot2, '%'
  print 'CD (3E-4) contact solution:', 100 * tt2[2] / tot2, '%'
