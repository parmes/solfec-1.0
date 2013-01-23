# 492 spheres falling into a container
import time

# global data

a = 0.13
b = 0.13
c = 0.010
r = 0.008

stop = 1.0

# simulation routine

def sphdrop (step, solver):

  solfec = SOLFEC ('DYNAMIC', step, 'out/sphdrop_' + str(step) + '_stop_' + str(solver) + '_solver')

  bulk = BULK_MATERIAL (solfec, model='KIRCHHOFF', young=1E9, poisson=0.25, density=1E3)

  surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3,
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
      for i in range (0,1):
	if s != t:
	  sph = SPHERE ((a*s/5 - i % 2 * 0.01, b*t/5 + i % 2 * 0.01, a + i*0.03), r, 1, 1)
	  if solver == 'penalty': BODY (solfec, 'RIGID', sph, bulk)
	  else:
	    bod = BODY (solfec, 'PSEUDO_RIGID', sph, bulk)
	    bod.scheme = 'DEF_LIM'
	    bod.damping = 0.0

  GRAVITY (solfec, (0, 0, -10))
  if solver == 'penalty':
    slv = PENALTY_SOLVER ('EXPLICIT')
  else:
    slv = NEWTON_SOLVER (1E-6, delta = 1E-6)

  OUTPUT (solfec, 0.01)

  t0 = time.time()

  RUN (solfec, slv, stop)

  t1 = time.time() - t0

  return t1

# main module

t1 = sphdrop (2E-5, 'penalty')
t2 = sphdrop (2E-4, 'newton')

print 'Penalty, step 2E-5, runtime = ', t1
print 'Newton, step 1E-4, runtime = ', t2
