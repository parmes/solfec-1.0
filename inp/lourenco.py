# Lourenco wall example
import math

step = 8E-6
NWIDTH = 5
NHEIGHT = 9
damping = 1.0

solfec = SOLFEC ('DYNAMIC', step, 'out/lourenco')

sur = SURFACE_MATERIAL (solfec,
                        model = 'SIGNORINI_COULOMB',
                        friction = 0.62)

sur = SURFACE_MATERIAL (solfec,
                        model = 'SIGNORINI_COULOMB',
			surf1 = 1,
			surf2 = 3,
                        friction = 0.62,
			cohesion = 0.3E6)

bulk = BULK_MATERIAL (solfec,
                      model = 'KIRCHHOFF',
		      young = 15.5E9,
		      poisson = 0.2,
		      density = 1E3)


GRAVITY (solfec, (0, 0, -1), 10)

a = 0.10
b = 0.10
c = 0.10

nodes = [-a, -b, 0,
          a, -b, 0,
          a,  b, 0,
         -a,  b, 0,
         -a, -b, c,
          a, -b, c,
          a,  b, c,
         -a,  b, c]

for j in range (NHEIGHT):
  if j % 2 == 0:
    point = (-a/2, 0, j*c)
    brick = HEX (nodes, 2, 2, 2, 1, [1, 1, 1, 1, 1, 1])
    SCALE (brick, (0.5, 1, 1))
    TRANSLATE (brick, point)
    b = BODY (solfec, 'FINITE_ELEMENT', brick, bulk)
    b.damping = damping
    point = ((NWIDTH-1)*a*2+a/2, 0, j*c)
    brick = HEX (nodes, 2, 2, 2, 1, [1, 1, 1, 1, 1, 1])
    SCALE (brick, (0.5, 1, 1))
    TRANSLATE (brick, point)
    b = BODY (solfec, 'FINITE_ELEMENT', brick, bulk)
    b.damping = damping
    for i in range (0, NWIDTH-1):
      point = (a+i*a*2, 0, j*c)
      brick = HEX (nodes, 2, 2, 2, 1, [1, 1, 1, 1, 1, 1])
      TRANSLATE (brick, point)
      b = BODY (solfec, 'FINITE_ELEMENT', brick, bulk)
      b.damping = damping
  else:
    for i in range (NWIDTH):
      point = (i*a*2, 0, j*c)
      brick = HEX (nodes, 2, 2, 2, 1, [1, 1, 1, 1, 1, 1])
      TRANSLATE (brick, point)
      b = BODY (solfec, 'FINITE_ELEMENT', brick, bulk)
      b.damping = damping

base = HEX (nodes, 1, 1, 1, 2, [2, 2, 2, 2, 2, 2])
TRANSLATE (base, (a, 0, -c))
SCALE (base, (NWIDTH + 1, 1, 1))
TRANSLATE (base, (-2*a, 0, 0))
BODY (solfec, 'OBSTACLE', base, bulk)

base = HEX (nodes, 1, 1, 1, 3, [3, 3, 3, 3, 3, 3])
TRANSLATE (base, (a, 0, -c))
SCALE (base, (NWIDTH + 1, 1, 1))
TRANSLATE (base, (-2*a, 0, c*(NHEIGHT+1)))
top = BODY (solfec, 'RIGID', base, bulk)
#SET_VELOCITY (solfec, top, (-a, 0, c * NHEIGHT + c/2), (1, 0, 0), 0.015)
tms = TIME_SERIES ([0, 0, 0.5, 0, 1.0, 20E3])
FORCE (top, 'SPATIAL', (-a, 0, c * NHEIGHT + c/2), (1, 0, 0), tms)
FORCE (top, 'SPATIAL', (2 * a * NWIDTH / 2 - a, 0, c * NHEIGHT + c/2), (0, 0, -1), 30E3)

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'CONTINUE', diagsolver = 'PROJECTED_GRADIENT')

OUTPUT (solfec, 1E-3, compression = 'FASTLZ')

RUN (solfec, gs, 1.0)

if not VIEWER() and solfec.mode == 'READ':

  timers = ['TIMINT', 'CONDET', 'LOCDYN', 'CONSOL']
  dur = DURATION (solfec)
  total = 0.0

  for timer in timers:
    th = TIMING_HISTORY (solfec, timer, dur[0], dur[1])
    sum = 0.0
    for tt in th [1]: sum += tt
    total += sum
    print timer, 'TIME:', sum

  print 'TOTAL TIME:', total
