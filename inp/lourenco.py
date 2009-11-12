# Lourenco wall example
import math

step = 1E-5
NWIDTH = 8
NHEIGHT = 8
damping = 1.0

solfec = SOLFEC ('DYNAMIC', step, 'out/lourenco')

sur = SURFACE_MATERIAL (solfec,
                        model = 'SIGNORINI_COULOMB',
                        friction = 0.3)

bulk = BULK_MATERIAL (solfec,
                      model = 'KIRCHHOFF',
		      young = 1E9,
		      poisson = 0.25,
		      density = 1E3)


gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'CONTINUE', diagsolver = 'PROJECTED_GRADIENT')

GRAVITY (solfec, (0, 0, -1), 10)

a = 0.10
b = 0.05
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
  if j % 2:
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
FORCE (top, 'SPATIAL', (a * (NWIDTH-1)/2, 0, c * NHEIGHT + c/2), (1, 0, 0), 1000)

OUTPUT (solfec, 10 * step)

RUN (solfec, gs, 10000 * step, 50 * step)

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
