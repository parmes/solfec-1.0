# Lourenco wall example
import math

step = 0.0001
stop = step * 500
NWIDTH = 5
NHEIGHT = 9
formul = 'TL'

solfec = SOLFEC ('QUASI_STATIC', step, 'out/lourenco')

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
    BODY (solfec, 'FINITE_ELEMENT', brick, bulk, form = formul)
    point = ((NWIDTH-1)*a*2+a/2, 0, j*c)
    brick = HEX (nodes, 2, 2, 2, 1, [1, 1, 1, 1, 1, 1])
    SCALE (brick, (0.5, 1, 1))
    TRANSLATE (brick, point)
    BODY (solfec, 'FINITE_ELEMENT', brick, bulk, form = formul)
    for i in range (0, NWIDTH-1):
      point = (a+i*a*2, 0, j*c)
      brick = HEX (nodes, 2, 2, 2, 1, [1, 1, 1, 1, 1, 1])
      TRANSLATE (brick, point)
      BODY (solfec, 'FINITE_ELEMENT', brick, bulk, form = formul)
  else:
    for i in range (NWIDTH):
      point = (i*a*2, 0, j*c)
      brick = HEX (nodes, 2, 2, 2, 1, [1, 1, 1, 1, 1, 1])
      TRANSLATE (brick, point)
      BODY (solfec, 'FINITE_ELEMENT', brick, bulk, form = formul)

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
#SET_VELOCITY (top, (-a, 0, c * NHEIGHT + c/2), (1, 0, 0), 0.015)
tms = TIME_SERIES ([0, 0, 50 * step, 0, stop, 20E3])
FORCE (top, 'SPATIAL', (-a, 0, c * NHEIGHT + c/2), (1, 0, 0), tms)
FORCE (top, 'SPATIAL', (2 * a * NWIDTH / 2 - a, 0, c * NHEIGHT + c/2), (0, 0, -1), 30E3)

#gs = GAUSS_SEIDEL_SOLVER (1E-5, 1000)
gs = NEWTON_SOLVER (1E-9, 1000)

OUTPUT (solfec, 1E-3, compression = 'ON')

RUN (solfec, gs, stop)

if not VIEWER() and solfec.mode == 'READ':

  timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL']
  dur = DURATION (solfec)
  th = HISTORY (solfec, timers, dur[0], dur[1])
  total = 0.0

  for i in range (0, len (timers)):
    sum = 0.0
    for tt in th [i+1]: sum += tt
    print timers [i], 'TIME:', sum
    total += sum

  print 'TOTAL TIME:', total
