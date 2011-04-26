# sheared stack of cubes

from math import pow

nsides = 10
lside = 0.1
fric = 0.25
kinem = 'FINITE_ELEMENT'
solv = 'ns'
sparea = 0.01 * lside * lside
gravacc = 10
step = 1E-3
nqst = 50
stop = 200 * step

def cube (x, y, z, a, b, c, div, sur, vol):

  nodes = [0, 0, 0,
	   a, 0, 0,
	   a, b, 0,
	   0, b, 0,
	   0, 0, c,
	   a, 0, c,
	   a, b, c,
	   0, b, c]

  shp = HEX (nodes, div, div, div, vol, [sur, sur, sur, sur, sur, sur])

  TRANSLATE (shp, (x, y, z))

  return shp

### main module ###
#import rpdb2; rpdb2.start_embedded_debugger('a')

outdir = 'out/shearcubes' + '_' + kinem + '_' + solv + '_f' + str(fric)

solfec = SOLFEC ('DYNAMIC', step, outdir)

CONTACT_SPARSIFY (solfec, 0.05, minarea = sparea)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = fric)

hardmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E10, poisson = 0.25, density = 2E3)

GRAVITY (solfec, (0, 0, -gravacc))

# create base
shp = cube (0, 0, -lside, nsides*lside, nsides*lside, lside, 1, 1, 1)
BODY (solfec, 'OBSTACLE', shp, hardmat)

# create the pushing wall
shp = cube (0, -lside, 0, nsides*lside, lside, nsides*lside, 1, 2, 2)
bod = BODY (solfec, 'RIGID', shp, hardmat)
FIX_POINT (bod, (nsides*lside, 0, 0))
FIX_POINT (bod, (0, 0, 0))
tms = TIME_SERIES ([0, 0, nqst * step, 0, (nqst+1) * step, nsides * lside, stop, nsides * lside])
SET_VELOCITY (bod, bod.center, (0, 1, 0), tms)

# create the remaining bricks
for i in range (nsides):
  for j in range (nsides):
    for k in range (nsides):
      x = i*lside
      y = j*lside
      z = k*lside
      shp = cube (x, y, z, lside, lside, lside, 2, 3, 3)
      bod = BODY (solfec, kinem, shp, hardmat)
      if bod.scheme != 'RIGID': bod.scheme = 'DEF_LIM'

# damping
for bod in solfec.bodies: bod.damping = 0.001

if solv == 'gs':
  sv = GAUSS_SEIDEL_SOLVER (1E-2, 1000, 1E-7)
else:
  sv = NEWTON_SOLVER (1E-7, 1000, theta = 0.25)

OUTPUT (solfec, step, 'ON')

RUN (solfec, sv, stop)

if not VIEWER() and solfec.mode == 'READ':

  timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL', 'PARBAL']
  dur = DURATION (solfec)
  th = HISTORY (solfec, timers, dur[0], dur[1])
  total = 0.0

  for i in range (0, len(timers)):
    sum = 0.0
    for tt in th [i+1]: sum += tt
    print timers [i], 'TIME:', sum
    if i < 6: total += sum

  print 'TOTAL TIME:', total

  try:
    import matplotlib.pyplot as plt
    th = HISTORY (solfec, [(solfec, 'KINETIC'), (solfec, 'INTERNAL'), (solfec, 'EXTERNAL')], 0, stop)
    plt.plot (th [0], th [1], label='kin')
    plt.plot (th [0], th [2], label='int')
    plt.plot (th [0], th [3], label='ext')
    plt.axis (xmin = 0, xmax = stop)
    plt.xlabel ('Time')
    plt.ylabel ('Energy')
    plt.legend(loc = 'upper right')
    plt.savefig (outdir + '/ene.eps')
  except ImportError:
    pass # no reaction
