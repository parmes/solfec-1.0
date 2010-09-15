# Newton cradle
from math import sin
from math import cos

theta = 3.14159 / 8.0
step = 0.001
stop = 150
ofrq = 0.01

solfec = SOLFEC ('DYNAMIC', step, 'out/newton-cradle')

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E7, poisson = 0.3, density = 1E3)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = 1.0)

GRAVITY (solfec, (0, 0, -9.81))

gs = GAUSS_SEIDEL_SOLVER (1E-10, 100)

GEOMETRIC_EPSILON (1E-15)

for i in range (0, 5):
  x = i * (0.1 + 1E-10)
  c = (x, 0, 0)
  p = (x, 0, 0.5)

  if x == 0:
    c = TRANSLATE (c, (-0.5 * sin (theta), 0, 0.5 * (1 - cos (theta))))
   
  sph = SPHERE (c, 0.05, 0, 0)
  bod = BODY (solfec, 'RIGID', sph, bulkmat)
  PUT_RIGID_LINK (bod, None, c, p)


OUTPUT (solfec, ofrq)
RUN (solfec, gs, stop)

if not VIEWER() and solfec.mode == 'READ':
  try:
    import matplotlib.pyplot as plt
    th = HISTORY (solfec, [(solfec, 'KINETIC'), (solfec, 'EXTERNAL')], 0, stop)
    plt.plot (th [0], th [1], label='KIN')
    plt.plot (th [0], th [2], label='EXT')
    tot = []
    for i in range(0, len (th[0])): tot.append (th[1][i] - th[2][i])
    plt.plot (th [0], tot, label='TOT')
    plt.axis (xmin = 0, xmax = stop)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Energy [J]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/newton-cradle/newton-cradle-ene.eps')
  except ImportError:
    pass # no reaction
