from random import randint
from random import random
from random import seed
from math import sin

def make_wee_box (x, y, z, r, kinem, scheme, material, solfec):
  nodes = [x - r, y - r, z - r,
           x + r, y - r, z - r,
	   x + r, y + r, z - r,
	   x - r, y + r, z - r,
           x - r, y - r, z + r,
           x + r, y - r, z + r,
	   x + r, y + r, z + r,
	   x - r, y + r, z + r]

  msh = HEX (nodes, 1, 1, 1, 0, [0, 0, 0, 0, 0, 0])
  bod = BODY (solfec, kinem, msh, material)
  bod.scheme = scheme

def make_big_box (x, y, z, wx, wy, wz, material, solfec):
  box = HULL ([x, y, z,
               x+wx, y, z,
	       x+wx, y+wy, z,
	       x, y+wy, z,
               x, y, z+wz,
               x+wx, y, z+wz,
	       x+wx, y+wy, z+wz,
	       x, y+wy, z+wz], 2, 2)

  wallx1 = SCALE (COPY (box), (1, 0.1, 1))
  wallx2 = TRANSLATE (COPY (wallx1), (0, 0.9*wy, 0))
  wally1 = SCALE (COPY (box), (0.1, 0.8, 1))
  wally2 = TRANSLATE (COPY (wally1), (0.9*wz, 0, 0))
  TRANSLATE ([wally1, wally2], (0, 0.1*wy, 0))
  wallz1 = SCALE (COPY (box), (1, 1, 0.1))
  TRANSLATE (wallz1, (0, 0, -0.1*wz))

  shape = [wallx1, wallx2, wally1, wally2, wallz1]
  BODY (solfec, 'OBSTACLE', shape, material)

#kinem = 'FINITE_ELEMENT'
kinem = 'PSEUDO_RIGID'
scheme = 'DEF_IMP'
#kinem = 'RIGID'
#scheme = 'RIG_NEG'
step = 0.003
stop = 10.0
n = 10

solfec = SOLFEC ('DYNAMIC', step, 'out/boxsoup')
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = 1.0)
bulk = BULK_MATERIAL (solfec, 'KIRCHHOFF', young = 15E8, poisson = 0.25, density = 1.8E3)
GRAVITY (solfec, (0, 0, -10))
gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000)

IMBALANCE_TOLERANCE (solfec, 2.0, lockdir = 'ON')

make_big_box (0, 0, 0, n, n, n, bulk, solfec)

for i in range (2,n-1, 2):
  for j in range (2,n-1, 2):
    for k in range (2, n-1, 2):
      make_wee_box (i, j, k, 0.5, kinem, scheme, bulk, solfec)

#OUTPUT (solfec, 0.01)
RUN (solfec, gs, stop)

if not VIEWER() and solfec.mode == 'READ':
  try:
    import matplotlib.pyplot as plt
    th = HISTORY (solfec, [(solfec, 'KINETIC'), (solfec, 'INTERNAL'), (solfec, 'EXTERNAL'), (solfec, 'CONTACT') ], 0, stop, progress = 'ON')
    plt.plot (th [0], th [1], label='KIN')
    plt.plot (th [0], th [2], label='INT')
    plt.plot (th [0], th [3], label='EXT')
    tot = []
    for i in range(0, len (th[0])): tot.append (th[1][i] + th[2][i] - th[3][i])
    plt.plot (th [0], tot, label='TOT')
    plt.axis (xmin = 0, xmax = stop)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Energy [J]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/boxsoup/boxsoup.eps')
  except ImportError:
    pass # no reaction
