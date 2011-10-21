# a drum mixing example

from random import randint
from random import random
from random import seed
from math import sin
from math import cos

def stone (x, y, z, r, kind, material, solfec):
  m = randint (8, 32)
  points = []
  for n in range (m):
    points.append (x + r * (1.0 - random()) * 2.0)
    points.append (y + r * (1.0 - random()) * 2.0)
    points.append (z + r * (1.0 - random()) * 2.0)

  hul = HULL (points, 1, 1)
  if kind == 'rig':
    BODY (solfec, 'RIGID', hul, material)
  elif kind == 'fem':
    m = ROUGH_HEX (hul, 2, 2, 2)
    b = BODY (solfec, 'FINITE_ELEMENT', hul, material, form = 'BC', mesh = m)
    b.scheme = 'DEF_LIM'
    b.damping = 1E-3

def wheel (x, y, z, r, t, v, s):
  a = 0.0
  points = []
  while a < 6.2:
    px = x + r * sin (a)
    py = y
    pz = z + r * cos (a)
    points.append (px)
    points.append (py)
    points.append (pz)
    points.append (px)
    points.append (py+t)
    points.append (pz)
    a += 0.2
  return HULL (points, v, s)

 
seed (1)
kin = 'rig'
step = 0.001
skip = 0.01
dura = 20.0

solfec = SOLFEC ('DYNAMIC', step, 'out/drum' + kin)
surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.4)
drumat = BULK_MATERIAL (solfec)
bodmat = BULK_MATERIAL (solfec, density = drumat.density / 100)
GRAVITY (solfec, (0, 0, -9.8))
sv = GAUSS_SEIDEL_SOLVER (1E-3, 20)
OUTPUT (solfec, skip)

# drum
pip = PIPE ((0, 0, 0), (0, 0.5, 0), 1, 0.05, 1, 32, 1, 1, [1, 2, 3, 4])
el1 = ELLIP ((-0.9, 0.25, 0), (0.1, 0.24, 0.1), 1, 1)
el2 = ELLIP ((0.9, 0.25, 0), (0.1, 0.24, 0.1), 1, 1)
sp1 = SPHERE ((0, -1, 0), 0.05, 1, 1)
sp2 = SPHERE ((0, 1.5, 0), 0.05, 1, 1)
shape = [pip, el1, el2, sp1, sp2]
bod = BODY (solfec, 'RIGID', shape, drumat)
FIX_POINT (bod, (0, -1, 0))
FIX_POINT (bod, (0, 1.5, 0))

# sides
el1 = wheel (0, -0.002, 0, 1.05, -0.1, 1, 1)
BODY (solfec, 'OBSTACLE', el1, drumat)
el2 = wheel (0, 0.502, 0, 1.05, 0.1, 1, 1)
BODY (solfec, 'OBSTACLE', el2, drumat)

# stones
for i in range (0, 20):
  y = 0.05
  while y < 0.4:
    x = -0.7
    while x < 0.7:
      stone (x, y, 0.1, 0.05, kin, bodmat, solfec)
      x += 0.1
    y += 0.1
  RUN (solfec, sv, 0.15)

# bar
len = 50
bar = HULL ([-len, 0, 1.05, -len, 0.5, 1.05, -len, 0.5, 2, -len, 0, 2,
             0, 0, 1.05, 0, 0.5, 1.05, 0, 0.5, 2, 0, 0, 2], 1, 4)
bod = BODY (solfec, 'RIGID', bar, drumat)
FIX_DIRECTION (bod, (-len, 0, 1.05), (0, 0, 1))
FIX_DIRECTION (bod, (-len, 0.5, 1.05), (0, 0, 1))
FIX_DIRECTION (bod, (-len, 0.5, 1.05), (0, 1, 0))
SET_VELOCITY (bod, (-len, 0.25, 1.5), (1, 0, 0), 1)

# solver and run
RUN (solfec, sv, dura)
