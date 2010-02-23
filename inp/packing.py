from random import randint
from random import random
from random import seed
from math import sin

def make_particle (x, y, z, r, kind, material, solfec):
  m = randint (8, 32)
  points = []
  r = 0.2 * r + 0.8 * r * random ()
  for n in range (m):
    points.append (x + r * (1.0 - random()) * 2.0)
    points.append (y + r * (1.0 - random()) * 2.0)
    points.append (z + r * (1.0 - random()) * 2.0)

  hull = HULL (points, 1, 1)
  BODY (solfec, kind, hull, material)


def make_time_history (step, time, per, amp):

  t = 0.0
  points = []
  for i in range (0, int (time/step)):
    points.append (t)
    points.append (amp * (1.0 + 9.0 * (time-t)/time) * sin (t/per))
    t = t + step

  return TIME_SERIES (points)


def make_box (x, y, z, wx, wy, wz, material, solfec, velo):
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
  bod = BODY (solfec, 'RIGID', shape, material)
  FIX_DIRECTION (solfec, bod, (x, y, z), (0, 0, 1))
  FIX_DIRECTION (solfec, bod, (x+wx, y, z), (0, 0, 1))
  FIX_DIRECTION (solfec, bod, (x+wx, y+wy, z), (0, 0, 1))
  FIX_DIRECTION (solfec, bod, (x, y+wy, z), (0, 0, 1))
  FIX_DIRECTION (solfec, bod, (x, y, z), (0, 1, 0))
  SET_VELOCITY (solfec, bod, (x+wx*0.5, y+wy*0.5, z), (1, 0, 0), velo)

step = 0.001
stop = 10.0
velo = make_time_history (step, stop+1, 0.01, 1)
n = 10

seed (1)
solfec = SOLFEC ('DYNAMIC', step, 'out/packing')
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5)
bulk = BULK_MATERIAL (solfec, 'KIRCHHOFF', young = 1E5, poisson = 0.25, density = 1E3)
GRAVITY (solfec, (0, 0, -1), 10)
gs = GAUSS_SEIDEL_SOLVER (1E-2, 10)

make_box (0, 0, 0, n, n, n, bulk, solfec, velo)

EXTENTS (solfec, (-1, -1, -1, n+1, n+1, n+1))

for i in range (1,n-1):
  for j in range (1,n-1):
    for k in range (1, 5*n):
      make_particle (i, j, k, 0.5, 'PSEUDO_RIGID', bulk, solfec)

OUTPUT (solfec, 0.01)
RUN (solfec, gs, stop)
