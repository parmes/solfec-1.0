from random import randint
from random import random
from random import seed

def make_rock (x, y, z, r, kind, material, solfec):
  m = randint (8, 64)
  points = []
  for n in range (m):
    points.append (x + r * (1.0 - random()) * 2.0)
    points.append (y + r * (1.0 - random()) * 2.0)
    points.append (z + r * (1.0 - random()) * 2.0)

  hull = HULL (points, 1, 1)
  BODY (solfec, kind, hull, material)
 
step = 0.001
stop = 3.0
n = 7
m = 7
l = 5

seed (1)
solfec = SOLFEC ('DYNAMIC', step, 'out/rocks')
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5)
bulk = BULK_MATERIAL (solfec, 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3)
GRAVITY (solfec, (0, 0, -1), 10)
gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000)

for i in range (n):
  for j in range (m):
    make_rock (i, j, 0, 1, 'OBSTACLE', bulk, solfec)

for i in range (2,n-1):
  for j in range (2,m-1):
    for k in range (2, 2+l):
      make_rock (i, j, k, 0.4, 'RIGID', bulk, solfec)

OUTPUT (solfec, 0.03)
RUN (solfec, gs, stop)
