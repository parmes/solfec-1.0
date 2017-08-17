# SOLFEC parallel scaling example
# -------------------------------------------------------
# Rotating drum filled with ellipsoids

from math import sqrt
from math import sin
from math import cos

def WHEEL (x, y, z, r, t, v, s):
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

def BLOCK (x, y, z, wx, wy, wz, v, s):
  points = [x-.5*wx, y-.5*wy, z-.5*wz,
            x+.5*wx, y-.5*wy, z-.5*wz,
            x+.5*wx, y+.5*wy, z-.5*wz,
            x-.5*wx, y+.5*wy, z-.5*wz,
            x-.5*wx, y-.5*wy, z+.5*wz,
            x+.5*wx, y-.5*wy, z+.5*wz,
            x+.5*wx, y+.5*wy, z+.5*wz,
            x-.5*wx, y+.5*wy, z+.5*wz]
  return HULL (points, v, s)

def ELLIPS (rx, ry, rz, x0, y0, z, wx, wy, v, s):
  global iell, nell, sol, mat
  a = x0 + .5*wx
  b = y0 + .5*wy
  x = x0 - .5*wx
  while x < a:
    y = y0 - .5*wy
    while y < b:
      if iell < nell:
        if RANK() == 0:
	  shp = ELLIP ((x, y, z), (rx, ry, rz), v, s)
	  BODY (sol, 'RIGID', shp, mat)
	iell = iell + 1
      y = y + 2.*ry
    x = x + 2.*rx

step = 0.001
outi = 0.05
stop = 5.0
nell = 100
iell = 0

sol = SOLFEC ('DYNAMIC', step, 'out/rotating-drum')
mat = BULK_MATERIAL (sol, young = 1E6, poisson = 0.25, density = 100)
SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.4)
sv = GAUSS_SEIDEL_SOLVER (1E-3, 100)
GRAVITY (sol, (0, 0, -10))
OUTPUT (sol, outi)

# drum
pip = PIPE ((0, 0, 0), (0, 0.5, 0), 1, 0.05, 1, 32, 1, 1, [1, 2, 3, 4])
bl1 = BLOCK(-0.9, 0.25, 0, 0.2, 0.5, 0.2, 1, 1)
bl2 = BLOCK (0.9, 0.25, 0, 0.2, 0.5, 0.2, 1, 1)
bod = BODY (sol, 'OBSTACLE', [pip, bl1, bl2], mat)

# side walls
wh1 = WHEEL (0, -0.002, 0, 1.05, -0.05, 1, 1)
BODY (sol, 'OBSTACLE', wh1, mat)
wh2 = WHEEL (0, 0.502, 0, 1.05, 0.05, 1, 1)
BODY (sol, 'OBSTACLE', wh2, mat)

# ellipsoid radii
rr = (.1/nell)**(1./3.)
rx = .9*rr
ry = .5*rr
rz = rr
# .5*acc*t^2 = 2*rz+eps --> t = sqrt(2*(2*rr+eps)/acc)
dt = sqrt(2*(2*rz+.1*rz)/10)

rotating = False
def callback (sol):
  if iell == nell and not rotating:
    INITIAL_VELOCITY (bod, (0, 0, 0), (0, 1, 0))
  else: ELLIPS (rx, ry, rz, 0, 0.25, 0.25, 1.4, 0.4, 2, 2)
  return 1

CALLBACK (sol, dt, sol, callback)
RUN (sol, sv, stop)
