# linked chain

from math import sin
from math import cos

r = 0.01
l = 0.05
d = 0.50
nlnk = 20
step = 0.0007
stop = 2
ofrq = step

GEOMETRIC_EPSILON (1E-7)

def single_link (x, y, z, r, l):
  da = 3.14159 / 4.0
  ver = []
  for i in range (0, 8):
    x0 = r * sin (da * i)
    y0 = r * cos (da * i) 
    ver.append (x + x0)
    ver.append (y + y0)
    ver.append (z)
    ver.append (x + x0)
    ver.append (y + y0)
    ver.append (z + l)
  ver.append (x)
  ver.append (y)
  ver.append (z - r)
  ver.append (x)
  ver.append (y)
  ver.append (z + l + r)
  hul = HULL (ver, 1, 1)
  return hul

sol = SOLFEC ('DYNAMIC', step, 'out/links')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E9,
		      poisson = 0.25,
		      density = 500)

surfmat = SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.5)

eps = 2E-2
shp = single_link (0, 0, r, r, l)
bd0 = BODY (sol, 'RIGID', shp, bulk)
for i in range (1, nlnk):
 shp = single_link (0, 0, r + i*(l+2*r+eps), r, l)
 bd1 = BODY (sol, 'RIGID', shp, bulk)
 PUT_RIGID_LINK (bd0, bd1, (0, 0, i*(l+2*r+eps)-eps), (0, 0, i*(l+2*r+eps)))
 bd0 = bd1

table = HULL ([0, 0, -.05,
               0, d, 0.05,
	       d, d, 0.05,
	       d, 0, -.05,
               0, 0, -0.1,
               0, d, -0.1,
	       d, d, -0.1,
	       d, 0, -0.1], 1, 1)

TRANSLATE (table, (-0.5*d, -0.5*d, 0))
BODY (sol, 'OBSTACLE', table, bulk)

gs = GAUSS_SEIDEL_SOLVER (1E-4, 100, 1E-4)

GRAVITY (sol, (0, 0, -10))
OUTPUT (sol, ofrq)
RUN (sol, gs, stop)
