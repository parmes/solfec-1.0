# simple core model
from math import sin
from math import cos
from math import sqrt
from solfec import *

# constants
PI = 3.14159265358979323846 
FIG9 = 1
FIG11 = 2

def gcore_brick_half (inr, outd, keyinw, keyoutw1, keyoutw3, keyh, height, hstep, type):

  faces = [4, 0, 1, 5, 4, 0,
           4, 1, 2, 6, 5, 0,
           4, 2, 3, 7, 6, 0,
           4, 3, 0, 4, 7, 0,
           4, 0, 3, 2, 1, 0,
           4, 4, 5, 6, 7, 0]


  step = 2.0 * PI / 16.0
  angle = PI / 2.0 + step / 2.0

  zero = (0, 0, 0)
  zet = (0, 0, 1)
  toph =  outd * 0.5
  outr = toph / sin (angle)
  keyoutw2 = keyoutw1 + (keyoutw3 - keyoutw1) * (toph - inr) / (toph - inr - keyh)

  a = [keyoutw1/2.0, inr, 0,
       keyoutw2/2.0, toph-keyh, 0,
      -keyoutw2/2.0, toph-keyh, 0,
      -keyoutw1/2.0, inr, 0,
       keyoutw1/2.0, inr, height,
       keyoutw2/2.0, toph-keyh, height,
      -keyoutw2/2.0, toph-keyh, height,
      -keyoutw1/2.0, inr, height]

  set = []
  list = []
  c1 = None
  c2 = None

  for i in range (8):

    cvy = CONVEX (a, faces, 0);
    ROTATE (cvy, zero, zet, i*45);

    if i == 1: c1 = cvy

    if type == FIG9:
      if i % 2 != 0: set.append (cvy)
    elif type == FIG11:
      if i % 2 == 0: set.append (cvy)

    list.append (cvy)


  b = [keyoutw2/2.0, toph-keyh, 0,
       keyoutw3/2.0, toph, 0,
       keyinw/2.0, toph, 0,
       keyinw/2.0, toph-keyh, 0,
       keyoutw2/2.0, toph-keyh, height,
       keyoutw3/2.0, toph, height,
       keyinw/2.0, toph, height,
       keyinw/2.0, toph-keyh, height]

  c = [-keyinw/2.0, toph-keyh, 0,
       -keyinw/2.0, toph, 0,
       -keyoutw3/2.0, toph, 0,
       -keyoutw2/2.0, toph-keyh, 0,
       -keyinw/2.0, toph-keyh, height,
       -keyinw/2.0, toph, height,
       -keyoutw3/2.0, toph, height,
       -keyoutw2/2.0, toph-keyh, height]


  for i in range (8):

    cvx = CONVEX (b, faces, 0)
    cvy = CONVEX (c, faces, 0)

    ROTATE (cvx, zero, zet, i*45)
    ROTATE (cvy, zero, zet, i*45)

    if i == 1: c2 = cvx

    if type == FIG9:
      if i % 2 != 0:
	set.append (cvx)
	set.append (cvy)
    elif type == FIG11:
      if i % 2 == 0:
	set.append (cvx)
	set.append (cvy)

    list.append (cvx)
    list.append (cvy)


  d = [a [9], a [10], 0,
       c [6], c [7], 0,
       outr * cos (angle), outr * sin (angle), 0,
       outr * cos (angle+step), outr * sin (angle+step), 0,
       c2.vertex(1)[0], c2.vertex(1)[1], 0,
       c1.vertex(0)[0], c1.vertex(0)[1], 0,
       a [9], a [10], height,
       c [6], c [7], height,
       outr * cos (angle), outr * sin (angle), height,
       outr * cos (angle+step), outr * sin (angle+step), height,
       c2.vertex(1)[0], c2.vertex(1)[1], height,
       c1.vertex(0)[0], c1.vertex(0)[1], height]

  g = [4, 0, 1, 7, 6, 0,
       4, 1, 2, 8, 7, 0,
       4, 2, 3, 9, 8, 0,
       4, 3, 4, 10, 9, 0,
       4, 4, 5, 11, 10, 0,
       4, 5, 0, 6, 11, 0,
       6, 0, 5, 4, 3, 2, 1, 0,
       6, 6, 7, 8, 9, 10, 11, 0]

  for i in range (8):

    cvy = CONVEX (d, g, 0)
    ROTATE (cvy, zero, zet, i*45);

    if type == FIG9: set.append (cvy)

    list.append (cvy)

  scal = (1, 1, (height - hstep) / height)

  for item in set:
    SCALE (item, scal)

  return list

def gcore_loose_key (pnt, lx, ly, lz, zrot, material, solfec):

  vertices  = [1, 0, 0,
               1, 1, 0,
               0, 1, 0,
               0, 0, 0,
               1, 0, 1,
               1, 1, 1,
               0, 1, 1,
               0, 0, 1]

  faces = [4, 0, 1, 5, 4, 1,
           4, 1, 2, 6, 5, 1,
           4, 2, 3, 7, 6, 1,
           4, 3, 0, 4, 7, 1,
           4, 0, 3, 2, 1, 1,
           4, 4, 5, 6, 7, 1]

  hex = HEX (vertices, 1, 1, 3, 1, [1, 1, 1, 1, 1, 1])

  scl = (lx, ly, lz)
  vec = (pnt [0] - 0.5*lx, pnt [1] - 0.5*ly, pnt [2])
  zet = (0, 0, 1)

  SCALE (hex, scl)
  TRANSLATE (hex, vec)
  if zrot != 0.0: ROTATE (hex, pnt, zet, zrot)

  BODY (solfec, 'RIGID', hex, material)

def gcore_integral_key (pnt, l, a, b, h, material, solfec, kinem_kind, integ_scheme):

  vertices = [1, 0, 0,
              1, 1, 0,
              0, 1, 0,
              0, 0, 0,
              1, 0, 1,
              1, 1, 1,
              0, 1, 1,
              0, 0, 1]

  faces = [4, 0, 1, 5, 4, 2,
           4, 1, 2, 6, 5, 2,
           4, 2, 3, 7, 6, 2,
           4, 3, 0, 4, 7, 2,
           4, 0, 3, 2, 1, 2,
           4, 4, 5, 6, 7, 2]

  cv1 = CONVEX (vertices, faces, 2)
  scl = (l, l, h)
  vec = (pnt [0] - 0.5*l, pnt [1] - 0.5*l, pnt [2])
  SCALE (cv1, scl)
  TRANSLATE (cv1, vec)

  cv2 = CONVEX (vertices, faces, 2)
  scl = (a, b, h)
  vec = (pnt [0] - 0.5*a, pnt [1] + 0.5*l, pnt [2])
  SCALE (cv2, scl)
  TRANSLATE (cv2, vec)

  cv3 = CONVEX (vertices, faces, 2)
  scl = (a, b, h)
  vec = (pnt [0] - 0.5*a, pnt [1] - b - 0.5*l, pnt [2])
  SCALE (cv3, scl)
  TRANSLATE (cv3, vec)

  cv4 = CONVEX (vertices, faces, 2)
  scl = (b, a, h)
  vec = (pnt [0] - b - 0.5*l, pnt [1] - 0.5*a, pnt [2])
  SCALE (cv4, scl)
  TRANSLATE (cv4, vec)

  cv5 = CONVEX (vertices, faces, 2)
  scl = (b, a, h)
  vec = (pnt [0] + 0.5*l, pnt [1] - 0.5*a, pnt [2])
  SCALE (cv5, scl)
  TRANSLATE (cv5, vec)

  zet = (0, 0, 1)
  shape = [cv1, cv2, cv3, cv4, cv5]
  ROTATE (shape, pnt, zet, 45.0)

  if kinem_kind != 'RIGID':
    bod = BODY (solfec, 'PSEUDO_RIGID', shape, material)
  else: bod = BODY (solfec, 'RIGID', shape, material)
  bod.scheme = integ_scheme

def gcore_brick (x, y, z):

  dfac = 0.015
  outd = 0.4598
  height = 0.225
  hstep = 0.0098
  keyw = 0.0381
  keyh = 0.0381

  cvx = gcore_brick_half (0.1315, outd, keyw, 0.05156, 0.05161, keyh, height, 0.0101, FIG11)
  zero = (0, 0, 0)
  yaxis =  (0, 1, 0)
  vec = (x, y, z + height)
  ROTATE (cvx, zero, yaxis, 180)
  TRANSLATE (cvx, vec)

  vec = (x, y, z + height)
  cvy = gcore_brick_half (0.1315, outd, keyw, 0.05075, 0.05080, keyh, height, hstep, FIG9)
  TRANSLATE (cvy, vec)

  cvx.extend (cvy)

  return cvx

def gcore_bricks_and_keys (loose_gap, integral_gap, material, solfec, kinem_kind, integ_scheme, N_BRICKS, M_BRICKS, N_LAYERS):

  dfac = 0.015
  outd = 0.4598
  height = 0.225
  hstep = 0.0098
  keyw = 0.0381
  keyh = 0.0381

  for k in range (N_LAYERS):

    z = k * (2*height - hstep)

    # outer bricks
    for i in range (N_BRICKS):
      for j in range (M_BRICKS):

	x = -(outd + dfac) + i * (outd + dfac)
	y = -(outd + dfac) + j * (outd + dfac)

	shp = gcore_brick (x, y, z)
	if kinem_kind == 'FINITE_ELEMENT':
          msh = PIPE ((x, y, z), (0, 0, 0.45), 0.125, 0.13, 3, 8, 2, 0, [0, 0, 0, 0])
	  bod = BODY (solfec , kinem_kind, shp, material, mesh = msh, form = 'BC')
	else: bod = BODY (solfec , kinem_kind, shp, material)
	bod.scheme = integ_scheme

    # loose keys
    lx = keyw - 2.0*loose_gap
    ly = (2.0 * keyh + dfac) - 2.0*loose_gap
    lz = 2*height - hstep

    for i in range (N_BRICKS):
      for j in range (M_BRICKS-1):

	pnt = (-(outd + dfac) + i * (outd + dfac), -0.5*(outd + dfac) + j * (outd + dfac), z)
	gcore_loose_key (pnt, lx, ly, lz, 0, material, solfec)

    for i in range (N_BRICKS-1):
      for j in range (M_BRICKS):

	pnt = (-0.5*(outd + dfac) + i * (outd + dfac), -(outd + dfac) + j * (outd + dfac), z)
	gcore_loose_key (pnt, lx, ly, lz, 90, material, solfec)

    # integral keys
    c = sqrt(2.0) * (outd + dfac)
    l = c - (outd + 2.0*dfac)
    a =  keyw - 2.0*integral_gap
    b = keyh + dfac - integral_gap

    for i in range (N_BRICKS-1):
      for j in range (M_BRICKS-1):

	pnt = (-0.5*(outd + dfac) + i * (outd + dfac), -0.5*(outd + dfac) + j * (outd + dfac), z)
	gcore_integral_key (pnt, l, a, b, lz, material, solfec, kinem_kind, integ_scheme)

def gcore_base (material, solfec, shake_base, N_BRICKS, M_BRICKS, N_LAYERS):

  vertices = [1, 0, 0,
              1, 1, 0,
              0, 1, 0,
              0, 0, 0,
              1, 0, 1,
              1, 1, 1,
              0, 1, 1,
              0, 0, 1]

  faces = [4, 0, 1, 5, 4, 3,
           4, 1, 2, 6, 5, 3,
           4, 2, 3, 7, 6, 3,
           4, 3, 0, 4, 7, 3,
           4, 0, 3, 2, 1, 3,
           4, 4, 5, 6, 7, 3]

  dfac = 0.015
  outd = 0.4598
  height = 0.225
  hstep = 0.0098
  margin = 0.05
  thick = 0.1
  sx = - (outd + dfac) - outd / 2 - margin - thick
  sy = - (outd + dfac) - outd / 2 - margin - thick
  lx = N_BRICKS * outd + (N_BRICKS-1) * dfac + 2 * (margin + thick)
  ly = M_BRICKS * outd + (M_BRICKS-1) * dfac + 2 * (margin + thick)
  lz = N_LAYERS * (2*height - hstep)
  shape = []

  cvx = CONVEX (vertices, faces, 3)
  scl = (lx,  thick,  lz)
  vec = (sx, sy, 0)
  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)
  shape.append (cvx)

  cvx = CONVEX (vertices, faces, 3)
  scl = (lx,  thick,  lz)
  vec = (sx, sy + ly - thick, 0)
  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)
  shape.append (cvx)

  cvx = CONVEX (vertices, faces, 3)
  scl = (thick,  ly - 2*thick,  lz)
  vec = (sx, sy + thick, 0)
  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)
  shape.append (cvx)

  cvx = CONVEX (vertices, faces, 3)
  scl = (thick,  ly - 2*thick,  lz)
  vec = (sx + lx - thick, sy + thick, 0)
  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)
  shape.append (cvx)

  wall = BODY (solfec, 'OBSTACLE', shape, material)

  for i in range (N_BRICKS):
    for j in range (M_BRICKS):

      x = -(outd + dfac) + i * (outd + dfac)
      y = -(outd + dfac) + j * (outd + dfac)
      sx = x - outd * 0.5 - dfac * 0.25
      sy = y - outd * 0.5 - dfac * 0.25
      lx = outd + dfac * 0.5
      ly = outd + dfac * 0.5

      shape = []
      cvx = CONVEX (vertices, faces, 3)
      scl = (lx,  ly,  thick)
      vec = (sx, sy, -thick)
      SCALE (cvx, scl)
      TRANSLATE (cvx, vec)
      shape.append (cvx)

      keyw = 0.0381
      keyh = 0.0381
      vec = (x, y, - 0.1 * height + hstep)
      shp = gcore_brick_half (0.1315, outd, keyw, 0.05075, 0.05080, keyh, 0.1 * height, hstep, FIG9)
      TRANSLATE (shp, vec)
      shape += shp

      if shake_base == 'TRUE':
	base = BODY (solfec, 'RIGID', shape, material)
	CONTACT_EXCLUDE_BODIES (base, wall)
	acc = TIME_SERIES ('inp/cores/inc/acc-0.dat')
	FIX_DIRECTION (base, (sx, sy, -thick), (0, 0, 1))
	FIX_DIRECTION (base, (sx + lx, sy, -thick), (0, 0, 1))
	FIX_DIRECTION (base, (sx, sy + ly, -thick), (0, 0, 1))
	FIX_DIRECTION (base, (sx + lx, sy + ly, -thick), (0, 0, 1))
	FIX_DIRECTION (base, (sx, sy, -thick), (0, 1, 0))
	FIX_DIRECTION (base, (sx, sy + ly, -thick), (0, 1, 0))
	SET_ACCELERATION (base, (sx + 0.5 * lx, sy + 0.5 * ly, - thick), (1, 0, 0), acc)
      else:
	BODY (solfec, 'OBSTACLE', shape, material)

def simple_core_create (loose_gap, integral_gap, material, solfec, kinem_kind, integ_scheme, shake_base, N_BRICKS, M_BRICKS, N_LAYERS):

  gcore_base (material, solfec, shake_base, N_BRICKS, M_BRICKS, N_LAYERS)
  gcore_bricks_and_keys (loose_gap, integral_gap, material, solfec, kinem_kind, integ_scheme, N_BRICKS, M_BRICKS, N_LAYERS)
