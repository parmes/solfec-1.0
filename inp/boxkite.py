# box kite example
from math import sin
from math import cos
from math import sqrt

PI = 3.14159265358979323846 
FIG9 = 1
FIG11 = 2
SHEAR = 3
SEPARATION = 4

def gcore_brick (inr, outd, keyinw, keyoutw1, keyoutw3, keyh, height, hstep, type):

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

  cvx = CONVEX (vertices, faces, 1)

  scl = (lx, ly, lz)
  vec = (pnt [0] - 0.5*lx, pnt [1] - 0.5*ly, pnt [2])
  zet = (0, 0, 1)

  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)
  if zrot != 0.0: ROTATE (cvx, pnt, zet, zrot)

  BODY (solfec, 'RIGID', cvx, material)

def gcore_integral_key (pnt, l, a, b, h, material, solfec):

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

  BODY (solfec, 'RIGID', shape, material)

def gcore_outer_bricks_and_keys (loose_gap, integral_gap, material, solfec):

  dfac = 0.015
  outd = 0.4598
  height = 0.025
  hstep = 0.0098
  keyw = 0.0381
  keyh = 0.0381

  # outer bricks
  for i in range (3):
    for j in range (3):

      if i == j and i == 1: continue

      x = -(outd + dfac) + i * (outd + dfac)
      y = -(outd + dfac) + j * (outd + dfac)
      z = 0.0
	  
      cvx = gcore_brick (0.1315, outd, keyw, 0.05075, 0.05080, keyh, height, hstep, FIG9)
      TRANSLATE (cvx, (x, y, z))

      BODY (solfec, 'RIGID', cvx, material)

      cvx = gcore_brick (0.1315, outd, keyw, 0.05156, 0.05161, keyh, height, 0.0101, FIG11)
      zero = (0, 0, 0)
      yaxis =  (0, 1, 0)
      vec = (x, y, z + (2*height-hstep))
      ROTATE (cvx, zero, yaxis, 180)
      TRANSLATE (cvx, vec)

      BODY (solfec , 'RIGID', cvx, material)

  # loose keys
  lx = keyw - 2.0*loose_gap
  ly = (2.0 * keyh + dfac) - 2.0*loose_gap
  lz = 2*height - hstep

  for i in range (3):
    for j in range (2):

      pnt = (-(outd + dfac) + i * (outd + dfac), -0.5*(outd + dfac) + j * (outd + dfac), 0.0)
      gcore_loose_key (pnt, lx, ly, lz, 0, material, solfec)

  for i in range (2):
    for j in range (3):

      pnt = (-0.5*(outd + dfac) + i * (outd + dfac), -(outd + dfac) + j * (outd + dfac), 0.0 )
      gcore_loose_key (pnt, lx, ly, lz, 90, material, solfec)

  # integral keys
  c = sqrt(2.0) * (outd + dfac)
  l = c - (outd + 2.0*dfac)
  a =  keyw - 2.0*integral_gap
  b = keyh + dfac - integral_gap

  for i in range (2):
    for j in range (2):

      pnt = (-0.5*(outd + dfac) + i * (outd + dfac), -0.5*(outd + dfac) + j * (outd + dfac), 0.0)
      gcore_integral_key (pnt, l, a, b, lz, material, solfec)

def gcore_base (material, solfec):

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
  height = 0.025
  hstep = 0.0098
  lz = 2*height - hstep
  shape = []

  cvx = CONVEX (vertices, faces, 3)
  scl = (1.50,  1.50,  0.1)
  vec = (-0.75, -0.75, -0.1)
  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)
  shape.append (cvx)

  cvx = CONVEX (vertices, faces, 3)
  len = 3 * outd + 4 * dfac
  h = (1.50 - len) / 2
  scl = (1.50,  h,  lz)
  vec = (-0.75, -0.75, 0)
  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)
  shape.append (cvx)

  cvx = CONVEX (vertices, faces, 3)
  len = 3 * outd + 4 * dfac
  h = (1.50 - len) / 2
  scl = ( 1.50,  h,  lz)
  vec = (-0.75, 0.75 - h, 0)
  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)
  shape.append (cvx)

  cvx = CONVEX (vertices, faces, 3)
  len = 3 * outd + 4 * dfac
  h = (1.50 - len) / 2
  scl = ( h,  1.50 - 2*h,  lz)
  vec = (-0.75, -0.75 + h, 0)
  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)
  shape.append (cvx)

  cvx = CONVEX (vertices, faces, 3)
  len = 3 * outd + 4 * dfac
  h = (1.50 - len) / 2
  scl = (h,  1.50 - 2*h,  lz)
  vec = (0.75-h, -0.75 + h, 0)
  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)
  shape.append (cvx)

  BODY (solfec, 'OBSTACLE', shape, material)

def gcore_central_bricks (loangle, hiangle, keywayangle, kind, tms, material, solfec):

  pt = (0, 0, 0)
  zero  = (0, 0, 0)
  ox = (1, 0, 0)
  oy = (0, 1, 0)
  oz = (0, 0, 1)

  outd = 0.4598
  height = 0.025
  hstep = 0.0098
  keyw = 0.0381
  keyh = 0.0381

  out = [None, None, None, None]

  # loangle rotation of low slice
  nl = ROTATE (ox, zero, oz, loangle)
  cvx = gcore_brick (0.1315, outd, keyw, 0.05075, 0.05080, keyh, height, hstep, FIG9)
  ROTATE (cvx, zero, oz, keywayangle)
  (out[0], out[1]) = SPLIT (cvx, pt, nl)

  # hiangle rotation of high slice
  nl = ROTATE (ox, zero, oz, hiangle)
  cvx = gcore_brick (0.1315, outd, keyw, 0.05156, 0.05161, keyh, height, 0.0101, FIG11)
  vec = (0, 0, (2*height-hstep))
  ROTATE (cvx, zero, oy, 180)
  ROTATE (cvx, zero, oz, keywayangle)
  TRANSLATE (cvx, vec)
  (out[2], out[3]) = SPLIT (cvx, pt, nl)

  label = ['HALF0', 'HALF1', 'HALF2', 'HALF3']

  for i in range (4):

    bod = BODY (solfec, 'RIGID', out [i], material, label [i])

    if i == 2:

      if kind == SHEAR:
	dir = ROTATE (nl, zero, oz, 90)
      else: dir = nl

      dir = SCALE (dir, (-1, -1, -1))
      FORCE (bod, 'SPATIAL', MASS_CENTER (bod), dir, tms)

    elif i == 3:

      if kind == SHEAR:
	dir = ROTATE (nl, zero, oz, 90)
      else: dir = nl

      FORCE (bod, 'SPATIAL', MASS_CENTER (bod), dir, tms)
    

def box_kite_create (loose_gap, integral_gap, high_angle, low_angle, keyway_angle, kind, tms, material, solfec):

  gcore_base (material, solfec)
  gcore_outer_bricks_and_keys (loose_gap, integral_gap, material, solfec)
  gcore_central_bricks (low_angle, high_angle, keyway_angle, kind, tms, material, solfec)

def box_kite_test (case_name, time_step, duration, friction_coef,
                  load_case, load_hist, solver,
                  loose_gap, integral_gap, high_angle, low_angle, keyway_angle):

  solfec = SOLFEC ('DYNAMIC', time_step, 'out/boxkite/' + case_name)
  surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = friction_coef)
  bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 5E6, poisson = 0.25, density = 1E3)
  GRAVITY (solfec, (0, 0, -10))
  box_kite_create (loose_gap, integral_gap,  high_angle, low_angle, keyway_angle, load_case, load_hist, bulkmat, solfec)
  RUN (solfec, solver, duration)

  return solfec

def box_kite_print_history (solfec, case_name, load_case, high_angle):

  zero = (0, 0, 0)
  ox = (1, 0, 0)
  oz = (0, 0, 1)

  if solfec.mode == 'WRITE': return

  bod0 = BYLABEL (solfec, 'BODY', 'HALF2')
  bod1 = BYLABEL (solfec, 'BODY', 'HALF3')
  pnt0 = MASS_CENTER (bod0)
  pnt1 = MASS_CENTER (bod1)
  dur = DURATION (solfec)

  th = HISTORY (solfec, [(bod0, pnt0, 'DX'), (bod0, pnt0, 'DY'),
                         (bod1, pnt1, 'DX'), (bod1, pnt1, 'DY')], dur[0], dur[1])

  if load_case == SHEAR:
    dir = ROTATE (ox, zero, oz, high_angle + 90)
  else:
    dir = ROTATE (ox, zero, oz, high_angle)

  time = th [0]
  disp = []
  for i in range (0, len (th[0])):
    disp.append (dir [0] * (th[1][i] - th[3][i]) + dir [1] * (th[2][i] - th[4][i]))
    disp [i] = abs (disp [i])

  file = open ('out/boxkite/' + case_name + '/disp.txt', 'w')

  for i in range (0, len (time)):
    file.write (str(time [i])  + '   ' +  str (disp [i]) + '\n')

### main module ###

save = []

GEOMETRIC_EPSILON (1E-6)

step = 0.001
duration = 0.01
load_hist = TIME_SERIES ([0, 10, 1, 10])
friction = 0.0

solver = GAUSS_SEIDEL_SOLVER (1E0, 50, 1E-6)

allcases = [
# large clearance
            ('29-LT', SHEAR, 0.0015, 0.0013, 0, 0, 0),
            ('29-LN', SEPARATION, 0.0015, 0.0013, 0, 0, 0), 
            ('31-LT', SHEAR, 0.0015, 0.0013, 45, 0, 45),
            ('31-LN', SEPARATION, 0.0015, 0.0013, 45, 0, 45),
            ('33-LT', SHEAR, 0.0015, 0.0013,  0,   0, 45),
            ('33-LN', SEPARATION, 0.0015, 0.0013,  0,   0, 45),
            ('39-LT', SHEAR, 0.0015, 0.0013, 45,   0,  0),
            ('39-LN', SEPARATION, 0.0015, 0.0013, 45,   0,  0),
            ('44-LT', SHEAR, 0.0015, 0.0013, 45, 135,  0),
            ('44-LN', SEPARATION, 0.0015, 0.0013, 45, 135,  0),
# small clearance
            ('46-ST', SHEAR, 0.0003, 0.0002, 0, 0, 0),
            ('46-SN', SEPARATION, 0.0003, 0.0002, 0, 0, 0), 
            ('48-ST', SHEAR, 0.0003, 0.0002, 45, 0, 45),
            ('48-SN', SEPARATION, 0.0003, 0.0002, 45, 0, 45),
            ('50-ST', SHEAR, 0.0015, 0.0002,  0,   0, 45),
            ('50-SN', SEPARATION, 0.0015, 0.0002,  0,   0, 45),
            ('56-ST', SHEAR, 0.0003, 0.0002, 45,   0,  0),
            ('56-SN', SEPARATION, 0.0003, 0.0002, 45,   0,  0),
            ('61-ST', SHEAR, 0.0003, 0.0002, 45, 135,  0),
            ('61-SN', SEPARATION, 0.0003, 0.0002, 45, 135,  0)
	    ]

for case in allcases:
  sol = box_kite_test (case [0], step, duration, friction, case [1], load_hist, solver, case [2], case [3], case [4], case [5], case [6])
  save.append (sol)
  box_kite_print_history (sol, case [0], case [1], case [4])

def one_test (number):
  case = allcases [number]
  sol = box_kite_test (case [0], step, duration, friction, case [1], load_hist, solver, case [2], case [3], case [4], case [5], case [6])
  save.append (sol)
  box_kite_print_history (sol, case [0], case [1], case [4])

#one_test (0)
