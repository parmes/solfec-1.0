# 2D puzzle example

KINEM = 'PSEUDO_RIGID'
SCHEME = 'DEF_LIM'

def create_big_piece (x, y, z, solfec, material):

  pattern = [[1, 1, 1, 1, 0, 1, 1, 1, 1],
             [1, 0, 1, 0, 0, 0, 1, 0, 1],
	     [0, 0, 1, 1, 1, 1, 1, 0, 0],
	     [1, 0, 1, 0, 0, 0, 1, 0, 1],
	     [1, 1, 1, 1, 0, 1, 1, 1, 1]]

  box = HULL ([x, y, z,
               x+1, y, z,
	       x+1, y+1, z,
	       x, y+1, z,
               x, y, z+1,
               x+1, y, z+1,
	       x+1, y+1, z+1,
	       x, y+1, z+1], 1, 1)


  shape = []
  x0 = x
  for row in pattern:
    x = x0
    for flag in row:
      if flag:
	shp = COPY (box)
	TRANSLATE (shp, (x, y, z))
	shape.append (shp)
      x = x + 1
    y = y + 1

  bod = BODY (solfec, KINEM, shape, material)
  bod.scheme = SCHEME

def create_vertical_piece (x, y, z, solfec, material):

  pattern = [[1, 1, 1],
             [0, 1, 0],
             [0, 1, 0],
             [1, 1, 1]]

  box = HULL ([x, y, z,
               x+1, y, z,
	       x+1, y+1, z,
	       x, y+1, z,
               x, y, z+1,
               x+1, y, z+1,
	       x+1, y+1, z+1,
	       x, y+1, z+1], 2, 2)

  eps = 1E-3
  TRANSLATE (box, (-x, -y, -z))
  veps = (1-2*eps, 1-2*eps, 1)
  SCALE (box, veps)
  TRANSLATE (box, (x+eps, y+eps, z))

  shape = []
  x0 = x
  for row in pattern:
    x = x0
    for flag in row:
      if flag:
	shp = COPY (box)
	TRANSLATE (shp, (x, y, z))
	shape.append (shp)
      x = x + 1
    y = y + 1

  bod = BODY (solfec, KINEM, shape, material)
  bod.scheme = SCHEME

def create_horizontal_piece (x, y, z, solfec, material):

  pattern = [[1, 0, 0, 1],
             [1, 1, 1, 1],
             [1, 0, 0, 1]]

  box = HULL ([x, y, z,
               x+1, y, z,
	       x+1, y+1, z,
	       x, y+1, z,
               x, y, z+1,
               x+1, y, z+1,
	       x+1, y+1, z+1,
	       x, y+1, z+1], 2, 2)

  eps = 1E-3
  TRANSLATE (box, (-x, -y, -z))
  veps = (1-2*eps, 1-2*eps, 1)
  SCALE (box, veps)
  TRANSLATE (box, (x+eps, y+eps, z))

  shape = []
  x0 = x
  for row in pattern:
    x = x0
    for flag in row:
      if flag:
	shp = COPY (box)
	TRANSLATE (shp, (x, y, z))
	shape.append (shp)
      x = x + 1
    y = y + 1

  bod = BODY (solfec, KINEM, shape, material)
  bod.scheme = SCHEME

def create_module (x, y, z, solfec, mat1, mat2):
  create_big_piece (x, y, z, solfec, mat1)
  if y > 0: create_vertical_piece (x+1.5, y-1.0, z, solfec, mat2)
  if x > 0: create_horizontal_piece (x-1, y+0.5, z, solfec, mat2)

def create_forced_cube (x, y, z, wx, wy, wz, solfec, material):

  pattern = [[1, 0, 0, 1],
             [1, 1, 1, 1],
             [1, 0, 0, 1]]

  shp = HULL ([x, y, z,
               x+wx, y, z,
	       x+wx, y+wy, z,
	       x, y+wy, z,
               x, y, z+wz,
               x+wx, y, z+wz,
	       x+wx, y+wy, z+wz,
	       x, y+wy, z+wz], 3, 3)

  point = (x+0.5*wx, y+0.5*wy, z+0.5*wz)
  bod = BODY (solfec, 'RIGID', shp, material)
  FORCE (bod, 'SPATIAL', point, (0, 1, 0), 1E8)

def create_obstacle_cube (x, y, z, wx, wy, wz, solfec, material):

  pattern = [[1, 0, 0, 1],
             [1, 1, 1, 1],
             [1, 0, 0, 1]]

  shp = HULL ([x, y, z,
               x+wx, y, z,
	       x+wx, y+wy, z,
	       x, y+wy, z,
               x, y, z+wz,
               x+wx, y, z+wz,
	       x+wx, y+wy, z+wz,
	       x, y+wy, z+wz], 4, 4)

  bod = BODY (solfec, 'OBSTACLE', shp, material)

# main module
step = 0.001
stop = 1

solfec = SOLFEC ('DYNAMIC', step, 'out/puzzle')
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.7, restitution = 0.0)
bulk1 = BULK_MATERIAL (solfec, 'KIRCHHOFF', young = 15E8, poisson = 0.25, density = 1.8E3)
bulk2 = BULK_MATERIAL (solfec, 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3)
GRAVITY (solfec, (0, 0, 0))

#sv = GAUSS_SEIDEL_SOLVER (1E-3, 10, 1E-3)
sv = NEWTON_SOLVER (1E-5, 100)

n_i = 8
n_j = 4

for i in range (0, n_i):
  for j in range (0, n_j):
    create_module (4.5 * i, 2.5 * j, 0, solfec, bulk1, bulk2)

create_forced_cube (1, -2, 0, 2, 2, 1, solfec, bulk2)
create_forced_cube (9 * n_i - 3, -2, 0, 2, 2, 1, solfec, bulk2)
create_obstacle_cube (4.5 * n_i - 1, 5*n_j, 0, 2, 2, 1, solfec, bulk2)

RUN (solfec, sv, stop)
