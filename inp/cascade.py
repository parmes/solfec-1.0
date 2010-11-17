
def read_shapes (path, volid, surid):

  inp = open (path, 'r')
  nodes = []
  shapes = []

  lin = inp.readline ()
  lst = lin.split ()

  # read nodes
  if lst [0] != 'Points':
    print 'Invalid cascade formar'
    return None

  m = int (lst [1])

  for i in range (0, m):

    lin = inp.readline ()
    lst = lin.split ()

    if len (lst) != 4:
      print 'Invalid node format'

    nodes.append (float (lst [1]))
    nodes.append (float (lst [2]))
    nodes.append (float (lst [3]))

  # skip spaces
  lin = inp.readline ()
  lst = lin.split ()
  while lin != "" and len(lst) == 0:
    lin = inp.readline ()
    lst = lin.split ()

  # read elements
  if lst [0] == 'Elements':

    m = int (lst [1])

    for i in range (0, m):

      lin = inp.readline ()
      lst = lin.split ()
      verts = []

      for i in range (1, len (lst)):
	k = 3 * (int (lst [i]) - 1)
	x = nodes [k]
	y = nodes [k+1]
	z = nodes [k+2]
	verts.append (x)
	verts.append (y)
	verts.append (0)
	verts.append (x)
	verts.append (y)
	verts.append (z)

      shp = HULL (verts, volid, surid)
      shapes.append (shp)

  else:
    print 'Invalid format: no elements defined'

  inp.close ()

  return shapes

def create_balls (shapes, numbers, radius, leyers, thickness):

  volid = 2
  surid = 2
  balls = []
  obstacles = []
  xmin = 1E6
  ymin = 1E6
  xmax = -xmin
  ymax = -ymin
  zmin = 1E9
  for num in numbers:
    shp = shapes [num - 1]
    w = []
    for i in range (0, shp.nver):
      v = shp.vertex (i)
      if v[0] < xmin: xmin = v[0]
      if v[1] < ymin: ymin = v[1]
      if v[0] > xmax: xmax = v[0]
      if v[1] > ymax: ymax = v[1]
      if v[2] > 0.0:
	w.append (v[0])
	w.append (v[1])
	w.append (v[2])
	w.append (v[0])
	w.append (v[1])
	w.append (v[2] + thickness)
	if v[2] < zmin: zmin = v[2]
    obs = HULL (w, volid, surid)
    obstacles.append (obs)

  n = int ((xmax - xmin) / (2*radius))
  m = int ((ymax - ymin) / (2*radius))

  for i in range (0, n):
    for j in range (0, m):
      for k in range (1, leyers+1):
	x = xmin + 2 * radius * i
	y = ymin + 2 * radius * j
	z = zmin + 2 * radius * k
	shp = SPHERE ((x, y, z), radius, volid, surid)
	balls.append (shp)

  up = OVERLAPPING (obstacles, balls)
  return OVERLAPPING (shapes, up, 'NOT')

# main module
step = 0.005
stop = 60

solfec = SOLFEC ('DYNAMIC', step, 'out/cascade')
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3)
bulk = BULK_MATERIAL (solfec, 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3)
GRAVITY (solfec, (0, 0, -10))
gs = GAUSS_SEIDEL_SOLVER (1E-3, 10)

shapes = read_shapes ('inp/mesh/cascade.dat', 1, 1)

balls = create_balls (shapes, [86], 0.5, 5, 1)

for shp in shapes:
  BODY (solfec, 'OBSTACLE', shp, bulk)

bv = []
id = 1
for shp in balls:
  b = BODY (solfec, 'RIGID', shp, bulk, label = str(id))
  bv.append (b)
  id = id + 1

RUN (solfec, gs, stop)

if not VIEWER() and solfec.mode == 'READ':

  dur = DURATION (solfec)
  while solfec.time < dur [1]:

    print '\nTime: %g' % solfec.time

    for b in bv:
      disp = DISPLACEMENT (b, b.center)
      pos = (b.center[0] + disp[0],
	     b.center[1] + disp[1],
	     b.center[2] + disp[2])
      print 'Body (%s) position = (%g, %g, %g)' % (b.label, pos[0], pos[1], pos[2])

    FORWARD (solfec, 1)
