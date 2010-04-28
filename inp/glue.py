# gluing example

def cube (x, y, z, a, b, c, sur, vol):

  nodes = [0, 0, 0,
	   a, 0, 0,
	   a, b, 0,
	   0, b, 0,
	   0, 0, c,
	   a, 0, c,
	   a, b, c,
	   0, b, c]

  shp = HEX (nodes, 1, 1, 1, vol, [sur, sur, sur, sur, sur, sur])

  TRANSLATE (shp, (x, y, z))

  return shp

def glued_shape (x, y, z, material, solfec):
  list = []

  shp = cube (x, y, z, 1, 1, 1, 2, 2)
  b1 = BODY (solfec, 'RIGID', shp, material)
  list.append (b1)
  shp = cube (x+1, y, z, 1, 1, 1, 2, 2)
  b2 = BODY (solfec, 'RIGID', shp, material)
  list.append (b2)
  GLUE_POINTS (b1, b2, (x+1, y, z))
  GLUE_POINTS (b1, b2, (x+1, y+1, z))
  GLUE_POINTS (b1, b2, (x+1, y, z+1))
  GLUE_POINTS (b1, b2, (x+1, y+1, z+1))
  shp = cube (x-1, y, z, 1, 1, 1, 2, 2)
  b2 = BODY (solfec, 'RIGID', shp, material)
  list.append (b2)
  GLUE_POINTS (b1, b2, (x, y, z))
  GLUE_POINTS (b1, b2, (x, y+1, z))
  GLUE_POINTS (b1, b2, (x, y, z+1))
  GLUE_POINTS (b1, b2, (x, y+1, z+1))
  shp = cube (x, y+1, z, 1, 1, 1, 2, 2)
  b2 = BODY (solfec, 'RIGID', shp, material)
  list.append (b2)
  GLUE_POINTS (b1, b2, (x, y+1, z))
  GLUE_POINTS (b1, b2, (x+1, y+1, z))
  GLUE_POINTS (b1, b2, (x, y+1, z+1))
  GLUE_POINTS (b1, b2, (x+1, y+1, z+1))
  shp = cube (x, y-1, z, 1, 1, 1, 2, 2)
  b2 = BODY (solfec, 'RIGID', shp, material)
  list.append (b2)
  GLUE_POINTS (b1, b2, (x, y, z))
  GLUE_POINTS (b1, b2, (x+1, y, z))
  GLUE_POINTS (b1, b2, (x, y, z+1))
  GLUE_POINTS (b1, b2, (x+1, y, z+1))
  shp = cube (x, y, z+1, 1, 1, 1, 2, 2)
  b2 = BODY (solfec, 'RIGID', shp, material)
  list.append (b2)
  GLUE_POINTS (b1, b2, (x, y, z+1))
  GLUE_POINTS (b1, b2, (x+1, y, z+1))
  GLUE_POINTS (b1, b2, (x, y+1, z+1))
  GLUE_POINTS (b1, b2, (x+1, y+1, z+1))
  shp = cube (x, y, z-1, 1, 1, 1, 2, 2)
  b2 = BODY (solfec, 'RIGID', shp, material)
  list.append (b2)
  GLUE_POINTS (b1, b2, (x, y, z))
  GLUE_POINTS (b1, b2, (x+1, y, z))
  GLUE_POINTS (b1, b2, (x, y+1, z))
  GLUE_POINTS (b1, b2, (x+1, y+1, z))
  for a in list:
    for b in list:
      CONTACT_EXCLUDE_BODIES (a, b)

def glued_scene (material, solfec):

  # create an obstacle base
  shp = cube (0, 0, -1, 1, 1, 1, 1, 1)
  BODY (solfec, 'OBSTACLE', shp, material)

  # create the remaining shapes
  glued_shape (0, 0, 2, material, solfec)

### main module ###
#import rpdb2; rpdb2.start_embedded_debugger('a')

step = 0.001

solfec = SOLFEC ('DYNAMIC', step, 'out/glue')

CONTACT_SPARSIFY (solfec, 0.005)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3)

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E15, poisson = 0.25, density = 1E1)

GRAVITY (solfec, (0, 0, -9.81))

glued_scene (bulkmat, solfec)

gs = GAUSS_SEIDEL_SOLVER (1E-3, 10000, failure = 'CONTINUE', diagsolver = 'PROJECTED_GRADIENT')

IMBALANCE_TOLERANCE (solfec, 1.1, 'ON', 2.0)

RUN (solfec, gs, 2 * step)

if not VIEWER() and solfec.mode == 'READ':

  timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL', 'PARBAL']
  dur = DURATION (solfec)
  th = HISTORY (solfec, timers, dur[0], dur[1])
  total = 0.0

  for i in range (0, len(timers)):
    sum = 0.0
    for tt in th [i+1]: sum += tt
    print timers [i], 'TIME:', sum
    if i < 6: total += sum

  print 'TOTAL TIME:', total
