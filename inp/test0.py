# test of Solfec functionality with geometry defined by convex polyhedra
import math

w = 10
l = 10
h = 1
a = 0.2
b = 0.1
c = 0.1
n = 4
m = 4
step = 1E-4
time = 1.0

base = HULL (
       [-w/2, -l/2, -h,
         w/2, -l/2, -h,
         w/2,  l/2, -h,
        -w/2,  l/2, -h,
        -w/2, -l/2,  0,
         w/2, -l/2,  0,
         w/2,  l/2,  0,
        -w/2,  l/2,  0], 1, 1)

def makebrick (sur, vol):
  brick0 = CONVEX (
	  [0, 0, 0,
	   a/2, 0, 0,
	   a/2, b/2, 0,
	   0, b/2, 0,
	   0, 0, c/2,
	   a/2, 0, c/2,
	   a/2, b/2, c/2,
	   0, b/2, c/2],
	  [4, 0, 3, 2, 1, sur,
	   4, 1, 2, 6, 5, sur,
	   4, 2, 3, 7, 6, sur,
	   4, 3, 0, 4, 7, sur,
	   4, 0, 1, 5, 4, sur,
	   4, 4, 5, 6, 7, sur], vol)

  brick1 = TRANSLATE (COPY (brick0), (a/2, 0, 0))
  brick2 = TRANSLATE (COPY (brick0), (a/2, b/2, 0))
  brick3 = TRANSLATE (COPY (brick0), (0, b/2, 0))
  brick4 = TRANSLATE (COPY (brick0), (0, 0, c/2))
  brick5 = TRANSLATE (COPY (brick4), (a/2, 0, 0))
  brick6 = TRANSLATE (COPY (brick4), (a/2, b/2, 0))
  brick7 = TRANSLATE (COPY (brick4), (0, b/2, 0))

  return [brick0, brick1, brick2, brick3, brick4, brick5, brick6, brick7]

nodes = [-0.05, -0.05, 0.0,
          0.05, -0.05, 0.0,
          0.05,  0.05, 0.0,
         -0.05,  0.05, 0.0,
         -0.05, -0.05, 1.0,
          0.05, -0.05, 1.0,
          0.05,  0.05, 1.0,
         -0.05,  0.05, 1.0]


elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, 0]

msh1 = MESH (nodes, elements, 0)

msh2 = HEX (nodes, 2, 3, 2, 0, [0, 1, 2, 3, 4, 5], dy = [1, 1, 2])

sph = SPHERE ((0, 0, 0), 1, 1, 1)

sph = SPHERE ((0, 0, 1), 1, 1, 1, sph)

sol = SOLFEC ('DYNAMIC', step, 'out/test0')

#sol2 = SOLFEC ('DYNAMIC', step, 'cvxwall.out')
#sol1 = sol #seems like reference counting results in calling destructor of sol1 here
#sol3 = SOLFEC ('DYNAMIC', step, 'cvxwall.out')

sur = SURFACE_MATERIAL (sol,
                        model = 'SIGNORINI_COULOMB',
                        friction = 0.3,
			spring = 1E3,
			dashpot = 0)

#sur.cohesion = 0.5 + 10 + math.sin (0.123)
#print 'cohesion = ', sur.cohesion

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E9,
		      poisson = 0.25,
		      density = 1E3)

#print bulk.label, 'young = ', bulk.young

base = BODY (sol, 'OBSTACLE', base, bulk)

hex = TRANSLATE (SCALE (msh2, (10, 10, 1)), (0, 3, 0))
hex = ROTATE (hex, (0, 3, 0), (0, 0, 1), 45)

BODY (sol, 'OBSTACLE', hex, bulk)

bod = BODY (sol, 'OBSTACLE', SPHERE ((0, -3, 1), 1, 4, 5), bulk)

#FORCE (bod, 'SPATIAL', MASS_CENTER (bod), (0, 1, 0), 1e6)

tms = TIME_SERIES ([0, 1, 3, 4])

def gscallback (gs):
  print gs.error
  return 0

gs = GAUSS_SEIDEL_SOLVER (1E-2, 1000, failure = 'CALLBACK', callback=gscallback)
gs.history = 'ON'

ex = EXPLICIT_SOLVER ()

GRAVITY (sol, (0, 0, -1), 10)

for i in range (-n/2, n/2):
  for j in range (m):

    point = (i*a, -b/2, j*c)
    brick = makebrick (j + i, i - j)
    copy = TRANSLATE (COPY (brick), point)
    point = (i*a+a/2, -b/2, j*c)
    shp = None

    if i == -n/2:
      (one, shp) = SPLIT (copy, point, (1, 0, 0))
    elif i == n/2-1:
      (shp, one) = SPLIT (copy, point, (1, 0, 0))
    else: shp = copy

    BODY (sol, 'RIGID', shp, bulk)

OUTPUT (sol, 10 * step)

EXTENTS (sol, (-w, -l, -2 * h, w, l, 2 * h))

def runcallback (gs):
  if gs.error != 0: print gs.error
  else: print 'gs iters = ', gs.iters

CALLBACK (sol, 0.01, gs, runcallback)

RUN (sol, ex, time)
