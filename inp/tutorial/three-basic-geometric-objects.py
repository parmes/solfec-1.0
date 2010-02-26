# This example illustrates using CONVEX, MESH and SPHERE objects

w = 10
l = 10
h = 1
floor_vid = 1
floor_sid = 1

floor = HULL (
       [-w/2, -l/2, -h,
         w/2, -l/2, -h,
         w/2,  l/2, -h,
        -w/2,  l/2, -h,
        -w/2, -l/2,  0,
         w/2, -l/2,  0,
         w/2,  l/2,  0,
        -w/2,  l/2,  0], floor_vid, floor_sid)

step = 1E-3

solfec = SOLFEC ('DYNAMIC', step, 'out/tutorail/three-basic-geometric-objects')

bulk = BULK_MATERIAL (solfec,
                      model = 'KIRCHHOFF',
		      young = 15E9,
		      poisson = 0.3,
		      density = 2E3)

BODY (solfec, 'OBSTACLE', floor, bulk)

a = 2
b = 2
c = 2
brick_vid = 2
brick_sid = 2

brick = CONVEX (
	  [0, 0, 0,
	   a, 0, 0,
	   a, b, 0,
	   0, b, 0,
	   0, 0, c,
	   a, 0, c,
	   a, b, c,
	   0, b, c],
	  [4, 0, 3, 2, 1, brick_sid,
	   4, 1, 2, 6, 5, brick_sid,
	   4, 2, 3, 7, 6, brick_sid,
	   4, 3, 0, 4, 7, brick_sid,
	   4, 0, 1, 5, 4, brick_sid,
	   4, 4, 5, 6, 7, brick_sid], brick_vid)

b1 = COPY (brick)
b2 = TRANSLATE (COPY (brick), (a, 0, 0))
b3 = TRANSLATE (COPY (brick), (0, b, 0))
b4 = TRANSLATE (COPY (brick), (a, b, 0))

shape = [b1, b2, b3, b4]
TRANSLATE (shape, (-a, -b, 0))
BODY (solfec, 'RIGID', shape, bulk)

nodes = [-1.0, -1.0, 0.0,
          1.0, -1.0, 0.0,
          1.0,  1.0, 0.0,
         -1.0,  1.0, 0.0,
         -1.0, -1.0, 2.0,
          1.0, -1.0, 2.0,
          1.0,  1.0, 1.0,
         -1.0,  1.0, 1.0]

mesh = HEX (nodes, 2, 3, 2, 3, [3, 3, 3, 3, 3, 3], dy = [1, 1, 2])
TRANSLATE (mesh, (0, 0, c))
BODY (solfec, 'PSEUDO_RIGID', mesh, bulk)

sphere = SPHERE ((0, 0, 5), 1, 1, 1)
body = BODY (solfec, 'RIGID', sphere, bulk)
INITIAL_VELOCITY (body, (0, 0, -10), (0, 0, 0))

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5, restitution = 0.0)

GRAVITY (solfec, (0, 0, -10))

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000)

RUN (solfec, gs, 1.0)
