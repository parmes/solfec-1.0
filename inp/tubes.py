# balls

step = 1E-3
stop = 5
outd = 0.02

solfec = SOLFEC ('DYNAMIC', step, 'out/domino')

GRAVITY (solfec, (0, 0, -10))

table = HULL ([0, 0, 0,
               0, 1, 0,
	       1, 1, 0,
	       1, 0, 0,
               0, 0, -0.1,
               0, 1, -0.1,
	       1, 1, -0.1,
	       1, 0, -0.1], 1, 1)

pipmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E5, poisson = 0.3, density = 300)

balmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E5, poisson = 0.3, density = 600)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5, restitution = 0)

SCALE (table, (10, 10, 1))
BODY (solfec, 'OBSTACLE', table, balmat)

ndir = 10
nrad = 10
nthi = 3

for i in range (0, 10):
  for j in range (0, 10):
    msh = PIPE ((0.5+i, 0.5+j, 0), (0, 0, 3), 0.2, 0.2, ndir, nrad, nthi, 2, [2, 3, 4, 5])
    pnt = []
    for k in range (0, nrad*(nthi+1)): pnt.append (msh.node (k))
    bod = BODY (solfec, 'FINITE_ELEMENT', msh, pipmat)
    bod.selfcontact = 'ON'
    for p in pnt: FIX_POINT (bod, p)

shp = SPHERE ((2, 2, 5), 2, 3, 3)
bod = BODY (solfec, 'RIGID', shp, balmat)
INITIAL_VELOCITY (bod, (1, 1, -1), (0, 0, 0))

CONTACT_EXCLUDE_SURFACES (solfec, 1, 2)

gs = GAUSS_SEIDEL_SOLVER (1E1, 1000, 1E-6)

OUTPUT (solfec, outd)

RUN (solfec, gs, stop)
