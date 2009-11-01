# test of Solfec functionality with geometry defined by convex polyhedra

w = 1.0
l = 1.0
h = 1.0
a = w/3.0
b = l/3.0
c = a + b
n = 4
m = 4
step = 1E-4

base = HULL (
       [-w/2, -l/2, -h,
         w/2, -l/2, -h,
         w/2,  l/2, -h,
        -w/2,  l/2, -h,
        -w/2, -l/2,  0,
         w/2, -l/2,  0,
         w/2,  l/2,  0,
        -w/2,  l/2,  0], 1, 1)

nodes = [-a, -b, 0,
          a, -b, 0,
          a,  b, 0,
         -a,  b, 0,
         -a, -b, c,
          a, -b, c,
          a,  b, c,
         -a,  b, c]

elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, 0]

#msh1 = MESH (nodes, elements, 0)

msh1 = HEX (nodes, 2, 1, 1, 0, [0, 1, 2, 3, 4, 5])

sol = SOLFEC ('DYNAMIC', step, 'out/fem0')

sur = SURFACE_MATERIAL (sol,
                        model = 'SIGNORINI_COULOMB',
                        friction = 0,
			spring = 1E3,
			dashpot = 0)

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E6,
		      poisson = 0.0,
		      density = 1E1)

#BODY (sol, 'OBSTACLE', base, bulk)
bod = BODY (sol, 'FINITE_ELEMENT', msh1, bulk)
#bod = BODY (sol, 'PSEUDO_RIGID', msh1, bulk)
FIX_POINT (sol, bod, (-a, -b, 0))
FIX_POINT (sol, bod, (-a, b, 0))

gs = GAUSS_SEIDEL_SOLVER (1E-5, 1000)

GRAVITY (sol, (0, 0, -1), 10)

OUTPUT (sol, 1 * step)

RUN (sol, gs, 10 * step)
