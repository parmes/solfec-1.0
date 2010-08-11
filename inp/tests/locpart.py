# locdyn from partitioned domain test

step = 0.002
stop = step

sol = SOLFEC ('DYNAMIC', step, 'out/tests/locpart')

if not VIEWER(): sol.verbose = 'OFF'

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E7,
		      poisson = 0.3,
		      density = 1E3)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

nodes = [0, 0, 0,
         1, 0, 0,
         1, 1, 0,
         0, 1, 0,
         0, 0, 1,
         1, 0, 1,
         1, 1, 1,
         0, 1, 1]

msh = HEX (nodes, 15, 15, 4, 0, [0, 0, 0, 0, 0, 0])
SCALE (msh, (15, 10, 1))
BEND (msh, (0, 0, -3), (-1, 0, 0), 270)
BEND (msh, (5, 7, 0), (0, 0, 1), 90)
bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)
PARTITION (bod, 8)

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000)

GRAVITY (sol, (0, 0, -10))

RUN (sol, gs, stop)

LOCDYN_DUMP (sol, 'out/tests/locpart/locpart-' + str (NCPU(sol)))
