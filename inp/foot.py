# tetrahedral mesh example

GEOMETRIC_EPSILON (1E-8)
step = 0.001
stop = 5
sv = GAUSS_SEIDEL_SOLVER (1E-4, 1000)

sol = SOLFEC ('DYNAMIC', step, 'out/foot')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E6,
		      poisson = 0.3,
		      density = 1E3)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

GRAVITY (sol, (0, 0, -10))

foot = TETRAHEDRALIZE ('inp/mesh/foot.poly', 'out/foot/mesh.dat', 0.2, 1.1)
TRANSLATE (foot, (0, 0, 1))
bod = BODY (sol, 'FINITE_ELEMENT', foot, bulk)
bod.damping = 0.001

base = HEX ([-1, -1, -1,
	     11, -1, -1,
	     11,  5, -1,
	     -1,  5, -1,
             -1, -1, 0,
	     11, -1, 0,
	     11,  5, 0,
	     -1,  5, 0], 1, 1, 1, 1, [2, 2, 2, 2, 2, 2])
bod = BODY (sol, 'OBSTACLE', base, bulk)

RUN (sol, sv, stop)
