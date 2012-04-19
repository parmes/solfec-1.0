# tetrahedral mesh example

#GEOMETRIC_EPSILON (1E-4)
step = 0.0005
stop = 5
sv = GAUSS_SEIDEL_SOLVER (1E-4, 1000)

sol = SOLFEC ('DYNAMIC', step, 'out/tetgen')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 15E6,
		      poisson = 0.3,
		      density = 1E3)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

GRAVITY (sol, (0, 0, -10))

shp = PIPE  ((0, 0, 0), (0, 0, 2), 0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])

(b, a) = SPLIT (COPY (shp), (0, 0, 1), (1, 1, 1))
b = TETRAHEDRALIZE (b, 'out/tetgen/tet1.dat', 0.2, 1.5)
a = TETRAHEDRALIZE (a, 'out/tetgen/tet2.dat', 0.2, 1.5)

bod = BODY (sol, 'FINITE_ELEMENT', a, bulk)
bod = BODY (sol, 'FINITE_ELEMENT', b, bulk)

shp = PIPE  ((0, 0, -1), (0, 0, -1), 0.5, 1, 2, 8, 1, 1, [1, 1, 1, 1, 1, 1])
bod = BODY (sol, 'OBSTACLE', shp, bulk)

OUTPUT (sol, step)
RUN (sol, sv, stop)
