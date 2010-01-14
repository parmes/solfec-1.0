# pinned fem box example 

a = 1.0
b = 1.0
c = 2.0
step = 1E-4

nodes = [-a, -b, 0,
          a, -b, 0,
          a,  b, 0,
         -a,  b, 0,
         -a, -b, c,
          a, -b, c,
          a,  b, c,
         -a,  b, c]

msh = HEX (nodes, 2, 2, 2, 0, [0, 1, 2, 3, 4, 5])

sol = SOLFEC ('DYNAMIC', step, 'out/pinned-fem-box')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E6,
		      poisson = 0.0,
		      density = 1E2)

bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)
FIX_POINT (sol, bod, (-a, -b, 0))
FIX_POINT (sol, bod, (-a, b, 0))

gs = GAUSS_SEIDEL_SOLVER (1E-5, 1000)

GRAVITY (sol, (0, 0, -1), 10)

OUTPUT (sol, 10 * step)

RUN (sol, gs, 4.0)
