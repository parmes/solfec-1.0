a = 1
b = 1
step = 0.001
stop = 1.0
sol = SOLFEC ('DYNAMIC', step, 'out/partition')
c = NCPU (sol)

nodes = [-a, -b, 0,
          a, -b, 0,
          a,  b, 0,
         -a,  b, 0,
         -a, -b, c,
          a, -b, c,
          a,  b, c,
         -a,  b, c]

msh = HEX (nodes, 1, 1, NCPU (sol), 0, [0, 1, 2, 3, 4, 5])

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 100,
		      poisson = 0,
		      density = 1)

bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)
FIX_POINT (bod, (0, 0, c))
PARTITION (bod, NCPU (sol))

GRAVITY (sol, (0, 0, -1))

gs = GAUSS_SEIDEL_SOLVER (1E-5, 1000)

RUN (sol, gs, stop)
