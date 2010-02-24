# pinned fem box example 

a = 1.0
b = 1.0
c = 2.0
step = 0.01 # let the critical step rule

nodes = [-a, -b, 0,
          a, -b, 0,
          a,  b, 0,
         -a,  b, 0,
         -a, -b, c,
          a, -b, c,
          a,  b, c,
         -a,  b, c]

msh = HEX (nodes, 1, 1, 1, 0, [0, 1, 2, 3, 4, 5])

sol = SOLFEC ('DYNAMIC', step, 'out/pinned-fem-box')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 15E9,
		      poisson = 0.25,
		      density = 1.8E3)

bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)
bod.scheme = 'DEF_IMP'
FIX_POINT (sol, bod, (-a, -b, 0))
FIX_POINT (sol, bod, (-a, b, 0))


gs = GAUSS_SEIDEL_SOLVER (1E-5, 1000)

GRAVITY (sol, (0, 0, -1), 10)

OUTPUT (sol, 0.01)

RUN (sol, gs, 4.0)
