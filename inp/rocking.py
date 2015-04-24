# pinned fem box example 

a = 0.1
b = 0.1
c = 0.1
step = 0.001 # let the critical step rule
stop = 10

nodes = [-a, -b, 0,
          a, -b, 0,
          a,  b, 0,
         -a,  b, 0,
         -a, -b, c,
          a, -b, c,
          a,  b, c,
         -a,  b, c]

msh = HEX (nodes, 3, 3, 3, 0, [0, 1, 2, 3, 4, 5])

sol = SOLFEC ('DYNAMIC', step, 'out/rocking')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E9,
		      poisson = 0.25,
		      density = 1E3)

ROTATE (msh, (0, 0, 0), (1, 0, 0), 45)
bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)
bod.scheme = 'DEF_LIM'
bod.damping = 1E-5


msh = HEX (nodes, 1, 1, 1, 0, [0, 1, 2, 3, 4, 5])
SCALE (msh, (10, 10, 1))
TRANSLATE (msh, (0, 0, -0.5))
bod = BODY (sol, 'OBSTACLE', msh, bulk)

gs = GAUSS_SEIDEL_SOLVER (1E-5, 1000)

GRAVITY (sol, (0, 0, -10))

OUTPUT (sol, step)

RUN (sol, gs, stop)
