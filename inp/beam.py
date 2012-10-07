# soft free supported beam bending

b = 0.2
h = 0.4
l = 4.0
step = 0.01 # large step (soft)
stop = 0.5

nodes = [0, 0, 0,
         1, 0, 0,
         1, 1, 0,
         0, 1, 0,
         0, 0, 1,
         1, 0, 1,
         1, 1, 1,
         0, 1, 1]

beam_mesh = HEX (nodes, 40, 4, 8, 1, [1, 1, 1, 1, 1, 3])
SCALE (beam_mesh, (l, b, h))

supp_mesh = HEX (nodes, 1, 1, 1, 2, [2, 2, 2, 2, 2, 2])
SCALE (supp_mesh, (b, b, b))

solfec = SOLFEC ('QUASI_STATIC', step, 'out/beam')

bulk = BULK_MATERIAL (solfec,
                      model = 'KIRCHHOFF',
		      young = 30E7, # softer to produce visible deformation
		      poisson = 0.25,
		      density = 2E3)

body = BODY (solfec, 'FINITE_ELEMENT', beam_mesh, bulk)
PRESSURE (body, 3, -100E3) # large pressure

supp = TRANSLATE (COPY (supp_mesh), (0, 0, -b))
BODY (solfec, 'OBSTACLE', supp, bulk)
supp = TRANSLATE (COPY (supp_mesh), (l-b, 0, -b))
BODY (solfec, 'OBSTACLE', supp, bulk)

solver = GAUSS_SEIDEL_SOLVER (1E-4, 100)

OUTPUT (solfec, step)

RUN (solfec, solver, stop)
