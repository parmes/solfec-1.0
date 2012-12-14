# spaghetti breaking example

b = 0.002
h = 0.002
l = 0.3
nl = 200
step = 0.001
stop = 5

nodes = [0, 0, 0,
         1, 0, 0,
         1, 1, 0,
         0, 1, 0,
         0, 0, 1,
         1, 0, 1,
         1, 1, 1,
         0, 1, 1]

pasta_mesh = HEX (nodes, nl, 1, 1, 1, [1]*6)
SCALE (pasta_mesh, (l, b, h))

supp_mesh = HEX (nodes, 1, 1, 1, 2, [2, 2, 2, 2, 2, 2])
SCALE (supp_mesh, (b, b, b))

solfec = SOLFEC ('DYNAMIC', step, 'out/spaghetti')

bulk = BULK_MATERIAL (solfec,
                      model = 'KIRCHHOFF',
		      young = 1E6, # softer to produce visible deformation
		      poisson = 0.25,
		      density = 500)

body = BODY (solfec, 'FINITE_ELEMENT', pasta_mesh, bulk)
for i in range (1, nl):
  x = i*(l/nl)
  SIMPLIFIED_CRACK (body, (x, b/2., h/2.), (1, 0, 0), 1, 'TENSILE', ft=1E3)
body.scheme = 'DEF_LIM'
body.damping = 0.1*step

supp = TRANSLATE (COPY (supp_mesh), (b, 0, -b))
BODY (solfec, 'OBSTACLE', supp, bulk)
supp = TRANSLATE (COPY (supp_mesh), (l-2*b, 0, -b))
BODY (solfec, 'OBSTACLE', supp, bulk)
supp = TRANSLATE (COPY (supp_mesh), (l/2-b/2, 0, b))
body = BODY (solfec, 'OBSTACLE', supp, bulk)
SET_VELOCITY (body, body.center, (0, 0, 1), -0.005)

solver = GAUSS_SEIDEL_SOLVER (1E-4, 100)
OUTPUT (solfec, step)
RUN (solfec, solver, stop)
