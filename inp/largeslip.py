# large slip example

step = 0.001
stop = 0.5

sol = SOLFEC ('DYNAMIC', step, 'out/largeslip')

bulk = BULK_MATERIAL (sol, model = 'KIRCHHOFF', young = 200E9, poisson = 0.25, density = 8E3)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

a = 0.1
b = 0.1
c = 0.03

nodes = [-a, -b, 0,
          a, -b, 0,
          a,  b, 0,
         -a,  b, 0,
         -a, -b, c,
          a, -b, c,
          a,  b, c,
         -a,  b, c]

num = 10
msh = HEX (nodes, num, num, 2, 0, [0, 1, 2, 3, 4, 5])

shp = COPY (msh)
bod = BODY (sol, 'FINITE_ELEMENT', shp, bulk, form = 'BC')
bod.scheme = 'DEF_LIM'
bod.damping = 1E-3
for i in range (0, (num+1)**2):
  FIX_POINT (bod, msh.node (i))

shp = COPY (msh)
TRANSLATE (shp, (0, 0, c))
bod = BODY (sol, 'FINITE_ELEMENT', shp, bulk, form = 'BC')
bod.scheme = 'DEF_LIM'
bod.damping = 1E-3
PRESSURE (bod, 5, -1E3)
PRESSURE (bod, 2, -1E4)

gs = NEWTON_SOLVER (locdyn = 'OFF')

RUN (sol, gs, stop)
