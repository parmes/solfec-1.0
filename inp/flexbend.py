# flexing bent mesh example (reduced order)
step = 0.002
stop = 10

sol = SOLFEC ('DYNAMIC', step, 'out/flexbend')
bulk = BULK_MATERIAL (sol, model = 'KIRCHHOFF', young = 1E8, poisson = 0.3, density = 1E3)
SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

nodes = [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1]
msh = HEX (nodes, 15, 10, 2, 0, [0, 0, 0, 0, 0, 0])
SCALE (msh, (15, 10, 1))
BEND (msh, (0, 0, -3), (-1, 0, 0), 270)
bod = BODY (sol, 'FINITE_ELEMENT', COPY (msh), bulk)
data = MODAL_ANALYSIS (bod, 12, 'out/flexbend/modal')
DELETE (sol, bod)
BODY (sol, 'FINITE_ELEMENT', msh, bulk, form = 'RO', modal = data)

shp = HEX (nodes, 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])
SCALE (shp, (15, 15, 1))
TRANSLATE (shp, (0, -2, -10))
BODY (sol, 'OBSTACLE', shp, bulk)

sv = NEWTON_SOLVER ()
GRAVITY (sol, (0, 0, -10))
OUTPUT (sol, 0.01)
RUN (sol, sv, stop)
