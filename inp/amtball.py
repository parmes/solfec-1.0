# Advanced Manufacturing Techniques contact example 2
step = 0.001
stop = 0.3

sol = SOLFEC ('DYNAMIC', step, 'out/amtball')
bulk = BULK_MATERIAL (sol, model = 'KIRCHHOFF', young = 1E5, poisson = 0.2, density = 1000)
SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.2)

nodes = [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1]
msh = HEX (nodes, 1, 1, 1, 0, [0, 0, 0, 0, 0, 0])
SCALE (msh, (0.02, 0.02, 0.01))
bod = BODY (sol, 'OBSTACLE', msh, bulk)

msh = HEX (nodes, 15, 15, 7, 0, [0, 0, 0, 0, 0, 0])
SCALE (msh, (0.02, 0.02, 0.01))
TRANSLATE (msh, (0, 0, 0.01))
bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)

sph = SPHERE ((0.01, 0.01, 0.025), 0.005, 1, 1)
bod = BODY (sol, 'RIGID', sph, bulk)
FIX_DIRECTION (bod, bod.center, (1, 0, 0))
SET_VELOCITY (bod, bod.center, (0, 1, 0), 0.005)
SET_VELOCITY (bod, bod.center, (0, 0, 1), -0.005)
TORQUE (bod, 'SPATIAL', (0, 0, 1), 1)

gs = GAUSS_SEIDEL_SOLVER (1E-4, 1000)

RUN (sol, gs, stop)
