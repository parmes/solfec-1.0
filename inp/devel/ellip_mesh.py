sol = SOLFEC ('DYNAMIC', 1E-3, 'out/ellip_mesh')

mat = BULK_MATERIAL (sol, model = 'KIRCHHOFF', young = 1E9, poisson = 0.25, density = 1E3)

msh = ELLIP_MESH ((0, 0, 0), (0.6, 0.8, 1), 0.01, 1, 1)

bod = BODY (sol, 'FINITE_ELEMENT', msh, mat)

OUTPUT (sol, 1E-2)

RUN (sol, NEWTON_SOLVER(), 1.0)
