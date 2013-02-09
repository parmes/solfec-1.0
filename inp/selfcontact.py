# slef-contact example
# FIXME: this example illustrates an issue with
#        detection of near contact points with
#        constradictory normals; this leads to
#        unbounded growth of reactions in GS solver

step = 0.005
stop = 1

sol = SOLFEC ('DYNAMIC', step, 'out/selfcontact')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E8,
		      poisson = 0.3,
		      density = 1E3)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

nodes = [0, 0, 0,
         1, 0, 0,
         1, 1, 0,
         0, 1, 0,
         0, 0, 1,
         1, 0, 1,
         1, 1, 1,
         0, 1, 1]

msh = HEX (nodes, 2, 30, 1, 0, [0, 0, 0, 0, 0, 0])
SCALE (msh, (1, 30, 1))
BEND (msh, (0, 10, 4), (1, 0, 0), 219)
bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)
bod.selfcontact = 'ON'

FIX_POINT (bod, (0, 0, 0))
FIX_POINT (bod, (1, 0, 0))

#shp = HEX (nodes, 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])
#SCALE (shp, (3, 15, 1))
#TRANSLATE (shp, (-1, -2, -1))
#bod = BODY (sol, 'OBSTACLE', shp, bulk)

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000)

GRAVITY (sol, (0, 0, -10))

OUTPUT (sol, step)

RUN (sol, gs, stop)
