step = 0.01

sol = SOLFEC ('DYNAMIC', step, 'out/hybrid-solver0')

GRAVITY (sol, (0, 0, -10))

mat = BULK_MATERIAL (sol, model = 'KIRCHHOFF',
       young = 1, poisson = 0.25, density = 1)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.1)

nodes = [0, 0, 0,
         1, 0, 0,
	 1, 1, 0,
	 0, 1, 0,
	 0, 0, 1,
	 1, 0, 1,
	 1, 1, 1,
	 0, 1, 1]

msh = HEX (nodes, 1, 1, 1, 0, [0, 1, 2, 3, 4, 5])
bod1 = BODY (sol, 'RIGID', msh, mat) # boundary bodies are rigid

msh = HEX (nodes, 2, 2, 2, 0, [0, 1, 2, 3, 4, 5])
TRANSLATE (msh, (0, 0, 1.1))
bod2 = BODY (sol, 'RIGID', msh, mat)

ns = NEWTON_SOLVER ()

# nubering of bodies in Parmec starts from 0 hence below we
# use dictionary {0 : bod1.id} as the parmec2solfec mapping
hs = HYBRID_SOLVER ('examples/hybrid-solver0/hs0-parmec.py',
                    step, {0 : bod1.id}, ns)

import solfec as solfec # this is needed since 'OUTPUT' in Solfec
solfec.OUTPUT (sol, 0.03) # collides with 'OUTPUT' in Parmec

RUN (sol, hs, 10)
