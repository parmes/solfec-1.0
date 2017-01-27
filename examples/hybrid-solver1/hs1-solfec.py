# TODO --> rewrite below

from shutil import copyfile

step = 0.01

sol = SOLFEC ('DYNAMIC', step, 'out/hybrid-solver0')

GRAVITY (sol, (0, 0, -10))

mat = BULK_MATERIAL (sol, model = 'KIRCHHOFF', young = 1, poisson = 0.25, density = 1)

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
bod1 = BODY (sol, 'OBSTACLE', msh, mat) # boundary bodies need to be obstacles

msh = HEX (nodes, 2, 2, 2, 0, [0, 1, 2, 3, 4, 5])
TRANSLATE (msh, (0, 0, 1.1))
bod2 = BODY (sol, 'RIGID', msh, mat)

ns = NEWTON_SOLVER ()

# parmec's output files are written to the same location as the input path
# for that to be the solfec's output directory, we copy parmec's input file there
copyfile('examples/hybrid-solver0/hs0-parmec.py', 'out/hybrid-solver0/hs0-parmec.py')

# nubering of bodies in Parmec starts from 0 while in Solfec from 1
# hence below we used dictionary {0 : 1} as the parmec2solfec mapping
hs = HYBRID_SOLVER ('out/hybrid-solver0/hs0-parmec.py', step, {0 : 1}, ns, step)

import solfec as solfec # we need to be specific when using the OUTPUT command
solfec.OUTPUT (sol, 0.03) # since 'OUTPUT' in Solfec collides with 'OUTPUT' in Parmec

RUN (sol, hs, 10)
