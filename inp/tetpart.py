# tetrahedral mesh example
import sys
sys.path.append ('inp/mesh')
from netgenread import *

step = 1E-3

stop = 0.5

solfec = SOLFEC ('DYNAMIC', step, 'out/tetpart')

solver = GAUSS_SEIDEL_SOLVER (1E-5, 100)

material = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5, restitution = 0.25)

shape = NETGEN_READ ('inp/mesh/part0.mesh', 1, 1)

BODY (solfec, 'RIGID', shape, material)

shape = HULL ([0, 0, 0,
               1, 0, 0,
	       1, 1, 0,
	       0, 1, 0,
	       0, 0, 1,
	       1, 0, 1,
	       1, 1, 1,
	       0, 1, 1], 1, 1)

oshp = TRANSLATE (SCALE (COPY (shape), (6, 6, 0.5)), (-2, -2, -0.5))
BODY (solfec, 'OBSTACLE', oshp, material)

GRAVITY (solfec, (0, 0, -10))

RUN (solfec, solver, stop)
