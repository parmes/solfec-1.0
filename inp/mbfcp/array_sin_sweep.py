# array of bricks
import sys
sys.path.append ('inp/mesh')
from abaqusread import *

#import rpdb2; rpdb2.start_embedded_debugger ('a')

step = 1E-3

solfec = SOLFEC ('DYNAMIC', step, 'out/mbfcp/array_sin_sweep')

ABAQUS_READ ('/Users/tomek/Desktop/A1.inp', solfec)

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.1, restitution = 0.0)
