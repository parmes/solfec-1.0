# array of bricks
import sys
sys.path.append ('inp/mesh')
from abaqusread import *

#import rpdb2; rpdb2.start_embedded_debugger('a')

solfec = ABAQUS_READ ('/Users/tomek/Desktop/A1.inp', 1E-3, 'out/mbfcp/array_sin_sweep')
