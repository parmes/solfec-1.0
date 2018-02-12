###
import os, sys
def where(program):
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, program)):
      return path
  return None
path = where('parmec4')
if path == None:
  print 'ERROR: parmec4 not found in PATH!'
  print '       Download and compile parmec;'
  print '       Add parmec directory to PATH variable.'
  sys.exit(1)
sys.path.append(os.path.join (path, 'python'))
from mesh_hex import *
###

dO0 = 0.2
d01 = 0.2
d12 = 0.05
d47 = 0.01
d34 = (d12-d47)/2.
d45 = 0.02
dOz = 0.5
gap = 0.001
step = 1E-3

nodes1 = [(0, d34+gap, 0),
          (dO0+d45-gap, d34+gap, 0),
	  (dO0+d45-gap, d34+d47-gap, 0),
	  (0, d34+d47-gap, 0),
          (0, d34+gap, dOz),
          (dO0+d45-gap, d34+gap, dOz),
	  (dO0+d45-gap, d34+d47-gap, dOz),
	  (0, d34+d47-gap, dOz)]

mat0 = MATERIAL(1E3, 1E9, 0.25)

par0 = MESH_HEX(nodes1, 1, 1, 2, mat0, [0]*6)

p0 = nodes1[0]
p1 = nodes1[4]

def LINEAR_SPRING_DAMPER(part, point, leeway, ratio, accmax=1.):
  mass = EQM(part, point)
  print 'EQM:', mass
  stiff = (accmax * mass) / leeway
  damp = ratio * 2.0 * (mass * stiff)**0.5
  step = (2.0 / (stiff/mass)**0.5) * ((1.0+damp*damp)**0.5 - damp)
  print 'LSD critical step:', step
  return ([-1.0, -stiff, 1.0, stiff], [-1.0, -damp, 1.0, damp], step)

argv = ARGV()
lway = gap
if len(argv) > 1 and argv[0] == '-leeway':
  lway = float(argv[1])
print 'LEEWAY:', lway

(spr0, dsh0, h0) = LINEAR_SPRING_DAMPER (par0, p0, lway, 1.0)
(spr1, dsh1, h1) = LINEAR_SPRING_DAMPER (par0, p1, lway, 1.0)

print spr0, dsh0
print spr1, dsh1

SPRING (par0, p0, -1, p0, spr0, dsh0)
SPRING (par0, p1, -1, p1, spr1, dsh1)

step = 0.9*min(h0,h1)
print 'STEP:', step

#VELOCITY (par0, angular = (0, 0, 1))
#DEM (5.0, step, 0.1)
