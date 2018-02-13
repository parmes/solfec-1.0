import os
d0 = os.path.dirname(os.path.realpath(__file__))
execfile (d0 + '/hs4-globals.py')

argv = ARGV()
leeway = gap
quiet = False
if argv <> None and '-quiet' in argv: quiet = True
try:
  i = argv.index('-leeway')
  if len(argv) > i+1: leeway = float(argv[i+1])
except: pass

if not quiet: print 'Parmec model leeway:', leeway

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
  stiff = (accmax * mass) / leeway
  damp = ratio * 2.0 * (mass * stiff)**0.5
  hcrit = (2.0 / (stiff/mass)**0.5) * ((1.0+damp*damp)**0.5 - damp)
  return ([-1.0, -stiff, 1.0, stiff], [-1.0, -damp, 1.0, damp], hcrit)

(spr0, dsh0, h0) = LINEAR_SPRING_DAMPER (par0, p0, leeway, 1.0)
(spr1, dsh1, h1) = LINEAR_SPRING_DAMPER (par0, p1, leeway, 1.0)

if not quiet:
  print 'Paremc model spring and damper curves:'
  print spr0, dsh0
  print spr1, dsh1

SPRING (par0, p0, -1, p0, spr0, dsh0)
SPRING (par0, p1, -1, p1, spr1, dsh1)

step = 0.9*min(h0,h1)
if not quiet: print 'Parmec model suggested step:', step

#VELOCITY (par0, angular = (0, 0, 1))
#DEM (5.0, step, 0.1)
