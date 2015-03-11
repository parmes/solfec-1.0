# Two fuel bricks normal contact restitution test
import sys
import math
import commands
sys.path.append('scripts/abaqusreader')
sys.path.append('inp/mbfcp')
from abaqusreader import AbaqusInput
from math import cos 

# Analysis paramters
step = 1E-4
stop = 0.5
damp = 1E-6
impactVelocity = 0.15 #15cm/s
fbmod = 12
afile = 'inp/mesh/81fbi.inp'
rest = 0.0
refine = 0

# Model
solfec = SOLFEC ('DYNAMIC', step, 'out/fratest')

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.1, restitution = rest)

# Create a new AbaqusInput object from the .inp deck:
model = AbaqusInput(afile, solfec)

# Create bulk material with critical energy set up
bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3, tensile = 5E5)

# Create a Finite Element body for each Instance in the Assembly:
for inst in model.assembly.instances.values():	# .instances is a dict
  label = inst.name	              # use Abaqus instance name
  mesh = inst.mesh	              # solfec MESH object at the instance position

  if refine:
    meshin = TETRAHEDRALIZE (mesh, solfec.outpath + '/' + label + '.meshdata', 0.2, 1.5, min_angle = 0, max_angle = 100, ref_length = 0.01)
  else:
    meshin = COPY (mesh)
  
  bdy = BODY(solfec, 'FINITE_ELEMENT', COPY (meshin), bulkmat, label)
  data = MODAL_ANALYSIS (bdy, fbmod, solfec.outpath + '/modal' + label, abstol = 1E-13)
  DELETE (solfec, bdy)
  bdy = BODY(solfec, 'FINITE_ELEMENT', meshin, bulkmat, label, form = 'RO', modal = data)

GEOMETRIC_EPSILON (1E-2)
# boundary conditions
for b in solfec.bodies:
  
  b.scheme = 'DEF_LIM'
  b.damping = damp
  b.fracturecheck = 'ON' # enable fracture check
  
  c = b.center
  t = b.tensor
  
  print "body:" , b.mass
  if b.label.startswith('FB'): # find FB
    # vertical constraints
    p1 = TRANSLATE (c, (-0.18985+0.005,0.0, 0.4535-0.005))
    FIX_DIRECTION (b, p1, (0.0, 0.0, 1.0))
    FIX_DIRECTION (b, p1, (0.0, 1.0, 0.0))
    p2 = TRANSLATE (c, (0.18985-0.005, 0.0, 0.4535-0.005))
    FIX_DIRECTION (b, p2, (0.0, 0.0, 1.0))
    FIX_DIRECTION (b, p2, (0.0, 1.0, 0.0))
    p3 = TRANSLATE (c, (-0.18985+0.005, 0.0, -0.4535+0.005))
    FIX_DIRECTION (b, p3, (0.0, 1.0, 0.0))
    p4 = TRANSLATE (c, (0.18985-0.005, 0.0, -0.4535+0.005))
    FIX_DIRECTION (b, p4, (0.0, 1.0, 0.0))
    
    #ASSIGN INPUT VELOCITY
    if b.label.startswith('FB1(0)(0)'):
        INITIAL_VELOCITY (b, (impactVelocity/2.0, 0.0, 0.0),(0.0, 0.0, 0.0))
    elif b.label.startswith('FB2(0)(0)'):
        INITIAL_VELOCITY (b, (0.0-(impactVelocity/2.0), 0.0, 0.0),(0.0, 0.0, 0.0))
        
if not VIEWER() and solfec.mode == 'READ':

  print 'Exportind data to YAFFEMS...'
  for b in solfec.bodies:
    # export data for YAFFEMS fracture analysis
    n = FRACTURE_EXPORT_YAFFEMS (b, 'out/fratest/fra_' + b.label, 0.5)
    print 'Body ', b.label, 'has ', n, 'instances of fracture analysis'
    
# solver and run
GEOMETRIC_EPSILON (1E-6)
slv = NEWTON_SOLVER() #All default values used
OUTPUT(solfec, 1E-3) 
        
if RANK() == 0 and solfec.mode == 'WRITE':
  print 'Running', stop, 'seconds of analysis with step', step, '...'

if solfec.mode <> 'READ':
  RUN (solfec, slv, stop)
