# Two fuel bricks normal contact restitution test
import sys
import matplotlib.pyplot as plt
import math
import commands
sys.path.append('inp/mesh/abaqusreader')
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

# User paramters

argv = NON_SOLFEC_ARGV()

if argv == None:
  print '------------------------------------------------------'
  print 'No user paramters passed! Possible paramters:'
  print '------------------------------------------------------'
  print '-fbmod num => fuel brick modes, num >= 6 and <= 64'
  print '-damp num => damping, >= 0.0'
  print '-step num => time step, > 0.0'
  print '-afile path => Abaqus 81 array file path'
  print '-rest num => impact restitution'
  print '------------------------------------------------------'

if argv != None and len (argv) > 1:
  for i in range (0, len(argv)-1):
    if argv [i] == '-fbmod':
      fbmod = max (min (64, long (argv [i+1])), 6)
    elif argv [i] == '-damp':
      damp = max (float (argv [i+1]), 0.0)
    elif argv [i] == '-step':
      step = float (argv [i+1])
      if step <= 0.0: step = 1E-4
    elif argv [i] == '-afile':
      afile = argv [i+1]
    elif argv [i] == '-rest':
      rest = max (min (1.0, float (argv [i+1])), 0.0)

print 'Using:'
print '%d modes per fuel brick'%fbmod
print '%g damping'%damp
print '%g restitution'%rest
print '%g step'%step
print '------------------------------------------------------'

# Model

ending = '%s_%d_s%.1e_d%.1e_r%g'%(afile [afile.rfind ('/'):len(afile)].replace ('.inp',''), fbmod, step, damp, rest)

solfec = SOLFEC ('DYNAMIC', step, 'out/mbfcp/' + ending)

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.1, restitution = rest)

# Create a new AbaqusInput object from the .inp deck:
model = AbaqusInput(solfec, afile)

# Create a Finite Element body for each Instance in the Assembly:
for inst in model.assembly.instances.values():	# .instances is a dict
  label = inst.name	              # use Abaqus instance name
  mesh = inst.mesh	              # solfec MESH object at the instance position
  bulkmat = inst.material	        # solfec BULK_MATERIAL object
  bdy = BODY(solfec, 'FINITE_ELEMENT', COPY (mesh), bulkmat, label)
  data = MODAL_ANALYSIS (bdy, fbmod, solfec.outpath + '/modal' + label, abstol = 1E-13)
  DELETE (solfec, bdy)
  bdy = BODY(solfec, 'FINITE_ELEMENT', mesh, bulkmat, label, form = 'RO', modal = data)

# boundary conditions
for b in solfec.bodies:
  
  b.scheme = 'DEF_LIM'
  b.damping = damp
  
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

  b1 = None
  b2 = None
  for b in solfec.bodies:
  
      if b.label.startswith('FB'):
	if b.label.endswith('1(0)(0)'): 
	    b1 = b
	    c1 = b.center
	elif b.label.endswith('2(0)(0)'):
	    b2 = b
	    c2 = b.center
              
  p1 = TRANSLATE (c1, (0.0,-0.18985+0.005, 0.0))
  p2 = TRANSLATE (c1, (0.0, 0.18985-0.005, 0.0))
  p3 = TRANSLATE (c2, (0.0,-0.18985+0.005, 0.0))
  p4 = TRANSLATE (c2, (0.0, 0.18985-0.005, 0.0))
  th = HISTORY(solfec,[(b1,p1,'VX'), (b1,p2,'VX'),(b2,p3,'VX'),(b2,p4,'VX')],0.0,stop,1)
  plt.plot (th[0],th[1],lw = 2, label = 'Left_Brick_Top_Point_Velocity_VX_(m/s)')
  plt.plot (th[0],th[2],lw = 2, label = 'Left_Brick_Bottom_Point_Velocity_VX_(m/s)')
  plt.plot (th[0],th[3],lw = 2, label = 'Right_Brick_Top_Point_Velocity_VX_(m/s)')
  plt.plot (th[0],th[4],lw = 2, label = 'Right_Brick_Bottom_Point_Velocity_VX_(m/s)')

  res_vel_temp = zip(th[3], th[1])
  res_vel = [r1-r2 for r1,r2 in res_vel_temp]
			    
  vrest = []
  for v in res_vel:
    vrest.append(v/impactVelocity)

  plt.plot (th[0],vrest,lw = 2, label = 'Brick_Restitution')
  plt.legend(loc = 'lower left')

  plt.title ('SUB02 FB-FB Normal Contact FB damping=%s'% b.damping)
  plt.savefig (solfec.outpath + '/Brick_Velocity_Damp%g.png'% b.damping) 
  plt.show()
  
  print "For FB damping =", b.damping, " coeficient of restitution = ", vrest[-1]
    
# solver and run
GEOMETRIC_EPSILON (1E-6)
slv = NEWTON_SOLVER() #All default values used
OUTPUT(solfec, 1E-3) 
        
if RANK() == 0 and solfec.mode == 'WRITE':
  print 'Running', stop, 'seconds of analysis with step', step, '...'

if solfec.mode <> 'READ':
  RUN (solfec, slv, stop)
