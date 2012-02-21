# Model name 81array.py
# Abaqus input file = 81_Brick_Model_14.inp - Revised pitch and gaps as per test report after verification completed 26/09/11
# 81_Brick_Model_14.inp also includes modified FB and IB Young's modulus FB = 0.138GPa, IB = 0.138GPa, LK = 11.8GPa.
# Array of bricks with brick normal gaps to match those in test, DYNAMIC run all PSEUDO-RIGID bodies - Loose Key facets removed
# VELOCITY Time-history used from 5101060/23/03/08 - 2s, 0.3g, 3Hz sine dwell, then 3Hz to 10Hz linear sweep at 0.1Hz/s

import sys
import math
import sys
import time
import commands
sys.path.append('inp/mesh/abaqusreader')
sys.path.append('inp/mbfcp')
from abaqusreader import AbaqusInput
from math import cos 

# Analysis inputs

input_bricks = (['FB1(0)(0)', 'FB1(0)(1)', 'FB1(0)(2)', 'FB1(0)(3)', 'FB1(0)(4)', 'FB1(0)(5)', #left side
                'FB1(1)(0)', 'FB1(2)(0)', 'FB1(3)(0)', 'FB1(4)(0)', #top
                'FB1(1)(5)', 'FB1(2)(5)', 'FB1(3)(5)', 'FB1(4)(5)', #bottom
                'FB1(5)(0)', 'FB1(5)(1)', 'FB1(5)(2)', 'FB1(5)(3)', 'FB1(5)(4)', 'FB1(5)(5)', #right
                'IB2(0)(0)', 'IB2(0)(1)', 'IB2(0)(2)', 'IB2(0)(3)', 'IB2(0)(4)', #Side 1 (left)
                'IB1(0)(5)', 'IB1(1)(5)', 'IB1(2)(5)', 'IB1(3)(5)', 'IB1(4)(5)', #Side 2 (top)
                'IB2(5)(0)', 'IB2(5)(1)', 'IB2(5)(2)', 'IB2(5)(3)', 'IB2(5)(4)', #Side 3 (right)
                'IB1(0)(0)', 'IB1(1)(0)', 'IB1(2)(0)', 'IB1(3)(0)', 'IB1(4)(0)']) #Side 4 (bottom)

step = 1E-4
stop = 72.0

solfec = SOLFEC ('DYNAMIC', step, 'out/mbfcp/81array')
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.1, restitution = 0.0)

if RANK () == 0:
  commands.getoutput ("bunzip2 inp/mbfcp/ts81.py.bz2") # only first CPU unpacks the input

BARRIER () # let all CPUs meet here

from ts81 import TS81 # import the time series

#This is a velocity-time series, for a 0.3g acceleration sine input at 3Hz, 2s dwell period, the 3Hz to 10Hz dwell sweep at 0.1Hz/s
vel = TS81()

if RANK () == 0:
  commands.getoutput ("bzip2 inp/mbfcp/ts81.py") # only first CPU pack the input

# Create a new AbaqusInput object from the .inp deck:
model = AbaqusInput(solfec, 'inp/mesh/81array.inp')

# Create a Finite Element body for each Instance in the Assembly:
for inst in model.assembly.instances.values():	# .instances is a dict
  label = inst.name	              # use Abaqus instance name
  mesh = inst.mesh	              # solfec MESH object at the instance position
  bulkmat = inst.material	        # solfec BULK_MATERIAL object
  bdy = BODY(solfec, 'FINITE_ELEMENT', COPY (mesh), bulkmat, label)
  data = MODAL_ANALYSIS (bdy, 12, 'out/mbfcp/81array/modal' + label)
  DELETE (solfec, bdy)
  bdy = BODY(solfec, 'FINITE_ELEMENT', mesh, bulkmat, label, form = 'RO', modal = data)


#----------------------------------------------------------------------

# boundary conditions and input accelerations

for b in solfec.bodies:

  #print "body:", b.label, b.center


  b.scheme = 'DEF_LIM'
  c = b.center

  print "body mass:", b.label, b.mass, "Kg"
  
  if b.label.startswith('FB'): # find FBs
  
    # Set fuel brick damping
  
    b.damping = 1E-5

    # vertical constraints

    p1 = TRANSLATE (c, (-0.1, 0.1, 0.0))
    FIX_DIRECTION (b, p1, (0.0, 0.0, -1.0))
    p2 = TRANSLATE (c, (-0.1, -0.1, 0.0))
    FIX_DIRECTION (b, p2, (0.0, 0.0, -1.0))
    p3 = TRANSLATE (c, (0.14, 0.0, 0.0))
    FIX_DIRECTION (b, p3, (0.0, 0.0, -1.0))

    #FBs on boundary assigned velocity inputs and fixed in y direction

    if b.label in input_bricks:

        SET_VELOCITY (b, p3, (1.0, 0.0, 0.0), vel) # added motion
        FIX_DIRECTION (b, p2, (0.0, 1.0, 0.0)) # Fix point p2 in y direction
        FIX_DIRECTION (b, p3, (0.0, 1.0, 0.0)) # Fix point p3 in y direction

  elif b.label.startswith('IB'):
  
    # Set interstitial brick damping
    
    b.damping = 1E-4

    # vertical constraints

    p4 = TRANSLATE (c, (-0.07, 0.07, 0.0))
    FIX_DIRECTION (b, p4, (0.0, 0.0, -1.0))
    p5 = TRANSLATE (c, (-0.07, -0.07, 0.0))
    FIX_DIRECTION (b, p5, (0.0, 0.0, -1.0))
    p6 = TRANSLATE (c, (0.098, 0.0, 0.0))
    FIX_DIRECTION (b, p6, (0.0, 0.0, -1.0))
    
    #IBs on boundary assigned velocity inputs and fixed in y direction

    if b.label in input_bricks:

        SET_VELOCITY (b, p6, (1, 0, 0), vel) # added motion
        FIX_DIRECTION (b, p5, (0, 1, 0)) # Fix point p5 in y direction
        FIX_DIRECTION (b, p6, (0, 1, 0)) # Fix point p6 in y direction

  elif b.label.startswith('LK'):
  
    # Set loose key damping
    
    b.damping = 1E-4

    p7 = TRANSLATE (c, (0, 0, 0))
    FIX_DIRECTION (b, p7, (0, 0, -1))
 
#----------------------------------------------------------------------

# solver and run

GEOMETRIC_EPSILON (1e-6) # Use 100 to 10000 times smaller than the smallest characteristic geometrical feature of a model

#(i.e. initial clearances) smallest initial clearances (loose key / keyway gaps) approx 1.0mm therefore 0.001/1000=1e-6

#slv = GAUSS_SEIDEL_SOLVER (1E-3, 300, 1E-6)

slv = NEWTON_SOLVER(delta = 5E-6) #All default values used

OUTPUT (solfec, 2E-3) # The physical tests recorded digitased outputs at 2E-3s intervals

if RANK() == 0 and solfec.mode == 'WRITE':

  print 'Running', stop, 'seconds of analysis with step', step, '...'

t0 = time.time()

if solfec.mode <> 'READ':

  RUN (solfec, slv, stop)

elapsed = time.time() - t0
   
if RANK() == 0 and solfec.mode == 'WRITE':   
      
    print "analysis run time =", elapsed/3600.0, "hours"
