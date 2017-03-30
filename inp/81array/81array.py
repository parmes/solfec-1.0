# Model name 81array.py
# Abaqus input file = 81_Brick_Model_14.inp - Revised pitch and gaps as per test report after verification completed 26/09/11
# 81_Brick_Model_14.inp also includes modified FB and IB Young's modulus FB = 0.138GPa, IB = 0.138GPa, LK = 11.8GPa.
# Array of bricks with brick normal gaps to match those in test, DYNAMIC run all FEM-RO bodies - Loose Key facets removed
# VELOCITY Time-history used from 5101060/23/03/08 - 2s, 0.3g, 3Hz sine dwell, then 3Hz to 10Hz linear sweep at 0.1Hz/s

import sys
import math
import sys
import time
import pickle
import commands
sys.path.append('scripts/abaqusreader')
sys.path.append('scripts')
sys.path.append('examples/81array')
from abaqusreader import AbaqusInput
from math import cos 

# User paramters

argv = NON_SOLFEC_ARGV()

if argv == None:
  print '----------------------------------------------------------'
  print 'No user paramters passed! Possible paramters:'
  print '----------------------------------------------------------'
  print '-form name => where name is TL, BC, RO, PR or RG'
  print '-fbmod num => fuel brick modes, num >= 6 and <= 64'
  print '-ibmod num => interstitial brick modes, num >= 6 and <= 64'
  print '-lkmod num => loose key modes, num >= 6 and <= 12'
  print '-afile path => Abaqus 81 array file path'
  print '-step num => time step'
  print '-damp num => damping value'
  print '-rest num => impact restitution'
  print '-outi num => output interval'
  print '----------------------------------------------------------'

formu = 'RO'
fbmod = 12
ibmod = 12
lkmod = 12
afile = 'examples/81array/81array.inp'
step = 1E-4
damp = 1E-5
rest = 0.0
outi = 2E-3 # The physical tests recorded digitased outputs at 2E-3s intervals
if argv != None and len (argv) > 1:
  for i in range (0, len(argv)-1):
    if argv [i] == '-fbmod':
      fbmod = max (min (64, long (argv [i+1])), 6)
    elif argv [i] == '-ibmod':
      ibmod = max (min (64, long (argv [i+1])), 6)
    elif argv [i] == '-lkmod':
      lkmod = max (min (12, long (argv [i+1])), 6)
    elif argv [i] == '-form':
      if argv [i+1] in ('TL', 'BC', 'RO', 'PR', 'RG'):
	formu = argv [i+1]
    elif argv [i] == '-afile':
      afile = argv [i+1]
    elif argv [i] == '-step':
      step = float (argv [i+1])
      if step <= 0.0: step = 1E-4
    elif argv [i] == '-damp':
      damp = max (float (argv [i+1]), 0.0)
    elif argv [i] == '-rest':
      rest = max (min (1.0, float (argv [i+1])), 0.0)
    elif argv [i] == '-outi':
      outi = float (argv [i+1])

print 'Using formulation: ', formu
if formu == 'RO':
  print '%d modes per fuel brick'%fbmod
  print '%d modes per interstitial brick'%ibmod
  print '%d modes per loose key'%lkmod
print 'Using %g time step and %g damping'%(step, damp)
print '----------------------------------------------------------'

if formu == 'RO':
  ending = 'RO_FB%d_IB%d_LK%d'%(fbmod,ibmod,lkmod)
else: ending = formu

ending = '%s_%s_s%.1e_d%.1e_r%g'%(afile [afile.rfind ('/'):len(afile)].replace ('.inp',''), ending, step, damp, rest)

# Analysis inputs

input_bricks = (['FB1(0)(0)', 'FB1(0)(1)', 'FB1(0)(2)', 'FB1(0)(3)', 'FB1(0)(4)', 'FB1(0)(5)', #left side
                'FB1(1)(0)', 'FB1(2)(0)', 'FB1(3)(0)', 'FB1(4)(0)', #top
                'FB1(1)(5)', 'FB1(2)(5)', 'FB1(3)(5)', 'FB1(4)(5)', #bottom
                'FB1(5)(0)', 'FB1(5)(1)', 'FB1(5)(2)', 'FB1(5)(3)', 'FB1(5)(4)', 'FB1(5)(5)', #right
                'IB2(0)(0)', 'IB2(0)(1)', 'IB2(0)(2)', 'IB2(0)(3)', 'IB2(0)(4)', #Side 1 (left)
                'IB1(0)(5)', 'IB1(1)(5)', 'IB1(2)(5)', 'IB1(3)(5)', 'IB1(4)(5)', #Side 2 (top)
                'IB2(5)(0)', 'IB2(5)(1)', 'IB2(5)(2)', 'IB2(5)(3)', 'IB2(5)(4)', #Side 3 (right)
                'IB1(0)(0)', 'IB1(1)(0)', 'IB1(2)(0)', 'IB1(3)(0)', 'IB1(4)(0)']) #Side 4 (bottom)

stop = 72.0
dwell = 2.0 # length in seconds of constant-frequency dwell at start of analysis

solfec = SOLFEC ('DYNAMIC', step, 'out/' + ending)

OUTPUT (solfec, outi) # The physical tests recorded digitased outputs at 2E-3s intervals

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.1, restitution = rest)

if RANK () == 0:
  commands.getoutput ("bunzip2 examples/81array/ts81.py.bz2") # only first CPU unpacks the input

BARRIER () # let all CPUs meet here

from ts81 import TS81 # import the time series

#This is a velocity-time series, for a 0.3g acceleration sine input at 3Hz, 2s dwell period, the 3Hz to 10Hz dwell sweep at 0.1Hz/s
vel = TS81()

if RANK () == 0:
  commands.getoutput ("bzip2 examples/81array/ts81.py") # only first CPU pack the input

# Create a new AbaqusInput object from the .inp deck:
model = AbaqusInput(afile, solfec)

# Create a Finite Element body for each Instance in the Assembly:
for inst in model.assembly.instances.values():	# .instances is a dict
  label = inst.name	              # use Abaqus instance name
  mesh = inst.mesh	              # solfec MESH object at the instance position
  bulkmat = inst.material	        # solfec BULK_MATERIAL object
  if formu == 'RG':
    bdy = BODY(solfec, 'RIGID', mesh, bulkmat, label)
  elif formu == 'PR':
    bdy = BODY(solfec, 'PSEUDO_RIGID', mesh, bulkmat, label)
  elif formu != 'RO':
    bdy = BODY(solfec, 'FINITE_ELEMENT', mesh, bulkmat, label, form = formu)
  else:
    bdy = BODY(solfec, 'FINITE_ELEMENT', COPY (mesh), bulkmat, label)
    if label.startswith ('FB'): nm = fbmod
    elif label.startswith ('IB'): nm = ibmod
    else: nm = lkmod
    data = MODAL_ANALYSIS (bdy, nm, solfec.outpath + '/modal' + label, abstol = 1E-13)
    DELETE (solfec, bdy)
    bdy = BODY(solfec, 'FINITE_ELEMENT', mesh, bulkmat, label, form = 'RO', modal = data)

#----------------------------------------------------------------------

# boundary conditions and input accelerations

for b in solfec.bodies:

  if b.kind != 'RIGID': b.scheme = 'DEF_LIM'
  c = b.center

  print "body mass:", b.label, b.mass, "Kg"
  
  if b.label.startswith('FB'): # find FBs
  
    # Set fuel brick damping
  
    b.damping = damp

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
    
    b.damping = damp

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
    
    b.damping = damp

    p7 = TRANSLATE (c, (0, 0, 0))
    FIX_DIRECTION (b, p7, (0, 0, -1))
 
#----------------------------------------------------------------------

# solver and run

GEOMETRIC_EPSILON (1e-6) # Use 100 to 10000 times smaller than the smallest characteristic geometrical feature of a model

#(i.e. initial clearances) smallest initial clearances (loose key / keyway gaps) approx 1.0mm therefore 0.001/1000=1e-6

#slv = GAUSS_SEIDEL_SOLVER (1E-3, 1000, 1E-8)

slv = NEWTON_SOLVER (delta = 5E-6) # Small diagonal regularisation for the pseudo-transient continuation

if RANK() == 0 and solfec.mode == 'WRITE':
  print 'Running', stop, 'seconds of analysis with step', step, '...'

t0 = time.time()

if solfec.mode <> 'READ':

  RUN (solfec, slv, stop)

elapsed = time.time() - t0
   
if RANK() == 0 and solfec.mode == 'WRITE':   
    print "analysis run time =", elapsed/3600.0, "hours"

if not VIEWER() and solfec.mode == 'READ': # extract and output time series

  # --- configurable parameters, which shouldn't be changed ---
  FBlabels = ['FB1(0)(0)',
  'FB1(3)(4)', # B
  'FB2(2)(3)', # E
  'FB2(1)(2)', # C
  'FB2(2)(2)', # T
  'FB1(3)(3)', # F
  'FB2(0)(2)', # G
  'FB2(1)(1)', # J
  'FB1(3)(1)', # Q
  'FB2(1)(3)', # K
  'FB1(3)(2)', # M
  ] # Fuel bricks results are required for

  IBlabels = ['IB2(3)(2)', # 7
  'IB2(1)(0)', # 18
  'IB2(4)(2)', # 13
  'IB2(4)(0)', # 17
  ] # Interstital bricks results are required for

  boundarylabel = 'FB1(0)(0)' # label of a boundary brick - NOTE this MUST also be the first entry in 'FBlabels'
  # --------------------------------------------------------------

  # build results request
  requests = []
  for bl in FBlabels:
    bdy = BYLABEL(solfec, 'BODY', bl)
    c = bdy.center

    # copied from input script
    p1 = TRANSLATE (c, (0.1,-0.1, 0))
    p2 = TRANSLATE (c, (-0.1, -0.1, 0))

    requests.append( (bdy, p1, 'VX') )
    requests.append( (bdy, p2, 'VX') )

  for bl in IBlabels:
    bdy = BYLABEL(solfec, 'BODY', bl)
    c = bdy.center

    # copied from input script
    p4 = TRANSLATE (c, (0.07, -0.07, 0))
    p5 = TRANSLATE (c, (-0.07, -0.07, 0))

    requests.append( (bdy, p4, 'VX') )
    requests.append( (bdy, p5, 'VX') )

  # extract velocity results from the end of the dwell period onwards
  thv = HISTORY(solfec, requests, dwell, stop, progress='ON')

  # save time series using pickle object
  f = open(solfec.outpath + '/%s.thv'%ending, 'w')
  pickle.dump (thv, f)
  f.close ()
