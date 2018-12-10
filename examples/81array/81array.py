# Model name 81array.py
# Abaqus input file = 81_Brick_Model_14.inp - Revised pitch and gaps as per test report after verification completed 26/09/11
# 81_Brick_Model_14.inp also includes modified FB and IB Young's modulus FB = 0.138GPa, IB = 0.138GPa, LK = 11.8GPa.
# Array of bricks with brick normal gaps to match those in test, DYNAMIC run all FEM bodies - Loose Key facets removed
# VELOCITY Time-history used from 5101060/23/03/08 - 2s, 0.3g, 3Hz sine dwell, then 3Hz to 10Hz linear sweep at 0.1Hz/s

import os
import sys
import gzip
import math
import time
import pickle
import commands
from math import cos 

def where(program):
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, program)):
      return path
  return None

path = where('solfec')

if path == None:
  print 'ERROR: solfec not found in PATH!'
  print '       Download and compile solfec; add solfec directory to PATH variable;'
  sys.exit(1)

print '(Found solfec at:', path + ')'

sys.path.append(os.path.join (path, 'scripts'))
from abaqusreader import AbaqusInput

dirpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dirpath)

# User paramters

argv = NON_SOLFEC_ARGV()

formu = 'BC'
fbmod = 64
ibmod = 64
afile = os.path.join (dirpath, '81array.inp')
step = 1E-4
damp = 1E-7
rest = 0.0
outi = 2E-2 # The physical tests recorded digitased outputs at 2E-3s intervals
stop = 72.0
genbase = False

if (argv == None and RANK() == 0) or (argv != None and '-help' in argv):
  print '------------------------------------------------------------------------'
  print '81 array example parameters:'
  print '------------------------------------------------------------------------'
  print '-form name => where name is TL, BC, RO, MODAL, PR or RG (default: %s)' % formu
  print '-fbmod num => fuel brick modes (default: %d)' % fbmod
  print '-ibmod num => interstitial brick modes (default: %d)' % ibmod
  print '-afile path => ABAQUS 81 array *.inp file path (default: %s)' % afile
  print '-step num => time step (default: %g)' % step
  print '-damp num => damping value (default: %g)' % damp
  print '-rest num => impact restitution, >= 0, <= 1 (default: %g)' % rest
  print '-outi num => output interval (default: %g)' % outi
  print '-stop num => sumulation end (default: %s)' % stop
  print '-genbase => generate RO (read mode) or MODAL (write mode) bases and stop'
  print '-help => show this help and exit'
  print '------------------------------------------------------------------------'

if argv != None:
  for i in range (0, len(argv)):
    if argv [i] == '-fbmod':
      fbmod = max (long (argv [i+1]), 6)
    elif argv [i] == '-ibmod':
      ibmod = max (long (argv [i+1]), 6)
    elif argv [i] == '-form':
      if argv [i+1] in ('TL', 'BC', 'RO', 'MODAL', 'PR', 'RG'):
	formu = argv [i+1]
    elif argv [i] == '-afile':
      afile = argv [i+1]
    elif argv [i] == '-step':
      step = float (argv [i+1])
    elif argv [i] == '-damp':
      damp = max (float (argv [i+1]), 0.0)
    elif argv [i] == '-rest':
      rest = max (min (1.0, float (argv [i+1])), 0.0)
    elif argv [i] == '-outi':
      outi = float (argv [i+1])
    elif argv [i] == '-stop':
      stop = float (argv [i+1])
    elif argv [i] == '-genbase':
      genbase = True
    elif argv [i] == '-help':
      sys.exit(0)
    elif argv [i][0] == '-':
      print 'INFO: invalid parameter: %s' % argv[i]
      print '                    try: -help'
      sys.exit(0)

if RANK() == 0:
  print 'ABAQUS file: %s'%afile
  print 'Formulation: %s'%formu
  if formu in ['MODAL', 'RO']:
    print '%d modes per fuel brick'%fbmod
    print '%d modes per interstitial brick'%ibmod
    print 'loose key are modeled as BC-FEM'
  print '%g damping'%damp
  print '%g restitution'%rest
  print '%g step'%step
  print '%g stop'%stop
  print '%g output interval'%outi
  if genbase and formu in ['TL', 'BC']:
    print 'generating RO or MODAL base enabled'
  print '----------------------------------------------------------'

if formu in ['RO', 'MODAL']:
  ending = '%s_FB%d_IB%d'%(formu,fbmod,ibmod)
else: ending = formu

ending = '%s_%s_s%.1e_d%.1e_r%g'%(afile [afile.rfind ('/')+1:len(afile)].replace ('.inp',''), ending, step, damp, rest)

# Analysis inputs

input_bricks = (['FB1(0)(0)', 'FB1(0)(1)', 'FB1(0)(2)', 'FB1(0)(3)', 'FB1(0)(4)', 'FB1(0)(5)', #left side
                'FB1(1)(0)', 'FB1(2)(0)', 'FB1(3)(0)', 'FB1(4)(0)', #top
                'FB1(1)(5)', 'FB1(2)(5)', 'FB1(3)(5)', 'FB1(4)(5)', #bottom
                'FB1(5)(0)', 'FB1(5)(1)', 'FB1(5)(2)', 'FB1(5)(3)', 'FB1(5)(4)', 'FB1(5)(5)', #right
                'IB2(0)(0)', 'IB2(0)(1)', 'IB2(0)(2)', 'IB2(0)(3)', 'IB2(0)(4)', #Side 1 (left)
                'IB1(0)(5)', 'IB1(1)(5)', 'IB1(2)(5)', 'IB1(3)(5)', 'IB1(4)(5)', #Side 2 (top)
                'IB2(5)(0)', 'IB2(5)(1)', 'IB2(5)(2)', 'IB2(5)(3)', 'IB2(5)(4)', #Side 3 (right)
                'IB1(0)(0)', 'IB1(1)(0)', 'IB1(2)(0)', 'IB1(3)(0)', 'IB1(4)(0)']) #Side 4 (bottom)

dwell = 2.0 # length in seconds of constant-frequency dwell at start of analysis

solfec = SOLFEC ('DYNAMIC', step, 'out/' + ending)

if solfec.mode == 'READ' and formu in ['MODAL', 'RO', 'RG', 'PR'] and genbase:
  print 'WARNING: -genbase was used with invalid formulation in READ mode --> exiting ...'
  sys.exit(0)

OUTPUT (solfec, outi) # The physical tests recorded digitased outputs at 2E-3s intervals

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.1, restitution = rest)

if RANK () == 0:
  commands.getoutput ('bunzip2 ' + os.path.join(dirpath,'ts81.py.bz2')) # only first CPU unpacks the input

BARRIER () # let all CPUs meet here

from ts81 import TS81 # import the time series

#This is a velocity-time series, for a 0.3g acceleration sine input at 3Hz, 2s dwell period, the 3Hz to 10Hz dwell sweep at 0.1Hz/s
vel = TS81()

if RANK () == 0:
  commands.getoutput ('bzip2 ' + os.path.join(dirpath,'ts81.py')) # only first CPU pack the input

# Create a new AbaqusInput object from the .inp deck:
model = AbaqusInput(afile, solfec)

# read RO bases once
if formu == 'RO':
  for label in ['FB1', 'FB2', 'IB1', 'IB2']:
    path = afile.replace ('.inp','_' + label + '_base.pickle.gz')
    try:
      robase = pickle.load(gzip.open(path, 'rb'))
      REGISTER_BASE (solfec, robase, label)
    except:
      print 'Reading %s failed --> run BC analysis in WRITE and READ modes first' % path
      sys.exit(0)

if formu == 'MODAL' and not genbase:
  modalbase = {}
  for label in ['FB1', 'FB2', 'IB1', 'IB2']:
    basepath = solfec.outpath + '/' + label + '_modalbase'
    try:
      print 'Opening:', basepath + '.h5'
      f = open(basepath + '.h5', 'rb')
      f.close()
      modalbase [label] = MODAL_ANALYSIS (path = basepath)
    except:
      print 'Reading: %s.h5 has failed:' % basepath
      print 'Run serial analysis with "-formu MODAL -genbase" switches first'
      sys.exit(0)

# Create a Finite Element body for each Instance in the Assembly:
for inst in model.assembly.instances.values():	# .instances is a dict
  label = inst.name	              # use Abaqus instance name
  mesh = inst.mesh	              # solfec MESH object at the instance position
  bulkmat = inst.material	        # solfec BULK_MATERIAL object
  if formu == 'RG':
    bdy = BODY(solfec, 'RIGID', mesh, bulkmat, label)
  elif formu == 'PR':
    bdy = BODY(solfec, 'PSEUDO_RIGID', mesh, bulkmat, label)
  elif formu in ['TL', 'BC']:
    bdy = BODY(solfec, 'FINITE_ELEMENT', mesh, bulkmat, label, form = formu)
  elif formu == 'MODAL':
    if NCPU(solfec) == 1 and genbase:
      path = solfec.outpath + '/' + label[0:3] + '_modalbase'
      nm = {'FB':fbmod, 'IB':ibmod}
      if label[0:2] != 'LK' and os.path.exists(path + '.h5') == False:
        bdy = BODY(solfec, 'FINITE_ELEMENT', COPY (mesh), bulkmat, label)
	MODAL_ANALYSIS (bdy, nm[label[0:2]], path, abstol = 1E-13)
	print 'Saved MODAL base as: %s.h5' % path
    elif label[0:2] == 'LK': bdy = BODY(solfec, 'FINITE_ELEMENT', mesh, bulkmat, label, form = 'BC')
    else: bdy = BODY(solfec, 'FINITE_ELEMENT', mesh, bulkmat, label, form = 'BC-MODAL', base = modalbase[label[0:3]])
  elif formu == 'RO':
    if label[0:2] == 'LK': bdy = BODY(solfec, 'FINITE_ELEMENT', mesh, bulkmat, label, form = 'BC')
    else: bdy = BODY(solfec, 'FINITE_ELEMENT', mesh, bulkmat, label, form = 'BC-RO', base = label[0:3])

if solfec.mode == 'WRITE' and genbase:
  if formu != 'MODAL': print 'INFO: nothing happened --> use [-form MODAL -genbase] to generate modal bases'
  else: print 'INFO: -genbase was used to generate modal bases --> now run without [-genbase]; exiting...'
  sys.exit(0)
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

# statistics
counts = {'FB':0, 'IB':0, 'LK':0}
dofs = {'FB':0, 'IB':0, 'LK':0}
total_dofs = 0
for b in solfec.bodies:
  counts[b.label[0:2]] += 1
  dofs[b.label[0:2]] = len(b.velo)
  total_dofs += len(b.velo)
print
print '==============================================='
print 'BODY COUNTS:', counts
print 'BODY DOFS:', dofs
print 'TOTAL DOFS:', total_dofs
print '==============================================='
print

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

  # read displacement snapshots from a 'TL' or 'BC' simulation
  if formu in ['TL', 'BC'] and genbase:
    from sys import stdout
    import modred
    import numpy

    if modred.__version__[0] != '1':
      print 'modred version 1 is needed!'
      import sys
      sys.exit(1)

    fb1_defo = COROTATED_DISPLACEMENTS (solfec, 'FB1(2)(2)')
    fb2_defo = COROTATED_DISPLACEMENTS (solfec, 'FB2(2)(2)')
    ib1_defo = COROTATED_DISPLACEMENTS (solfec, 'IB1(2)(2)')
    ib2_defo = COROTATED_DISPLACEMENTS (solfec, 'IB2(2)(2)')
    dur = DURATION (solfec)
    print 'Sampling FEM displacements ...', '    ' , 
    SEEK (solfec, dur[0])
    skip = int((1.0/outi)/10.0)
    while solfec.time < dur[1]:
      FORWARD (solfec, skip, corotated_displacements='TRUE')
      print '\b\b\b\b\b%2d %%' % (100*solfec.time/dur[1]),
      stdout.flush ()
    print

    fb1 = BYLABEL (solfec, 'BODY', 'FB1(2)(2)')
    fb2 = BYLABEL (solfec, 'BODY', 'FB2(2)(2)')
    ib1 = BYLABEL (solfec, 'BODY', 'IB1(2)(2)')
    ib2 = BYLABEL (solfec, 'BODY', 'IB2(2)(2)')
    fb1_rig = RIGID_DISPLACEMENTS (fb1)
    fb2_rig = RIGID_DISPLACEMENTS (fb2)
    ib1_rig = RIGID_DISPLACEMENTS (ib1)
    ib2_rig = RIGID_DISPLACEMENTS (ib2)

    pod_input = [(fb1_rig, fb1_defo, 'FB1', fbmod),
                 (fb2_rig, fb2_defo, 'FB2', fbmod),
                 (ib1_rig, ib1_defo, 'IB1', ibmod),
                 (ib2_rig, ib2_defo, 'IB2', ibmod)]

    for (rig, defo, label, num_modes) in pod_input:
      vecs = numpy.transpose(numpy.array(rig+defo))
      svec = vecs.shape[0]
      nvec = vecs.shape[1]
      print '%s:' % label, 'calculating', num_modes, 'POD modes from', nvec, 'input vectors of size', svec, '...'
      try:
        modes, vals = modred.compute_POD_matrices_snaps_method(vecs, list(range(num_modes)))
      except:
	print 'POD generation has failed --> it is possible that you tried to extract too many modes'
	if label[0:2] == 'FB':
	  print '                              try using [-fbmod a_smaller_number] and re-run'
	else:
	  print '                              try using [-ibmod a_smaller_number] and re-run'
	sys.exit ()
      mod = numpy.transpose(modes).tolist()
      val = vals.tolist()
      basevec = [x for vec in mod for x in vec]
      podbase = (val[0:len(mod)], basevec)
      path = afile.replace ('.inp','_' + label + '_base.pickle.gz')
      print 'Saving: %s' % path
      pickle.dump(podbase, gzip.open(path,'wb'))
