# Two fuel bricks normal contact restitution test
import os
import sys
import math
import gzip
import numpy
import modred
import pickle
import commands
sys.path.append('scripts')
sys.path.append('examples/81array')
from abaqusreader import AbaqusInput
from math import cos 

# User paramters

argv = NON_SOLFEC_ARGV()

formu = 'BC'
fbmod = 24
afile = 'examples/81array/81fbi.inp'
step = 1E-4
damp = 1E-7
impactVelocity = 0.15 #15cm/s
rest = 0.0
stop = 0.5
genbase = False

if argv == None or '-help' in argv:
  print '---------------------------------------------------------------------'
  print '81 fuel brick impact example parameters:'
  print '---------------------------------------------------------------------'
  print '-form name => where name is TL, BC, RO, MODAL, PR or RG (default: %s)' % formu
  print '-fbmod num => fuel brick modes (default: %d)' % fbmod
  print '-damp num => damping, >= 0.0 (default: %g)' % damp
  print '-step num => time step, > 0.0 (default: %g)' % step
  print '-afile path => ABAQUS fuel brick impact *.inp file path (default: %s)' % afile
  print '-rest num => impact restitution, >= 0, <= 1 (default: %g)' % rest
  print '-genbase => generate RO base (default: same as for 81array.py)'
  print '-help => show this help and exit'
  print '---------------------------------------------------------------------'

if argv != None:
  for i in range (0, len(argv)):
    if argv [i] == '-fbmod':
      fbmod = max (long (argv [i+1]), 6)
    elif argv [i] == '-form':
      if argv [i+1] in ('TL', 'BC', 'RO', 'MODAL', 'PR', 'RG'):
	formu = argv [i+1]
    elif argv [i] == '-damp':
      damp = max (float (argv [i+1]), 0.0)
    elif argv [i] == '-step':
      step = float (argv [i+1])
      if step <= 0.0: step = 1E-4
    elif argv [i] == '-afile':
      afile = argv [i+1]
    elif argv [i] == '-rest':
      rest = max (min (1.0, float (argv [i+1])), 0.0)
    elif argv [i] == '-genbase':
      genbase = True
    elif argv [i] == '-help':
      sys.exit(0)
    elif argv [i][0] == '-':
      print 'INFO: invalid parameter: %s' % argv[i]
      print '                    try: -help'
      sys.exit(0)

print 'ABAQUS file: %s'%afile
print 'formulation: %s'%formu
if formu in ['MODAL', 'RO']:
  print '%d modes per fuel brick'%fbmod
print '%g damping'%damp
print '%g restitution'%rest
print '%g step'%step
if genbase and formu in ['TL', 'BC']:
  print 'generating RO base enabled'
print '------------------------------------------------------'

# Model

if formu in ['RO', 'MODAL']:
  ending = '%s_FB%d'%(formu, fbmod)
else: ending = formu

ending = '%s_%s_s%.1e_d%.1e_r%g'%(afile [afile.rfind ('/')+1:len(afile)].replace ('.inp',''), ending, step, damp, rest)

solfec = SOLFEC ('DYNAMIC', step, 'out/' + ending)

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.1, restitution = rest)

# Create a new AbaqusInput object from the .inp deck:
model = AbaqusInput(afile, solfec)

# read RO bases once
arraybase = False
if formu == 'RO':
  robase = {}
  for label in ['FB1', 'FB2']:
    path0 = solfec.outpath + '/' + label + '_base.pickle.gz'
    try:
      robase[label] = pickle.load(gzip.open(path0, 'rb'))
      print 'RO --> using %s base' % path0
    except:
      print 'Reading %s failed --> you can run -form BC -genbase to genrate this file' % path0
      path1 = afile.replace('fbi','array').replace ('.inp','_' + label + '_base.pickle.gz')
      print 'Now trying to use %s instead ...' % path1
      try:
        robase[label] = pickle.load(gzip.open(path1, 'rb'))
	arraybase = True
      except:
	print 'Reading %s failed --> run 81array.py -afile %s analysis in WRITE and READ modes first' % (path1, afile.replace('fbi', 'array'))
	sys.exit(0)

if arraybase: # 81 array base is rotated by 45 deg around the (0, 0, 1) direction
  for label in robase:
    base = robase[label]
    val0 = base[0]
    vec0 = base[1]
    vec1 = []
    for i in range(0, len(vec0), 3):
      x = (vec0[i], vec0[i+1], vec0[i+2])
      y = ROTATE (x, (0, 0, 0), (0, 0, 1), -45)
      vec1 += y
    robase[label] = (val0, vec1) # POD base rotated by -45 deg about (0, 0, 1)

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
    bdy = BODY(solfec, 'FINITE_ELEMENT', COPY (mesh), bulkmat, label)
    data = MODAL_ANALYSIS (bdy, fbmod, solfec.outpath + '/modal' + label, abstol = 1E-13)
    DELETE (solfec, bdy)
    bdy = BODY(solfec, 'FINITE_ELEMENT', mesh, bulkmat, label, form = 'BC-MODAL', base = data)
  elif formu == 'RO':
    bdy = BODY(solfec, 'FINITE_ELEMENT', mesh, bulkmat, label, form = 'BC-RO', base = robase[label[0:3]])

# boundary conditions
for b in solfec.bodies:
  
  if b.kind != 'RIGID':
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

  res_vel_temp = zip(th[3], th[1])
  res_vel = [r1-r2 for r1,r2 in res_vel_temp]
			    
  vrest = []
  for v in res_vel:
    vrest.append(v/impactVelocity)

  print "For FB damping =", b.damping, " coeficient of restitution = ", vrest[-1]
    
  try:
    import matplotlib.pyplot as plt

    plt.plot (th[0],th[1],lw = 2, label = 'Left_Brick_Top_Point_Velocity_VX_(m/s)')
    plt.plot (th[0],th[2],lw = 2, label = 'Left_Brick_Bottom_Point_Velocity_VX_(m/s)')
    plt.plot (th[0],th[3],lw = 2, label = 'Right_Brick_Top_Point_Velocity_VX_(m/s)')
    plt.plot (th[0],th[4],lw = 2, label = 'Right_Brick_Bottom_Point_Velocity_VX_(m/s)')
    plt.plot (th[0],vrest,lw = 2, label = 'Brick_Restitution')
    plt.legend(loc = 'lower right')
    plt.title ('SUB02 FB-FB Normal Contact FB damping=%s'% b.damping)
    plt.savefig (solfec.outpath + '/Brick_Velocity_Damp%g.png'% b.damping) 

  except ImportError:
    pass # no reaction
  
# solver and run
GEOMETRIC_EPSILON (1E-6)
slv = NEWTON_SOLVER() #All default values used
OUTPUT(solfec, 1E-3) 

# sample displacements
if formu in ['TL', 'BC'] and genbase:
  fb1 = BYLABEL (solfec, 'BODY', 'FB1(0)(0)')
  fb2 = BYLABEL (solfec, 'BODY', 'FB2(0)(0)')
  fb1_defo = COROTATED_DISPLACEMENTS (solfec, 'FB1(0)(0)')
  fb2_defo = COROTATED_DISPLACEMENTS (solfec, 'FB2(0)(0)')
  fb1_rig = RIGID_DISPLACEMENTS (fb1)
  fb2_rig = RIGID_DISPLACEMENTS (fb2)
        
if RANK() == 0 and solfec.mode == 'WRITE':
  print 'Running', stop, 'seconds of analysis with step', step, '...'

if solfec.mode <> 'READ':
  RUN (solfec, slv, stop)

if solfec.mode == 'READ' and formu in ['TL', 'BC'] and genbase:
  from sys import stdout
  print 'Sampling FEM displacements ...', '    ' , 
  SEEK (solfec, 0.0)
  while solfec.time < stop:
    FORWARD (solfec, 1, corotated_displacements='TRUE')
    print '\b\b\b\b\b%2d %%' % (100*solfec.time/stop),
    stdout.flush ()
  print

if formu in ['TL', 'BC'] and genbase:

  pod_input = [(fb1_rig, fb1_defo, 'FB1', fbmod),
	       (fb2_rig, fb2_defo, 'FB2', fbmod)]

  for (rig, defo, label, num_modes) in pod_input:
    vecs = numpy.transpose(numpy.array(rig+defo))
    svec = vecs.shape[0]
    nvec = vecs.shape[1]
    print '%s:' % label, 'calculating', num_modes, 'POD modes from', nvec, 'input vectors of size', svec, '...'
    try:
      modes, vals = modred.compute_POD_matrices_snaps_method(vecs, list(range(num_modes)))
    except:
      print 'POD generation has failed --> it is possible that you tried to extract too many modes'
      print '                              try using [-fbmod a_smaller_number] and re-run'
      sys.exit ()
    mod = numpy.transpose(modes).tolist()
    val = vals.tolist()
    basevec = [x for vec in mod for x in vec]
    podbase = (val[0:len(mod)], basevec)
    try: os.makedirs(solfec.outpath.replace(formu,'RO_FB%d'%fbmod))
    except: pass
    path = solfec.outpath.replace(formu,'RO_FB%d'%fbmod) + '/' + label + '_base.pickle.gz'
    pickle.dump (podbase, gzip.open(path,'wb'))
