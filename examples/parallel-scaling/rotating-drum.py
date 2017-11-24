# SOLFEC parallel scaling example
# -------------------------------------------------------
# Rotating drum filled with ellipsoids

from math import sqrt
from math import sin
from math import cos
import sys

# create drum side "wheel" wall
def WHEEL (x, y, z, r, t, v, s):
  a = 0.0
  points = []
  while a < 6.2:
    px = x + r * sin (a)
    py = y
    pz = z + r * cos (a)
    points.append (px)
    points.append (py)
    points.append (pz)
    points.append (px)
    points.append (py+t)
    points.append (pz)
    a += 0.2
  return HULL (points, v, s)

# create drum block
def BLOCK (x, y, z, wx, wy, wz, v, s):
  points = [x-.5*wx, y-.5*wy, z-.5*wz,
            x+.5*wx, y-.5*wy, z-.5*wz,
            x+.5*wx, y+.5*wy, z-.5*wz,
            x-.5*wx, y+.5*wy, z-.5*wz,
            x-.5*wx, y-.5*wy, z+.5*wz,
            x+.5*wx, y-.5*wy, z+.5*wz,
            x+.5*wx, y+.5*wy, z+.5*wz,
            x-.5*wx, y+.5*wy, z+.5*wz]
  return HULL (points, v, s)

# create particles
def PARTICLES (rx, ry, rz, x0, y0, z, nz, wx, wy, v, s):
  global ipar, npar, sol, mat, sphs, kifo, step
  if ipar == npar: return
  rr = (rx+ry+rz)/3.
  if sphs == 'ON':
    rx = rr
    ry = rr
    rz = rr
  a = x0 + .5*wx
  b = y0 + .5*wy
  iz = 0
  while iz < nz:
    iz = iz + 1
    x = x0 - .5*wx
    while x < a:
      y = y0 - .5*wy
      while y < b:
	if ipar < npar:
	  if RANK() == 0:
	    lb = 'PAR%d' % ipar
	    if kifo in ('PR','RG'):
	      if sphs == 'ON':
		shp = SPHERE ((x, y, z), rr, v, s)
	      else: shp = ELLIP ((x, y, z), (rx, ry, rz), v, s)
	      if kifo == 'RG': BODY (sol, 'RIGID', shp, mat, lb)
	      else:
		bod = BODY (sol, 'PSEUDO_RIGID', shp, mat, lb)
		bod.scheme = 'DEF_LIM' # semi-implicit time integration
		bod.damping = step # damp out free vibrations
	    else:
	      msh = ELLIP_MESH ((x, y, z), (rx, ry, rz), rr*0.1, v, s)
	      bod = BODY (sol, 'FINITE_ELEMENT', msh, mat, lb, form = 'BC')
	      bod.scheme = 'DEF_LIM' # semi-implicit time integration
	      bod.damping = step # damp out free vibrations
	  ipar = ipar + 1
	y = y + 2.*ry
      x = x + 2.*rx
    z = z - 2.*rz

# paramters 
npar = 100 # number of particles
step = 0.001 # time step
outi = 0.03 # output interval
stop = 10 # duration
fric = 0.3 # friction coefficient
angv = 1.0 # drum angular velocity
ipar = 0 # partical counter (used internally)
kifo = 'PR' # kinematics
solv = 'NS' # constraint solver
maxi = 1000 # constraint solver maximum iterations
nsdl = 0.0 # Newton solver delta
leps = 0.25 # Newton solver epsilon
lmxi = 10 # Newton solver linmaxiter
weak = 'OFF' # weak scaling test flag
sphs = 'OFF' # use sphereical particles
prfx = '' # predix string
xdmf = 'OFF' # export XDMF (-kifo FE)
argv = NON_SOLFEC_ARGV()

# print help
if argv != None and ('-help' in argv or '-h' in argv):
  print '------------------------------------------------------------------------'
  print 'Rotating drum with ellipsoidal or sphereical particles:'
  print '------------------------------------------------------------------------'
  print '-npar number --> number of particles (default: %d)' % npar
  print '-kifo name --> kinematics in {FE, PR, RG} (default: %s)' % kifo
  print '               where: FE -- Finite Element'
  print '                      PR -- Pseudo-rigid'
  print '                      RG -- Rigid'
  print '-solv name --> Constraint solver in {NS, GS} (default: %s)' % solv
  print '               where: NS -- Projected Newton solver'
  print '               where: GS -- Gauss-Seidel solver'
  print '-maxi number --> Constraint solver max. iter. (default: %d)' % maxi
  print '-nsdl number --> Newton solver delta (default: %g)' % nsdl
  print '-leps number --> Newton linear solver epsilon (default: %g)' % leps
  print '-lmxi number --> Newton linear solver linmaxiter (default: %d)' % lmxi
  print '-outi number --> output interval (default: %g)' % outi
  print '-weak --> enable weak scaling test (default: %s)' % weak
  print '          in this mode a "-npar" particles per'
  print '          MPI rank is approximately maintained'
  print '-step number --> time step (default: %g)' % step
  print '-stop number --> duration (default: %g)' % stop
  print '-fric number --> friction coefficient (default %g)' % fric
  print '-angv number --> drum angular velocity [rad/s] (default %g)' % angv
  print '-sphs --> use spherical particles (default %s)' % sphs
  print '-prfx string --> include a prefix string into the output path'
  print '-xdmf --> export XDMF in READ mode (-kifo FE; default: %s)' % xdmf
  print '-help or -h --> show this help and exit'
  print 
  print 'FYI: smaller time step may be appropriate for larger perticle numbers'
  print '     since the size of individual particles gets proportionally smaller'
  print
  print 'NOTE: because the output path depends on the input parameters'
  print '      use the same parameters to access results in READ mode, e.g.'
  print '      solfec -v path/to/rotating-drum.py {same parameters}, or'
  print '      use the output directory as an input path instead, e.g.'
  print '      solfec -v path/to/results/directory'
  print '------------------------------------------------------------------------'

# parse command line switches
if argv != None:
  for i in range (0,len(argv)):
    if argv [i] == '-npar':
      npar = int (argv [i+1])
    elif argv [i] == '-kifo':
      if argv [i+1] in ('FE', 'PR', 'RG'):
	kifo = argv [i+1]
    elif argv [i] == '-solv':
      if argv [i+1] in ('NS', 'GS'):
	solv = argv [i+1]
    elif argv [i] == '-maxi':
      maxi = int (argv [i+1])
    elif argv [i] == '-nsdl':
      nsdl = float (argv [i+1])
    elif argv [i] == '-leps':
      leps = float (argv [i+1])
    elif argv [i] == '-lmxi':
      lmxi = int (argv [i+1])
    elif argv [i] == '-outi':
      outi = float (argv [i+1])
    elif argv [i] == '-step':
      step = float (argv [i+1])
    elif argv [i] == '-stop':
      stop = float (argv [i+1])
    elif argv [i] == '-fric':
      fric = float (argv [i+1])
    elif argv [i] == '-angv':
      angv = float (argv [i+1])
    elif argv [i] == '-weak':
      weak = 'ON'
    elif argv [i] == '-sphs':
      sphs = 'ON'
    elif argv [i] == '-prfx':
      prfx = str (argv [i+1])
    elif argv [i] == '-xdmf':
      xdmf = 'ON'
    elif argv [i] in ('-help', '-h'):
      sys.exit(0)
    elif argv [i][0] == '-':
      print 'INFO: invalid parameter: %s' % argv[i]
      print '                    try: -help'
      sys.exit(0)

# number of CPUs
ncpu = NCPU ()

# output path components
begining = 'out/rotating-drum/'
if len(prfx) > 0: prfx += '_'
solvstr = 'NS_MAXI%d_NSDL%g_LEPS%g_LMXI%d' % (maxi, nsdl, leps, lmxi) if solv == 'NS' else 'GS_MAXI%d' % maxi
ending = prfx + 'STE%g_DUR%g_%s_%s_%s_FRI%g_ANG%g_N%d_%s%d' % \
  (step,stop,kifo,solvstr,'ELL' if sphs == 'OFF' else 'SPH',\
  fric,angv,npar,'S' if weak == 'OFF' else 'W',ncpu)
outpath = begining + ending

# create solfec object
sol = SOLFEC ('DYNAMIC', step, outpath)

# test whether command line switches need to
# be deduced from the output path in READ mode
if sol.mode == 'READ' and sol.outpath != outpath:
  ending = sol.outpath[sol.outpath.rfind('/')+1:]
  for x in ending.split('_'):
    if x[0:3] == 'STE': step = float(x[3:])
    elif x[0:3] == 'DUR': stop = float(x[3:])
    elif x[0:3] == 'FRI': fric = float(x[3:])
    elif x[0:3] == 'ANG': angv = float(x[3:])
    elif x in ('FE','PR','RG'): kifo = x
    elif x in ('NS','GS'): solv = x
    elif x[0:4] == 'MAXI': maxi = int(x[4:])
    elif x[0:4] == 'NSDL': nsdl = float(x[4:])
    elif x[0:4] == 'LEPS': leps = float(x[4:])
    elif x[0:4] == 'LMXI': lmxi = int(x[4:])
    elif x == 'ELL': sphs = 'OFF'
    elif x == 'SPH': sphs = 'ON'
    elif x[0] == 'N': N = int(x[1:])
    elif x[0] == 'S': ncpu = int(x[1:])
    elif x[0] == 'W':
      ncpu = int(x[1:])
      weak = 'ON'
    elif x != ending[:ending.index("_STE")]:
      # or if not a prefix report an error
      print 'ERROR: path ending', ending
      print '       format is invalid'
      sys.exit(1)
  print 'From', ending, 'read:'
  print '    ',
  print '(step, stop, kifo, solv, sphs, fric, npar, weak, ncpu) =',
  print '(%g, %g, %s, %s, %s, %g, %d, %s, %d)' % \
         (step, stop, kifo, solv, sphs, fric, npar, weak, ncpu), '\n'

# bulk and surface materials
mat = BULK_MATERIAL (sol, young = 1E6, poisson = 0.25, density = 100)
SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = fric)

# modify npar if week scaling is 'ON'
if weak == 'ON':
  npar = int((ncpu * npar))

# set gravity
GRAVITY (sol, (0, 0, -10))

# create solver
if solv == 'NS': slv = NEWTON_SOLVER (maxiter = maxi,\
  delta = nsdl, epsilon = leps, linmaxiter = lmxi, maxmatvec = 100000)
else: slv = GAUSS_SEIDEL_SOLVER (1.0, maxi, meritval=1E-8)

# output interval
OUTPUT (sol, outi)

# rotating drum
pip = PIPE ((0, 0, 0), (0, 0.5, 0), 1, 0.05, 1, 32, 1, 1, [1, 1, 1, 1])
bl1 = BLOCK(-0.9, 0.25, 0, 0.2, 0.5, 0.2, 1, 1)
bl2 = BLOCK (0.9, 0.25, 0, 0.2, 0.5, 0.2, 1, 1)
bod = BODY (sol, 'OBSTACLE', [pip, bl1, bl2], mat)

# side walls
wh1 = WHEEL (0, -0.002, 0, 1.05, -0.05, 1, 1)
BODY (sol, 'OBSTACLE', wh1, mat)
wh2 = WHEEL (0, 0.502, 0, 1.05, 0.05, 1, 1)
BODY (sol, 'OBSTACLE', wh2, mat)

# ellipsoid radii
rr = (.1/npar)**(1./3.)
rx = .9*rr
ry = .5*rr
rz = rr
nn = int(0.75/(2.*rz))
# .5*acc*t^2 = 2*rz+eps --> t = sqrt(2*(2*rr+eps)/acc)
dt = nn*sqrt(2*(2*rz+.1*rz)/10)

# simulation callback
rotating = False
def callback (bod):
  global rotating, ipar, npar
  if ipar < npar: # insert particles
    PARTICLES (rx, ry, rz, 0, 0.25, 0.25, nn, 1.4, 0.4, 2, 2)
  elif not rotating: # rotate drum
    INITIAL_VELOCITY (bod, (0, 0, 0), (0, angv, 0))
    rotating = True
  return 1

# set callback
CALLBACK (sol, dt, bod, callback)

# run simulation
import time
slv.itershist = 'ON'
start = time.clock()
RUN (sol, slv, stop)
if RANK() == 0 and sol.mode == 'WRITE':
  dt = time.clock() - start
  ln = '%s = %g' % (ending,dt)
  fp = begining + 'RUNTIMES'
  with open(fp, "a") as f: f.write(ln+'\n')
  print 'Runtime line:',  ln, 'appended to:', fp
  itershist = slv.itershist
  itsum = sum(itershist)
  itup = sum([x for x in itershist if x > 0 and x < maxi])
  itlo = sum([1 for x in itershist if x > 0 and x < maxi])
  itavg = itup / itlo if itup > 0 else 1
  nmaxi = itershist.count(maxi)
  ln = '%s = %d, %d, %d' % (ending,itsum,itavg,nmaxi)
  fp = begining + 'ITERS'
  with open(fp, "a") as f: f.write(ln+'\n')
  print 'Runtime line:',  ln, 'appended to:', fp

# READ mode post-processing
if sol.mode == 'READ' and not VIEWER():
  if kifo == 'FE' and xdmf == 'ON':
    xdmf_path = 'out/rotating-drum/'+ 'XDMF_' + ending if SUBDIR() == None \
		else 'out/rotating-drum/' + 'XDMF_'+ ending + '/' + SUBDIR()
    XDMF_EXPORT (sol, (0.0, stop), xdmf_path, [bod, 'PAR'])

  timer = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL', 'PARBAL']
  dur = DURATION (sol)
  th = HISTORY (sol, timer, dur[0], dur[1])
  timing = {'TOTAL':0.0}
  for i in range (0, len(timer)):
    ss = 0.0
    for tt in th [i+1]: ss += tt
    timing[timer[i]] = ss
    timing['TOTAL'] += ss
  ln = '%s = %s' % (ending, timing)
  fp = begining + 'TIMINGS'
  with open(fp, "a") as f: f.write(ln+'\n')
  print 'Timing line:',  ln, 'appended to:', fp
