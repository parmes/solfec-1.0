# SOLFEC parallel scaling example
# -------------------------------------------------------
# Array of elastic cubes undergoing sine sweep excitation

import os, sys
def where(program):
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, program)):
      return path
  return None
d0 = where('solfec')
if d0 == None:
  print 'ERROR: solfec not found in PATH!'
  print '       Download and compile solfec; add solfec directory to PATH variable;'
  sys.exit(1)
sys.path.append(d0+'/scripts')
from acc_sweep import *

M = 5 # array edge size
N = 1 # cube mesh size
gap = 0.001 # between cubes
step = 5E-4 # time step
stop = 5 # duration
kifo = 'BC' # kinematics
outi = 0.003 # output interval
weak = 'OFF' # week scaling
solv = 'NS' # solver
lofq = 1 # acc sweep low freq.
hifq = 5 # acc sweep high freq.
amag = 1 # acc sweep magnitude
xdmf = 'OFF' # export XMDF
argv = NON_SOLFEC_ARGV()

# print help
if argv != None and ('-help' in argv or '-h' in argv):
  print '------------------------------------------------------------------------'
  print 'Array of cubes excitation parameters:'
  print '------------------------------------------------------------------------'
  print '-M number --> array edge size (default: %d)' % M
  print '-N number --> cube mesh edge size (default: %d)' % N
  print '-kifo name --> kinematics in {TL, BC, PR, RG} (default: %s)' % kifo
  print '               where: TL -- Total Lagrangian FEM'
  print '                      BC -- Body Co-rotational FEM'
  print '                      PR -- Pseudo-rigid'
  print '                      RG -- Rigid'
  print '-solv name --> solver in {NS, GS} (default: %s)' % solv
  print '               where: NS -- Projected Newton solver'
  print '               where: GS -- Gauss-Seidel solver'
  print '-outi number --> output interval (default: %g)' % outi
  print '-weak --> enable weak scaling test (default: %s)' % weak
  print '          in this mode a constant MxMxM array size'
  print '          per MPI rank is approximately maintained'
  print '-step number --> time step (default: %g)' % step
  print '-stop number --> duration (default: %g)' % stop
  print '-xdmf --> export XDMF in READ mode (default: %s)' % xdmf
  print '-help or -h --> show this help and exit'
  print 
  print 'NOTE: becasue the output path depends on the input parameters'
  print '      use the same parameters to access results in READ mode, e.g.'
  print '      solfec -v path/to/array-of-cubes.py {same parameters}, or'
  print '      use the output directory as an input path instead, e.g.'
  print '      solfec -v path/to/results/directory'
  print '------------------------------------------------------------------------'

# parse command line switches
if argv != None:
  for i in range (0,len(argv)):
    if argv [i] == '-M':
      M = int (argv [i+1])
    elif argv [i] == '-N':
      N = int (argv [i+1])
    elif argv [i] == '-kifo':
      if argv [i+1] in ('TL', 'BC', 'PR', 'RG'):
	kifo = argv [i+1]
    elif argv [i] == '-solv':
      if argv [i+1] in ('NS', 'GS'):
	solv = argv [i+1]
    elif argv [i] == '-outi':
      outi = float (argv [i+1])
    elif argv [i] == '-step':
      step = float (argv [i+1])
    elif argv [i] == '-stop':
      stop = float (argv [i+1])
    elif argv [i] == '-weak':
      weak = 'ON'
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
begining = 'out/array-of-cubes/'
ending = 'STE%g_DUR%g_%s_%s_M%d_N%d_%s%d' % \
  (step,stop,kifo,solv,M,N,'S' if weak == 'OFF' else 'W',ncpu)
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
    elif x in ('TL','BC','PR','RG'): kifo = x
    elif x in ('NS','GS'): solv = x
    elif x[0] == 'M': M = int(x[1:])
    elif x[0] == 'N': N = int(x[1:])
    elif x[0] == 'S': ncpu = int(x[1:])
    elif x[0] == 'W':
      ncpu = int(x[1:])
      weak = 'ON'
    else:
      print 'ERROR: path ending', ending
      print '       format is invalid'
      sys.exit(1)
  print 'From', ending, 'read:'
  print '    ',
  print '(step, stop, kifo, solv, M, N, weak, ncpu) =',
  print '(%g, %g, %s, %s, %d, %d, %s, %d)' % \
         (step, stop, kifo, solv, M, N, weak, ncpu), '\n'

# bulk and surface materials
mat = BULK_MATERIAL (sol, young = 1E6, poisson = 0.25, density = 100)
SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.1)

# .1-wide cube corner nodes
nodes = [0.0, 0.0, 0.0,
         0.1, 0.0, 0.0,
	 0.1, 0.1, 0.0,
	 0.0, 0.1, 0.0,
	 0.0, 0.0, 0.1,
	 0.1, 0.0, 0.1,
	 0.1, 0.1, 0.1,
	 0.0, 0.1, 0.1]

# modify M if week scaling is 'ON'
if weak == 'ON':
  M = int((ncpu * M**3)**(1./3.))+1

# create an array of cubes
lst = []
for i in range (0,M):
  for j in range (0,M):
    for k in range (0,M):
      msh = HEX (nodes, N, N, N, 0, [1] * 6)
      TRANSLATE (msh, (i*(0.1+gap), j*(0.1+gap), k*(0.1+gap)))
      if kifo in ('TL', 'BC', 'PR'):
        if kifo == 'PR': bod = BODY (sol, 'PSEUDO_RIGID', msh, mat)
        else: bod = BODY (sol, 'FINITE_ELEMENT', msh, mat, form = kifo)
	bod.scheme = 'DEF_LIM' # semi-implicit time integration
	bod.damping = step # damp out free vibrations
      else: bod = BODY (sol, 'RIGID', msh, mat)
      lst.append(bod.id) # append list of body ids

# acceleration sweep signal
(vt, vd, vv, va) = acc_sweep (.2*stop, step, .8*stop, lofq, hifq, amag)
tsv = [None]*(len(vt)+len(vd))
tsv[::2] = vt
tsv[1::2] = vv
tsv = TIME_SERIES (tsv)

# create outer moving walls
for i in (0, 1, 2):
  for dst in (-0.1, M*(0.1+gap)-gap):
    msh = HEX (nodes, 1, 1, 1, 0, [0] * 6)
    vec = (M*(1+gap),M*(1+gap),1)
    SCALE (msh, tuple(vec[i:]+vec[:i]))
    vec = (0, 0, dst)
    TRANSLATE (msh, tuple(vec[i:]+vec[:i]))
    bod = BODY (sol, 'OBSTACLE', msh, mat)
    if HERE(sol, bod):
      SET_VELOCITY (bod, bod.center, (1, 1, 0), tsv)

# set gravity
GRAVITY (sol, (0, 0, -10))

# create solver
if solv == 'NS': slv = NEWTON_SOLVER ()
else: slv = GAUSS_SEIDEL_SOLVER (1.0, 1000, meritval=1E-8)

# output interval
OUTPUT (sol, outi)

# run simulation
import time
start = time.clock()
RUN (sol, slv, stop)
if RANK() == 0 and sol.mode == 'WRITE':
  dt = time.clock() - start
  ln = '%s = %g' % (ending,dt)
  fp = begining + 'RUNTIMES'
  with open(fp, "a") as f: f.write(ln+'\n')
  print 'Runtime line:',  ln, 'appended to:', fp

# READ mode post-processing
if sol.mode == 'READ' and not VIEWER():
  if xdmf == 'ON':
    xdmf_path = 'out/array-of-cubes/'+ 'XDMF_' + ending if SUBDIR() == None \
		else 'out/array-of-cubes/' + 'XDMF_'+ ending + '/' + SUBDIR()
    XDMF_EXPORT (sol, (0.0, stop), xdmf_path, lst)

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
