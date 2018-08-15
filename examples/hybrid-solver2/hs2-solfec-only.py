M = 5 # outer layers
N = 3 # inner layers
gap = 0.001 # betweeb bodies
step = 5E-4 # time step
stop = 5.0 # duration
lofq = 1 # low excitation frequency
hifq = 1 # high excitation freqency
amag = 1 # acceleration magnitude

# find path to parmec source directory in order
# to load the acceleration sweep signal script
import os, sys
def where(program):
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, program)):
      return path
  return None
path = where('parmec4')
if path == None:
  print 'ERROR: parmec4 not found in PATH!'
  print '       Download and compile parmec;',
  print 'add parmec directory to PATH variable;'
  sys.exit(1)
print '(Found parmec4 at:', path + ')'
sys.path.append(os.path.join (path, 'python'))

# generate acceleration sweep signal;
# note that parmec will use the associated
# velocity signal, rather than the acceleration itself
from acc_sweep import *
(vt, vd, vv, va) = acc_sweep (step, stop, lofq, hifq, amag)
tsv = [None]*(len(vt)+len(vd))
tsv[::2] = vt
tsv[1::2] = vv
tsv = TIME_SERIES (tsv) # velocity time series

# create solfec object
sol = SOLFEC ('DYNAMIC', step, 'out/hs2-solfec-only')

# bulk and surface materials
mat = BULK_MATERIAL (sol, model = 'KIRCHHOFF',
  young = 1E6, poisson = 0.25, density = 100)
SURFACE_MATERIAL (sol,
  model = 'SIGNORINI_COULOMB', friction = 0.1)

# template cube nodes
nodes = [0.0, 0.0, 0.0,
         0.1, 0.0, 0.0,
	 0.1, 0.1, 0.0,
	 0.0, 0.1, 0.0,
	 0.0, 0.0, 0.1,
	 0.1, 0.0, 0.1,
	 0.1, 0.1, 0.1,
	 0.0, 0.1, 0.1]

# create the array of cubes
outer = [0, M+N+M-1]
ijmap = {}
for i in range (0,M+N+M):
  for j in range (0,M+N+M):
    msh = HEX (nodes, 2, 2, 2, 0, [0, 1, 2, 3, 4, 5])
    TRANSLATE (msh, (i*(0.1+gap), j*(0.1+gap), 0))
    p1 = msh.node(0)
    p2 = msh.node(2)
    p3 = msh.node(8)
    if i in outer or j in outer:
      bod = BODY (sol, 'OBSTACLE', msh, mat)
    else:
      bod = BODY (sol, 'FINITE_ELEMENT', msh, mat)
      bod.scheme = 'DEF_LIM' # semi-implicit time integration
      bod.damping = 1E-4 # damping out free vibrations
      FIX_DIRECTION (bod, p1, (0, 0, 1))
      FIX_DIRECTION (bod, p2, (0, 0, 1))
      FIX_DIRECTION (bod, p3, (0, 0, 1))
    ijmap[(i,j)] = bod # map bodies to (i,j)-grid

# prescribe sine dwell motion 
# of the outer-most shell of bodies
for (i,j) in ijmap:
  if i in outer or j in outer:
    bod = ijmap[(i,j)]
    SET_VELOCITY (bod, bod.center, (1., 1., 0.), tsv)

# exclude contact detection between
# the outer-most shell of bodies
for (i,j) in ijmap:
  if i in outer or j in outer:
    bod0 = ijmap[(i,j)]
    bod1 = ijmap.get((i+1,j), None) if j in outer else None
    bod2 = ijmap.get((i,j+1), None) if i in outer else None
    if bod1 is not None:
      CONTACT_EXCLUDE_BODIES (bod0, bod1)
    if bod2 is not None:
      CONTACT_EXCLUDE_BODIES (bod0, bod2)

# create Newton solver
ns = NEWTON_SOLVER ()

# output interval
OUTPUT (sol, 0.03)

# run simulation
import time
start_time = time.time()
RUN (sol, ns, stop)
if RANK() == 0: print("--- %s seconds ---" % (time.time() - start_time))

# XDMF export
if sol.mode == 'READ' and not VIEWER():
  XDMF_EXPORT (sol, (0.0, stop), 'out/hs2-solfec-only/xdmf')
