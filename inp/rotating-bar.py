# rotating bar example
# testing BC and RO formulations
# ------------------------------
from scipy.io import mmread
from scipy.linalg import *

PoissonRatio = 0.26
MassDensity = 7.8E3

nodes = [-0.05, -0.05, -0.5,
          0.05, -0.05, -0.5,
          0.05,  0.05, -0.5,
         -0.05,  0.05, -0.5,
         -0.05, -0.05,  0.5,
          0.05, -0.05,  0.5,
          0.05,  0.05,  0.5,
         -0.05,  0.05,  0.5]

# here is a 2x2x10 mesh of a 0.1x0.1x0.5 rod
mesh = HEX (nodes, 2, 2, 20, 0, [0, 1, 2, 3, 4, 5])

# solver (not used)
sv = NEWTON_SOLVER ()

# compute all eigenvalues and eigenvectors
sl0 = SOLFEC ('DYNAMIC', 1E-3, 'out/rotating-bar/MK')
bl0 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = 200E4, poisson = 0.26, density = 7.8E3)
bod = BODY (sl0, 'FINITE_ELEMENT', COPY (mesh), bl0)
BODY_MM_EXPORT (bod, 'out/rotating-bar/MK/M.mtx', 'out/rotating-bar/MK/K.mtx')
M = mmread ('out/rotating-bar/MK/M.mtx').todense()
K = mmread ('out/rotating-bar/MK/K.mtx').todense()
for j in range (0, K.shape[1]):
  for i in range (j+1, K.shape[0]):
    K [j, i] = K [i, j] # above diagonal = below diagonal

x, y = eigh (K, M) # this produces y.T M y = 1 and y.T K y = x */
eval = [] # eigenvalue list
evec = [] # eigenvector list (BODY command takes a tuple (eval, evec) argument for the RO formulation)
for j in (0,1,2,3,4,5,13,18,25,33,38,55,144,146,147,150,151,153,154,155,157,167,169,171,172): #range (0, K.shape[0]):
  eval.append (x[j].real)
  for z in y[:,j]:
    evec.append (z.real)
data = (eval, evec)

# compare with MODAL_ANALYSIS results
'''
sum = 0.0
for i in range (M.shape[1]*6, M.shape[1]*12):
  a = evec [i]
  b = data [1][i]
  sum += (a-b)**2.0
print 'Evev diff norm = ', sum**0.5 # should be zero
'''

# two extreme points
p0 = (nodes[0], nodes[1], nodes[2])
p1 = (nodes[12], nodes[13], nodes[14])

# rotation: BC
h1 = 1./64.
d1 = 1.
sl1 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/BC1')
bl1 = BULK_MATERIAL (sl1, model = 'KIRCHHOFF', young = 200E4, poisson = 0.26, density = 7.8E3)
bd1 = BODY (sl1, 'FINITE_ELEMENT', COPY (mesh), bl1, form = 'BC')
bd1.scheme = 'DEF_LIM'
#bd1.damping = h1
INITIAL_VELOCITY (bd1, (0, 0, 0), (1, 0, 0))
RUN (sl1, sv, d1)

# rotation: RO
sl2 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/RO1')
bl2 = BULK_MATERIAL (sl2, model = 'KIRCHHOFF', young = 200E4, poisson = 0.26, density = 7.8E3)
bd2 = BODY (sl2, 'FINITE_ELEMENT', COPY (mesh), bl2, form = 'RO', modal = data)
bd2.scheme = 'DEF_LIM'
#bd2.damping = h1
INITIAL_VELOCITY (bd2, (0, 0, 0), (1, 0, 0))
RUN (sl2, sv, d1)

# rotation: RIG
sl3 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/RIG1')
bl3 = BULK_MATERIAL (sl3, model = 'KIRCHHOFF', young = 200E4, poisson = 0.26, density = 7.8E3)
bd3 = BODY (sl3, 'RIGID', COPY (mesh), bl3)
INITIAL_VELOCITY (bd3, (0, 0, 0), (1, 0, 0))
RUN (sl3, sv, d1)

if not VIEWER() and sl1.mode == 'READ' and sl2.mode == 'READ':
  th1 = HISTORY (sl1, [(bd1, p0, 'DX'), (bd1, p0, 'DY'), (bd1, p0, 'DZ'), (bd1, p1, 'DX'), (bd1, p1, 'DY'),(bd1, p1, 'DZ')], 0, d1)
  th2 = HISTORY (sl2, [(bd2, p0, 'DX'), (bd2, p0, 'DY'), (bd2, p0, 'DZ'), (bd2, p1, 'DX'), (bd2, p1, 'DY'),(bd2, p1, 'DZ')], 0, d1)
  lh1 = []
  lh2 = []
  for (dx0,dy0,dz0,dx1,dy1,dz1) in zip (th1[1], th1[2], th1[3], th1[4], th1[5], th1[6]):
    l = ((p0[0]+dx0-p1[0]-dx1)**2 + (p0[1]+dy0-p1[1]-dy1)**2 + (p0[2]+dz0-p1[2]-dz1)**2)**0.5
    lh1.append (l)
  for (dx0,dy0,dz0,dx1,dy1,dz1) in zip (th2[1], th2[2], th2[3], th2[4], th2[5], th2[6]):
    l = ((p0[0]+dx0-p1[0]-dx1)**2 + (p0[1]+dy0-p1[1]-dy1)**2 + (p0[2]+dz0-p1[2]-dz1)**2)**0.5
    lh2.append (l)

  try:
    import matplotlib.pyplot as plt
    plt.clf ()
    plt.title ('Bar rotation length history')
    plt.plot (th1 [0], lh1, label='BC')
    plt.plot (th2 [0], lh2, label='RO')
    plt.xlabel ('Time [s]')
    plt.ylabel ('Length [m]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/rotating-bar/length.eps')
  except ImportError:
    pass # no reaction

'''
# quasistatic test   
h1 = 1./64.
d1 = 1.
sl1 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/QS/RO')
bl1 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = 200E4, poisson = 0.26, density = 7.8E3)
bod = BODY (sl1, 'FINITE_ELEMENT', COPY (mesh), bl1, form = 'RO', modal = (eval, evec))
bod.scheme = 'DEF_LIM'
INITIAL_VELOCITY (bod, (0, 0, 0), (1, 0, 0))
#PRESSURE (bod, 0, 10)
#PRESSURE (bod, 5, 10)
bod.damping = h1# O(h) damping
RUN (sl1, sv, d1)
'''

def runtest (formulation, step, duration):
  num = len (solfec_stack)
  solfec = SOLFEC ('DYNAMIC', step, 'out/rotating-bar/' + str (num))
  solfec_stack.append (solfec)
  bod = BODY (solfec, 'FINITE_ELEMENT', COPY (mesh), bulk, form = formulation, modal = (eval, evec))
  bod.scheme = 'DEF_LIM'
  INITIAL_VELOCITY (bod, (0, 0, 0), (1, 0, 0))
  #PRESSURE (bod, 0, 10)
  #PRESSURE (bod, 5, 10)
  #bod.damping = step # O(h) damping stabilizes the convergence rate (removes oscilations on the step level => investigate)
  sv = NEWTON_SOLVER ()
  RUN (solfec, sv, duration)
  if solfec.mode == 'READ':
    SEEK (solfec, duration)
  return bod.conf

def norm(q):
  sum = 0.0
  for x in q: sum += x**2
  return sum**0.5

def diff(a, b):
  c = []
  for x, y in zip (a, b): c.append (x-y)
  return c
   
'''
h0 = 1.0 / (2.0 ** 16)
d0 = 1.0 / (2.0 ** 4)
q0 = runtest ('RO', h0, d0)

dq = []
for i in range (7, 14):
  h =  1.0 / (2.0 ** i)
  q = runtest ('RO', h, d0)
  dq.append (norm(diff(q0, q)))

print dq

if not VIEWER():
  for i in range (0, len(dq)-1):
    print dq[i] / dq[i+1]
'''
