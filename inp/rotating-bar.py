# rotating bar example
# testing BC and RO formulations
# ------------------------------
from scipy.io import mmread
from scipy.linalg import *
from math import log
import numpy as np

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

# here is a 2x2x20 mesh of a 0.1x0.1x0.5 rod
mesh = HEX (nodes, 2, 2, 20, 0, [0, 1, 2, 3, 4, 5])

# two extreme points
p0 = mesh.node (4)
p1 = mesh.node (3*3*21-5)

# solver (not used)
sv = NEWTON_SOLVER ()

# compute all eigenvalues and eigenvectors
sl0 = SOLFEC ('DYNAMIC', 1E-3, 'out/rotating-bar/MK')
bl0 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = 200E4, poisson = 0.26, density = 7.8E3)
bod = BODY (sl0, 'FINITE_ELEMENT', COPY (mesh), bl0)
eval = [] # selected eigenvalue list
evec = [] # selected eigenvector list (BODY command takes a tuple (eval, evec) argument for the RO formulation)
vsel = (0,1,2,3,4,5,13,18,25,33,38)
if 0:
  BODY_MM_EXPORT (bod, 'out/rotating-bar/MK/M.mtx', 'out/rotating-bar/MK/K.mtx')
  M = mmread ('out/rotating-bar/MK/M.mtx').todense()
  K = mmread ('out/rotating-bar/MK/K.mtx').todense()
  for j in range (0, K.shape[1]):
    for i in range (j+1, K.shape[0]):
      K [j, i] = K [i, j] # above diagonal = below diagonal
  x, y = eigh (K, M) # this produces y.T M y = 1 and y.T K y = x */
  for j in vsel:
    eval.append (x[j].real)
    for z in y[:,j]:
      evec.append (z.real)
else:
  data0 = MODAL_ANALYSIS (bod, 45, 'out/rotating-bar/MK/modal.data', verbose = 'ON')
  ndofs = mesh.nnod * 3
  for j in vsel:
    eval.append (data0[0][j])
    for k in range (j*ndofs,(j+1)*ndofs):
      evec.append (data0[1][k])
data = (eval, evec)

# undamped rotatin run
def undamped_rotation (h1, d1, TL, RG):

  # rotation: TL
  if TL:
    sl0 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/TL1_' + str(long(1./h1)) + '_' + str(long(d1)))
    bl0 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = 200E4, poisson = 0.26, density = 7.8E3)
    bd0 = BODY (sl0, 'FINITE_ELEMENT', COPY (mesh), bl0, form = 'TL')
    bd0.scheme = 'DEF_LIM'
    INITIAL_VELOCITY (bd0, (0, 0, 0), (1, 0, 0))
    RUN (sl0, sv, d1)

  # rotation: BC
  sl1 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/BC1_' + str(long(1./h1)) + '_' + str(long(d1)))
  bl1 = BULK_MATERIAL (sl1, model = 'KIRCHHOFF', young = 200E4, poisson = 0.26, density = 7.8E3)
  bd1 = BODY (sl1, 'FINITE_ELEMENT', COPY (mesh), bl1, form = 'BC')
  bd1.scheme = 'DEF_LIM'
  INITIAL_VELOCITY (bd1, (0, 0, 0), (1, 0, 0))
  RUN (sl1, sv, d1)

  # rotation: RO
  sl2 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/RO1_' + str(long(1./h1)) + '_' + str(long(d1)))
  bl2 = BULK_MATERIAL (sl2, model = 'KIRCHHOFF', young = 200E4, poisson = 0.26, density = 7.8E3)
  bd2 = BODY (sl2, 'FINITE_ELEMENT', COPY (mesh), bl2, form = 'RO', modal = data)
  bd2.scheme = 'DEF_LIM'
  INITIAL_VELOCITY (bd2, (0, 0, 0), (1, 0, 0))
  RUN (sl2, sv, d1)

  # rotation: RG
  if RG:
    sl3 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/RG1_' + str(long(1./h1)) + '_' + str(long(d1)))
    bl3 = BULK_MATERIAL (sl3, model = 'KIRCHHOFF', young = 200E4, poisson = 0.26, density = 7.8E3)
    bd3 = BODY (sl3, 'RIGID', COPY (mesh), bl3)
    INITIAL_VELOCITY (bd3, (0, 0, 0), (1, 0, 0))
    RUN (sl3, sv, d1)

  if not VIEWER() and sl1.mode == 'READ' and sl2.mode == 'READ':
    if TL:
      th0 = HISTORY (sl0, [(bd0, p0, 'DX'), (bd0, p0, 'DY'), (bd0, p0, 'DZ'), (bd0, p1, 'DX'), (bd0, p1, 'DY'),(bd0, p1, 'DZ'), (bd0, 'KINETIC'), (bd0, 'INTERNAL'), 'TIMINT'], 0, d1)
    th1 = HISTORY (sl1, [(bd1, p0, 'DX'), (bd1, p0, 'DY'), (bd1, p0, 'DZ'), (bd1, p1, 'DX'), (bd1, p1, 'DY'),(bd1, p1, 'DZ'), (bd1, 'KINETIC'), (bd1, 'INTERNAL'), 'TIMINT'], 0, d1)
    th2 = HISTORY (sl2, [(bd2, p0, 'DX'), (bd2, p0, 'DY'), (bd2, p0, 'DZ'), (bd2, p1, 'DX'), (bd2, p1, 'DY'),(bd2, p1, 'DZ'), (bd2, 'KINETIC'), (bd2, 'INTERNAL'), 'TIMINT'], 0, d1)
    if RG:
      th3 = HISTORY (sl3, ['TIMINT'], 0, d1)
    lh0 = []
    lh1 = []
    lh2 = []
    if TL:
      for (dx0,dy0,dz0,dx1,dy1,dz1) in zip (th0[1], th0[2], th0[3], th0[4], th0[5], th0[6]):
	l = ((p0[0]+dx0-p1[0]-dx1)**2 + (p0[1]+dy0-p1[1]-dy1)**2 + (p0[2]+dz0-p1[2]-dz1)**2)**0.5
	lh0.append (l)
    for (dx0,dy0,dz0,dx1,dy1,dz1) in zip (th1[1], th1[2], th1[3], th1[4], th1[5], th1[6]):
      l = ((p0[0]+dx0-p1[0]-dx1)**2 + (p0[1]+dy0-p1[1]-dy1)**2 + (p0[2]+dz0-p1[2]-dz1)**2)**0.5
      lh1.append (l)
    for (dx0,dy0,dz0,dx1,dy1,dz1) in zip (th2[1], th2[2], th2[3], th2[4], th2[5], th2[6]):
      l = ((p0[0]+dx0-p1[0]-dx1)**2 + (p0[1]+dy0-p1[1]-dy1)**2 + (p0[2]+dz0-p1[2]-dz1)**2)**0.5
      lh2.append (l)

    tot0 = []
    tot1 = []
    tot2 = []
    if TL:
      for (kin, int) in zip (th0[7], th0[8]):
	tot0.append (kin+int)
    for (kin, int) in zip (th1[7], th1[8]):
      tot1.append (kin+int)
    for (kin, int) in zip (th2[7], th2[8]):
      tot2.append (kin+int)

    tt0 = 0.0
    tt1 = 0.0
    tt2 = 0.0
    tt3 = 0.0
    runts = []
    runls = []
    colors = []
    hatchs = []
    if TL:
      for tt in th0[9]:
	tt0 += tt
      runts.append (tt0)
      runls.append ('TL')
      colors.append ('r')
      hatchs.append ('//')
    for tt in th1[9]:
      tt1 += tt
    runts.append (tt1)
    runls.append ('BC')
    colors.append ('g')
    hatchs.append ('\\\\')
    for tt in th2[9]:
      tt2 += tt
    runts.append (tt2)
    runls.append ('RO')
    colors.append ('b')
    hatchs.append ('o')
    if RG:
      for tt in th3[1]:
	tt3 += tt
      runts.append (tt3)
      runls.append ('RG')
      colors.append ('y')
      hatchs.append ('.')

    try:
      import matplotlib.pyplot as plt

      hstr = ' [h = 1/' + str(long(1/h1)) + ']'
      hdstr = ' [h = 1/' + str(long(1/h1)) + ', ' + str(long(d1)) + ']'
      plt.clf ()
      plt.title ('Rotating bar: length history' + hstr)
      if TL: plt.plot (th0 [0], lh0, label='TL', marker = 's')
      plt.plot (th1 [0], lh1, label='BC')
      plt.plot (th2 [0], lh2, label='RO', ls = '--', marker = 'o')
      plt.xlabel ('Time [s]')
      plt.ylabel ('Length [m]')
      plt.legend(loc = 'upper right')
      plt.savefig ('out/rotating-bar/undamp_length' + str(long(1/h1)) + '_' + str(long(d1)) + '.png')

      plt.clf ()
      plt.title ('Rotating bar: total energy' + hstr)
      if TL: plt.plot (th0 [0], tot0, label='TL', ls = '-.')
      plt.plot (th1 [0], tot1, label='BC')
      plt.plot (th2 [0], tot2, label='RO', ls = '--')
      plt.xlabel ('Time [s]')
      plt.ylabel ('Energy [J]')
      plt.legend(loc = 'upper right')
      plt.savefig ('out/rotating-bar/undamp_energy' + str(long(1/h1)) + '_' + str(long(d1)) + '.png')

      plt.clf ()
      plt.title ('Rotating bar: runtimes' + hdstr)
      for (p, r, c, t) in zip (np.arange (len(runts)), runts, colors, hatchs):
        plt.bar (p, r, 0.7, color = c, hatch = t)
      plt.xlabel ('Formulation')
      plt.ylabel ('Runtime [s]')
      plt.xticks (np.arange(len(runts))+0.35, runls)
      plt.savefig ('out/rotating-bar/undamp_runtimes' + str(long(1/h1)) + '_' + str(long(d1)) + '.png')

    except ImportError:
      pass # no reaction

undamped_rotation (1./64., 1, 1, 0)
undamped_rotation (1./256., 1, 1, 0)
undamped_rotation (1./64., 10, 1, 0)
undamped_rotation (1./64., 100, 1, 1)

# manual bar plot of output size
try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.title ('Rotating bar: output sizes [h=1/64, 100s]')
  plt.bar ([0], [572], 0.7, color = 'r', hatch = '//')
  plt.bar ([1], [572], 0.7, color = 'b', hatch = '\\\\')
  plt.bar ([2], [33], 0.7, color = 'g', hatch = 'o')
  plt.bar ([3], [23], 0.7, color = 'y', hatch = '.')
  plt.xlabel ('Formulation')
  plt.ylabel ('Output size [MB]')
  plt.xticks (np.arange(4)+0.35, ['TL','BC','RO','RG'])
  plt.savefig ('out/rotating-bar/undamp_outsizes_64_100.png')

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

'''
def runtest (formulation, step, duration):
  solfec = SOLFEC ('DYNAMIC', step, 'out/rotating-bar/CR/' + formulation + str(step))
  bulk = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 200E4, poisson = 0.26, density = 7.8E3)
  bod = BODY (solfec, 'FINITE_ELEMENT', COPY (mesh), bulk, form = formulation, modal = data)
  bod.scheme = 'DEF_LIM'
  INITIAL_VELOCITY (bod, (0, 0, 0), (1, 0, 0))
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
