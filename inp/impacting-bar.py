# impacting bar example
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
p1 = mesh.node (3*3*21-5)

# obstacle mesh
obsm = HEX (nodes, 1, 1, 1, 0, [0, 1, 2, 3, 4, 5])
SCALE (obsm, (2, 2, 0.2))
TRANSLATE (obsm, (0, 0, -0.6))

# solver (not used)
sv = NEWTON_SOLVER ()

# compute all eigenvalues and eigenvectors
sl0 = SOLFEC ('DYNAMIC', 1E-3, 'out/impacting-bar/MK')
bl0 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
bod = BODY (sl0, 'FINITE_ELEMENT', COPY (mesh), bl0)
eval = [] # selected eigenvalue list
evec = [] # selected eigenvector list (BODY command takes a tuple (eval, evec) argument for the RO formulation)
vsel = (0,1,2,3,4,5,13,18,25,33,38)
if 0:
  BODY_MM_EXPORT (bod, 'out/impacting-bar/MK/M.mtx', 'out/impacting-bar/MK/K.mtx')
  M = mmread ('out/impacting-bar/MK/M.mtx').todense()
  K = mmread ('out/impacting-bar/MK/K.mtx').todense()
  for j in range (0, K.shape[1]):
    for i in range (j+1, K.shape[0]):
      K [j, i] = K [i, j] # above diagonal = below diagonal
  x, y = eigh (K, M) # this produces y.T M y = 1 and y.T K y = x */
  for j in vsel:
    eval.append (x[j].real)
    for z in y[:,j]:
      evec.append (z.real)
else:
  data0 = MODAL_ANALYSIS (bod, 45, 'out/impacting-bar/MK/modal.data', verbose = 'ON')
  ndofs = mesh.nnod * 3
  for j in vsel:
    eval.append (data0[0][j])
    for k in range (j*ndofs,(j+1)*ndofs):
      evec.append (data0[1][k])
data = (eval, evec)

# impact comparison
def impact_comparison (h1, d1, E, pow0, pow1):

  toplot = []
  for p in range (pow0, pow1):

    damping = 1.0 / (2.0**p)

    # rotation: BC
    sl1 = SOLFEC ('DYNAMIC', h1, 'out/impacting-bar/BC_%g_%g_%g_%g'%(h1, d1, E, damping))
    bl1 = BULK_MATERIAL (sl1, model = 'KIRCHHOFF', young = E, poisson = PoissonRatio, density = MassDensity)
    bd1 = BODY (sl1, 'FINITE_ELEMENT', COPY (mesh), bl1, form = 'BC')
    bd1.scheme = 'DEF_LIM'
    bd1.damping = damping
    INITIAL_VELOCITY (bd1, (0, 0, -1), (0, 0, 0))
    BODY (sl1, 'OBSTACLE', COPY (obsm), bl1)
    RUN (sl1, sv, d1)
    if not VIEWER() and sl1.mode == 'READ':
      th1 = HISTORY (sl1, [(bd1, p1, 'DZ')], 0, d1)

    # rotation: RO
    sl2 = SOLFEC ('DYNAMIC', h1, 'out/impacting-bar/RO_%g_%g_%g_%g'%(h1, d1, E, damping))
    bl2 = BULK_MATERIAL (sl2, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
    bd2 = BODY (sl2, 'FINITE_ELEMENT', COPY (mesh), bl2, form = 'RO', modal = data)
    bd2.scheme = 'DEF_LIM'
    bd2.damping = damping
    INITIAL_VELOCITY (bd2, (0, 0, -1), (0, 0, 0))
    BODY (sl2, 'OBSTACLE', COPY (obsm), bl2)
    RUN (sl2, sv, d1)
    if not VIEWER() and sl2.mode == 'READ':
      th2 = HISTORY (sl2, [(bd2, p1, 'DZ')], 0, d1)
    
    if not VIEWER() and sl1.mode == 'READ' and sl2.mode == 'READ':
      toplot.append (('$\eta=1/' + str(2**p) + '$', th1[0], th1[1], th2[1]))

  if not VIEWER() and sl1.mode == 'READ' and sl2.mode == 'READ':

    try:
      import matplotlib.pyplot as plt

      hstr = ' ($E = %g, h = %g'%(E, h1) + '$)'

      plt.clf ()
      plt.title ('Rotating bar: BC top point DZ' + hstr)
      for dat in toplot:
        plt.plot (dat[1], dat[2], label=dat[0])
      plt.xlabel ('Time [s]')
      plt.ylabel ('DZ [m]')
      plt.semilogy (10)
      plt.legend(loc = 'best')
      plt.savefig ('out/impacting-bar/ib_dz_BC_h%g_d%g_E%g.eps'%(h1,d1,E))

      plt.clf ()
      plt.title ('Rotating bar: RO top point DZ' + hstr)
      for dat in toplot:
        plt.plot (dat[1], dat[3], label=dat[0])
      plt.xlabel ('Time [s]')
      plt.ylabel ('DZ [m]')
      plt.semilogy (10)
      plt.legend(loc = 'best')
      plt.savefig ('out/impacting-bar/ib_dz_RO_h%g_d%g_E%g.eps'%(h1,d1,E))

    except ImportError:
      pass # no reaction

# impacting testes
impact_comparison (1/256., 1/4., 200E9, 8, 15)
