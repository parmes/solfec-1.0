# 16 bars array example
# testing BC and RO formulations (W assembling)
# ----------------------------------------------
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

# base middle point
p0 = MASS_CENTER (mesh)

# number of bars
nw = 4

# obstacle mesh
obsm = HEX (nodes, 1, 1, 1, 0, [0, 1, 2, 3, 4, 5])
SCALE (obsm, (nw, nw, 0.2))
TRANSLATE (obsm, (0, 0, -0.6))

# solver
sv = GAUSS_SEIDEL_SOLVER (1E-3, 1000)

# impact comparison
def create_bars (h1, E, frict, damp, formulation):

  # compute all eigenvalues and eigenvectors
  if formulation == 'RO':
    pt0 = 'out/16-bars/MK_%g_%g_%g_%g'%(h1, E, frict, damp)
    sl0 = SOLFEC ('DYNAMIC', 1E-3, pt0)
    bl0 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = E, poisson = PoissonRatio, density = MassDensity)
    bod = BODY (sl0, 'FINITE_ELEMENT', COPY (mesh), bl0)
    eval = [] # selected eigenvalue list
    evec = [] # selected eigenvector list (BODY command takes a tuple (eval, evec) argument for the RO formulation)
    vsel = range (0, 32)

    if 0:
      BODY_MM_EXPORT (bod, pt0+'/M.mtx', pt0+'/K.mtx')
      M = mmread (pt0+'/M.mtx').todense()
      K = mmread (pt0+'/K.mtx').todense()
      for j in range (0, K.shape[1]):
	for i in range (j+1, K.shape[0]):
	  K [j, i] = K [i, j] # above diagonal = below diagonal
      x, y = eigh (K, M) # this produces y.T M y = 1 and y.T K y = x */
      for j in vsel:
	eval.append (x[j].real)
	for z in y[:,j]:
	  evec.append (z.real)
    else:
      data0 = MODAL_ANALYSIS (bod, 45, pt0 + '/modal.data', verbose = 'ON', abstol = 1E-14)
      ndofs = mesh.nnod * 3
      for j in vsel:
	eval.append (data0[0][j])
	for k in range (j*ndofs,(j+1)*ndofs):
	  evec.append (data0[1][k])
    data = (eval, evec)

  # 16 bars domain
  sl2 = SOLFEC ('DYNAMIC', h1, 'out/16-bars/%s_%g_%g_%g_%g'%(formulation, h1, E, frict, damp))
  SURFACE_MATERIAL (sl2, model = 'SIGNORINI_COULOMB', friction = frict, restitution = 0.0)
  bl2 = BULK_MATERIAL (sl2, model = 'KIRCHHOFF', young = E, poisson = PoissonRatio, density = MassDensity)
  GRAVITY (sl2, (0, 0, -9.8))
  for i in range (0, nw):
    for j in range (0, nw):
      shp = COPY (mesh)
      TRANSLATE (shp, ((1-nw)*0.05+0.1*i, (1-nw)*0.05+0.1*j, 0))
      if formulation == 'RO':
	bd2 = BODY (sl2, 'FINITE_ELEMENT', shp, bl2, form = formulation, modal = data)
	bd2.scheme = 'DEF_LIM'
	bd2.damping = damp
      elif formulation == 'BC':
	bd2 = BODY (sl2, 'FINITE_ELEMENT', shp, bl2, form = formulation)
	bd2.scheme = 'DEF_LIM'
	bd2.damping = damp
      else: bd2 = BODY (sl2, 'RIGID', shp, bl2)
  BODY (sl2, 'OBSTACLE', COPY (obsm), bl2)

  return sl2

# juxtaposed runtime bar plot
def juxtapose_plot (input_data, xdata, title, output_file):
  try:
    import matplotlib.pyplot as plt
    import matplotlib.figure as fig

    dat = [[],[],[],[]]
    for x in input_data:
      for i in range (0, 4):
	dat [i].append (x [i])

    N = len (input_data)
    TIMINT = np.array (dat [0])
    CONDET = np.array (dat [1])
    LOCDYN = np.array (dat [2])
    CONSOL = np.array (dat [3])

    ind = np.arange(N)
    width = 0.6

    p1 = plt.bar (ind, TIMINT, width, color='r')
    p2 = plt.bar (ind, CONDET, width, color='g', bottom=TIMINT, hatch='\\\\')
    p3 = plt.bar (ind, LOCDYN, width, color='y', bottom=TIMINT+CONDET)
    p4 = plt.bar (ind, CONSOL, width, color='c', bottom=TIMINT+CONDET+LOCDYN, hatch='//')
    tt = TIMINT+CONDET+LOCDYN+CONSOL
    yy = tt[0]

    for (ii, ti, cd, ld) in zip (ind, TIMINT, CONDET, LOCDYN):
      x = ii + width/2.
      y = ti + cd + ld/2.
      plt.text (x, y, '%.3g s'%(ld), ha='center', va='center')

    plt.title (title)
    plt.ylabel ('Time [s]')
    plt.xticks (ind+width/2., xdata)
    plt.legend ((p4[0], p3[0], p2[0], p1[0]), ('Contact solution', '$\mathbf{W}$ assembling', 'Contact detection', 'Time integration'), loc = "upper right")
    plt.axis (xmin = -0.2, xmax = ind[N-1]+width+0.2, ymax = 1.05 * yy)
    plt.savefig (output_file)
    plt.clf ()
  except ImportError:
    pass # no reaction

# create and run tests
h = 1. / 8192.
d = h * 10
E = 200E9
damp = h
fric = 0.1
sl0 = create_bars (h, E, fric, damp, 'BC')
sl1 = create_bars (h, E, fric, damp, 'RO')
sl2 = create_bars (h, E, fric, damp, 'RG')
RUN (sl0, sv, d)
RUN (sl1, sv, d)
RUN (sl2, sv, d)

# plot runtimes
dat0 = []
dat1 = []
dat2 = []
if not VIEWER() and sl0.mode == 'READ' and sl1.mode == 'READ' and sl2.mode == 'READ':
  th0 = HISTORY (sl0, ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL'], 0, d)
  th1 = HISTORY (sl1, ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL'], 0, d)
  th2 = HISTORY (sl2, ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL'], 0, d)

  dat0 = [0.0, 0.0, 0.0, 0.0]
  for (ti, cu, cd, ld, sv) in zip (th0[1], th0[2], th0[3], th0[4], th0[5]):
    dat0 [0] += ti
    dat0 [1] += cu + cd
    dat0 [2] += ld
    dat0 [3] += sv
  dat1 = [0.0, 0.0, 0.0, 0.0]
  for (ti, cu, cd, ld, sv) in zip (th1[1], th1[2], th1[3], th1[4], th1[5]):
    dat1 [0] += ti
    dat1 [1] += cu + cd
    dat1 [2] += ld
    dat1 [3] += sv
  dat2 = [0.0, 0.0, 0.0, 0.0]
  for (ti, cu, cd, ld, sv) in zip (th2[1], th2[2], th2[3], th2[4], th2[5]):
    dat2 [0] += ti
    dat2 [1] += cu + cd
    dat2 [2] += ld
    dat2 [3] += sv

  # make plot
  juxtapose_plot ([dat0, dat1, dat2], ('BC', 'RO', 'RG'), '16 bars, 10 steps runtimes', 'out/16-bars/16_bars_10_steps.eps')
