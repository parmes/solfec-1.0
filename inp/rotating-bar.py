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
bl0 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
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
  data0 = MODAL_ANALYSIS (bod, 45, 'out/rotating-bar/MK/modal.data', verbose = 'ON', abstol = 1E-13)
  ndofs = mesh.nnod * 3
  for j in vsel:
    eval.append (data0[0][j])
    for k in range (j*ndofs,(j+1)*ndofs):
      evec.append (data0[1][k])
data = (eval, evec)

# undamped rotatin run
def rotation_test (h1, d1, damping, TL, LEN, ENE):

  # rotation: TL
  if TL:
    sl0 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/TL1_%g_'%damping + str(long(1./h1)) + '_' + str(long(d1)))
    bl0 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
    bd0 = BODY (sl0, 'FINITE_ELEMENT', COPY (mesh), bl0, form = 'TL')
    bd0.scheme = 'DEF_LIM'
    bd0.damping = damping
    INITIAL_VELOCITY (bd0, (0, 0, 0), (1, 0, 0))
    RUN (sl0, sv, d1)

  # rotation: BC
  sl1 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/BC1_%g_'%damping + str(long(1./h1)) + '_' + str(long(d1)))
  bl1 = BULK_MATERIAL (sl1, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
  bd1 = BODY (sl1, 'FINITE_ELEMENT', COPY (mesh), bl1, form = 'BC')
  bd1.scheme = 'DEF_LIM'
  bd1.damping = damping
  INITIAL_VELOCITY (bd1, (0, 0, 0), (1, 0, 0))
  RUN (sl1, sv, d1)

  # rotation: RO
  sl2 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/RO1_%g_'%damping + str(long(1./h1)) + '_' + str(long(d1)))
  bl2 = BULK_MATERIAL (sl2, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
  bd2 = BODY (sl2, 'FINITE_ELEMENT', COPY (mesh), bl2, form = 'RO', modal = data)
  bd2.scheme = 'DEF_LIM'
  bd2.damping = damping
  INITIAL_VELOCITY (bd2, (0, 0, 0), (1, 0, 0))
  RUN (sl2, sv, d1)

  if not VIEWER() and sl1.mode == 'READ' and sl2.mode == 'READ':
    if TL:
      th0 = HISTORY (sl0, [(bd0, p0, 'DX'), (bd0, p0, 'DY'), (bd0, p0, 'DZ'), (bd0, p1, 'DX'), (bd0, p1, 'DY'),(bd0, p1, 'DZ'), (bd0, 'KINETIC'), (bd0, 'INTERNAL'), 'TIMINT'], 0, d1)
    th1 = HISTORY (sl1, [(bd1, p0, 'DX'), (bd1, p0, 'DY'), (bd1, p0, 'DZ'), (bd1, p1, 'DX'), (bd1, p1, 'DY'),(bd1, p1, 'DZ'), (bd1, 'KINETIC'), (bd1, 'INTERNAL'), 'TIMINT'], 0, d1)
    th2 = HISTORY (sl2, [(bd2, p0, 'DX'), (bd2, p0, 'DY'), (bd2, p0, 'DZ'), (bd2, p1, 'DX'), (bd2, p1, 'DY'),(bd2, p1, 'DZ'), (bd2, 'KINETIC'), (bd2, 'INTERNAL'), 'TIMINT'], 0, d1)
    lh0 = []
    lh1 = []
    lh2 = []
    if TL:
      for (dx0,dy0,dz0,dx1,dy1,dz1) in zip (th0[1], th0[2], th0[3], th0[4], th0[5], th0[6]):
	l = ((p0[0]+dx0-p1[0]-dx1)**2 + (p0[1]+dy0-p1[1]-dy1)**2 + (p0[2]+dz0-p1[2]-dz1)**2)**0.5
	lh0.append (l-1.0)
    for (dx0,dy0,dz0,dx1,dy1,dz1) in zip (th1[1], th1[2], th1[3], th1[4], th1[5], th1[6]):
      l = ((p0[0]+dx0-p1[0]-dx1)**2 + (p0[1]+dy0-p1[1]-dy1)**2 + (p0[2]+dz0-p1[2]-dz1)**2)**0.5
      lh1.append (l-1.0)
    for (dx0,dy0,dz0,dx1,dy1,dz1) in zip (th2[1], th2[2], th2[3], th2[4], th2[5], th2[6]):
      l = ((p0[0]+dx0-p1[0]-dx1)**2 + (p0[1]+dy0-p1[1]-dy1)**2 + (p0[2]+dz0-p1[2]-dz1)**2)**0.5
      lh2.append (l-1.0)

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

    try:
      import matplotlib.pyplot as plt

      if d1 > 10.: format = '.png'
      else: format = '.eps'

      if damping > 0.0: head = 'damped_%g'%damping
      else: head = 'undamp'

      hstr = ' ($h = 1/' + str(long(1/h1)) + ', \eta = %g'%damping + '$)'

      if LEN:
	plt.clf ()
	plt.title ('Rotating bar: elongation ' + hstr)
	if TL: plt.plot (th0 [0], lh0, label='TL', marker = 's')
	plt.plot (th1 [0], lh1, label='BC')
	plt.plot (th2 [0], lh2, label='RO', ls = '--', marker = 'o')
	plt.xlabel ('Time [s]')
	plt.ylabel ('Elongation [m]')
	plt.legend(loc = 'upper right')
	plt.gcf().subplots_adjust(left=0.15)
	plt.savefig ('out/rotating-bar/rb_' + head + '_length' + str(long(1/h1)) + '_' + str(long(d1)) + format)

      if ENE:
	plt.clf ()
	plt.title ('Rotating bar: total energy' + hstr)
	if TL: plt.plot (th0 [0], tot0, label='TL', ls = '-.')
	plt.plot (th1 [0], tot1, label='BC')
	plt.plot (th2 [0], tot2, label='RO', ls = '--')
	plt.xlabel ('Time [s]')
	plt.ylabel ('Energy [J]')
	plt.legend(loc = 'upper right')
	plt.savefig ('out/rotating-bar/rb_' + head + '_energy' + str(long(1/h1)) + '_' + str(long(d1)) + format)

    except ImportError:
      pass # no reaction

# rotatin comparison
def rotation_comparison (n, pow0, pow1, d1, E, damping):

  # nxnx10x mesh
  mesh = HEX (nodes, n, n, 10*n, 0, [0, 1, 2, 3, 4, 5])

  # top point
  p1 = mesh.node ((n+1)*(n+1)*(10*n+1)-((n+1)*(n+1)/2+1))

  # modal analysis
  path = 'out/rotating-bar/modal' + str(n)
  sl0 = SOLFEC ('DYNAMIC', 1.0, path)
  bl0 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = E, poisson = PoissonRatio, density = MassDensity)
  bd0 = BODY (sl0, 'FINITE_ELEMENT', COPY (mesh), bl0, form = 'TL')
  data0 = MODAL_ANALYSIS (bd0, 45, path + '/data', verbose = 'ON', abstol = 1E-13)
  ndofs = mesh.nnod * 3
  eval = []
  evec = []
  for j in vsel:
    eval.append (data0[0][j])
    for k in range (j*ndofs,(j+1)*ndofs):
      evec.append (data0[1][k])
  data = (eval, evec)

  toplot = []
  for p in range (pow0, pow1):

    h1 = 1.0 / (2.**p)

    # rotation: RIGID
    sl0 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/RG3_%d_%g_%g_'%(n, E, damping) + str(long(1./h1)) + '_' + str(long(d1)))
    bl0 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = E, poisson = PoissonRatio, density = MassDensity)
    bd0 = BODY (sl0, 'RIGID', COPY (mesh), bl0)
    INITIAL_VELOCITY (bd0, (0, 0, 0), (1, 0, 0))
    RUN (sl0, sv, d1)
    if not VIEWER() and sl0.mode == 'READ':
      th0 = HISTORY (sl0, [(bd0, p1, 'DZ')], 0, d1)

    # rotation: BC
    sl1 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/BC3_%d_%g_%g_'%(n, E, damping) + str(long(1./h1)) + '_' + str(long(d1)))
    bl1 = BULK_MATERIAL (sl1, model = 'KIRCHHOFF', young = E, poisson = PoissonRatio, density = MassDensity)
    bd1 = BODY (sl1, 'FINITE_ELEMENT', COPY (mesh), bl1, form = 'BC')
    bd1.scheme = 'DEF_LIM'
    bd1.damping = damping
    INITIAL_VELOCITY (bd1, (0, 0, 0), (1, 0, 0))
    RUN (sl1, sv, d1)
    if not VIEWER() and sl1.mode == 'READ':
      th1 = HISTORY (sl1, [(bd1, p1, 'DZ')], 0, d1)

    # rotation: RO
    sl2 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/RO3_%d_%g_%g_'%(n, E, damping) + str(long(1./h1)) + '_' + str(long(d1)))
    bl2 = BULK_MATERIAL (sl2, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
    bd2 = BODY (sl2, 'FINITE_ELEMENT', COPY (mesh), bl2, form = 'RO', modal = data)
    bd2.scheme = 'DEF_LIM'
    bd2.damping = damping
    INITIAL_VELOCITY (bd2, (0, 0, 0), (1, 0, 0))
    RUN (sl2, sv, d1)
    if not VIEWER() and sl2.mode == 'READ':
      th2 = HISTORY (sl2, [(bd2, p1, 'DZ')], 0, d1)
    
    if not VIEWER() and sl0.mode == 'READ' and sl1.mode == 'READ' and sl2.mode == 'READ':
      toplot.append (('h=1/' + str(2**p), th0[0], th0[1], th1[1], th2[1]))

  if not VIEWER() and sl0.mode == 'READ' and sl1.mode == 'READ' and sl2.mode == 'READ':

    DBC = []
    DRO = []
    for dat in toplot:
      vro = []
      vbc = []
      for (rg, bc, ro) in zip (dat[2], dat[3], dat[4]):
	vbc.append (abs(bc-rg))
	vro.append (abs(ro-rg))
      DBC.append ((dat[1], vbc, dat[0]))
      DRO.append ((dat[1], vro, dat[0]))

    try:
      import matplotlib.pyplot as plt

      hstr = ' ($E = %g, \eta = %g'%(E, damping) + '$)'

      plt.clf ()
      plt.title ('Rotating bar: BC top point $\|dz_{BC} - dz_{RIGID}\|$' + hstr)
      for dat in DBC:
        plt.plot (dat[0], dat[1], label=dat[2])
      plt.xlabel ('Time [s]')
      plt.ylabel ('$\|dz_{BC} - dz_{RIGID}\|$ [m]')
      plt.semilogy (10)
      plt.legend(loc = 'best')
      plt.savefig ('out/rotating-bar/rb_dz_BC_n%d_E%g_'%(n,E) + str(long(d1)) + '.eps')

      plt.clf ()
      plt.title ('Rotating bar: RO top point $\|dz_{RO} - dz_{RIGID}\|$' + hstr)
      for dat in DRO:
        plt.plot (dat[0], dat[1], label=dat[2])
      plt.xlabel ('Time [s]')
      plt.ylabel ('$\|dz_{RO} - dz_{RIGID}\|$ [m]')
      plt.semilogy (10)
      plt.legend(loc = 'best')
      plt.savefig ('out/rotating-bar/rb_dz_RO_n%d_E%g_'%(n,E) + str(long(d1)) + '.eps')

    except ImportError:
      pass # no reaction
 
def barlabel (plt, rects):
  for rect in rects:
    height = rect.get_height()
    width = rect.get_width()
    plt.text(rect.get_x()+width/2., height+0.001*width, '%.3g'%height, ha='center', va='bottom')

def undamped_rotation_runtimes (h1, d1):

  # 4x4x40 mesh
  mesh = HEX (nodes, 4, 4, 40, 0, [0, 1, 2, 3, 4, 5])

  # rotation: TL
  TLpath = 'out/rotating-bar/TL2_' + str(long(1./h1)) + '_' + str(long(d1))
  sl0 = SOLFEC ('DYNAMIC', h1, TLpath)
  bl0 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
  bd0 = BODY (sl0, 'FINITE_ELEMENT', COPY (mesh), bl0, form = 'TL')
  bd0.scheme = 'DEF_LIM'
  INITIAL_VELOCITY (bd0, (0, 0, 0), (1, 0, 0))
  RUN (sl0, sv, d1)

  # modal analysis
  data0 = MODAL_ANALYSIS (bd0, 45, TLpath + '/modal.data', verbose = 'ON', abstol = 1E-13)
  ndofs = mesh.nnod * 3
  eval = []
  evec = []
  for j in vsel:
    eval.append (data0[0][j])
    for k in range (j*ndofs,(j+1)*ndofs):
      evec.append (data0[1][k])
  data = (eval, evec)

  # rotation: BC
  sl1 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/BC2_' + str(long(1./h1)) + '_' + str(long(d1)))
  bl1 = BULK_MATERIAL (sl1, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
  bd1 = BODY (sl1, 'FINITE_ELEMENT', COPY (mesh), bl1, form = 'BC')
  bd1.scheme = 'DEF_LIM'
  INITIAL_VELOCITY (bd1, (0, 0, 0), (1, 0, 0))
  RUN (sl1, sv, d1)

  # rotation: RO
  sl2 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/RO2_' + str(long(1./h1)) + '_' + str(long(d1)))
  bl2 = BULK_MATERIAL (sl2, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
  bd2 = BODY (sl2, 'FINITE_ELEMENT', COPY (mesh), bl2, form = 'RO', modal = data)
  bd2.scheme = 'DEF_LIM'
  INITIAL_VELOCITY (bd2, (0, 0, 0), (1, 0, 0))
  RUN (sl2, sv, d1)

  # rotation: RG
  sl3 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/RG2_' + str(long(1./h1)) + '_' + str(long(d1)))
  bl3 = BULK_MATERIAL (sl3, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
  bd3 = BODY (sl3, 'RIGID', COPY (mesh), bl3)
  INITIAL_VELOCITY (bd3, (0, 0, 0), (1, 0, 0))
  RUN (sl3, sv, d1)

  # numerb of steps
  nstp = d1 / h1

  # bar plot comparison on 4x4x40 mesh
  if not VIEWER() and sl1.mode == 'READ' and sl2.mode == 'READ':
    th0 = HISTORY (sl0, ['TIMINT'], 0, d1)
    th1 = HISTORY (sl1, ['TIMINT'], 0, d1)
    th2 = HISTORY (sl2, ['TIMINT'], 0, d1)
    th3 = HISTORY (sl3, ['TIMINT'], 0, d1)

    tt0 = 0.0
    tt1 = 0.0
    tt2 = 0.0
    tt3 = 0.0
    runts = []
    runls = []
    colors = []
    hatchs = []
    for tt in th0[1]:
      tt0 += tt
    runts.append (tt0/nstp)
    runls.append ('TL')
    colors.append ('g')
    hatchs.append ('.')
    for tt in th1[1]:
      tt1 += tt
    runts.append (tt1/nstp)
    runls.append ('BC')
    colors.append ('r')
    hatchs.append ('//')
    for tt in th2[1]:
      tt2 += tt
    runts.append (tt2/nstp)
    runls.append ('RO')
    colors.append ('y')
    hatchs.append ('\\\\')
    for tt in th3[1]:
      tt3 += tt
    runts.append (tt3/nstp)
    runls.append ('RG')
    colors.append ('b')
    hatchs.append ('o')

    try:
      import matplotlib.pyplot as plt

      plt.clf ()
      plt.title ('Rotating bar: avg. runtime / time step (4x4x40 mesh)')
      for (p, r, c, t) in zip (np.arange (len(runts)), runts, colors, hatchs):
        r = plt.bar (p, r, 0.7, color = c, hatch = t)
        barlabel (plt, r)
      plt.xlabel ('Formulation')
      plt.ylabel ('Runtime [s]')
      plt.xticks (np.arange(len(runts))+0.35, runls)
      plt.savefig ('out/rotating-bar/rb_undamp_avgrun_4x4x40.eps')

    except ImportError:
      pass # no reaction

  # scaling runtimes
  scaling_data = []
  for n in range (2, 10):

    # nxnx10x mesh
    mesh = HEX (nodes, n, n, 10*n, 0, [0, 1, 2, 3, 4, 5])

    # modal analysis
    path = 'out/rotating-bar/modal' + str(n)
    sl0 = SOLFEC ('DYNAMIC', h1, path)
    bl0 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
    bd0 = BODY (sl0, 'FINITE_ELEMENT', COPY (mesh), bl0, form = 'TL')
    data0 = MODAL_ANALYSIS (bd0, 45, path + '/data', verbose = 'ON', abstol = 1E-13)
    ndofs = mesh.nnod * 3
    eval = []
    evec = []
    for j in vsel:
      eval.append (data0[0][j])
      for k in range (j*ndofs,(j+1)*ndofs):
	evec.append (data0[1][k])
    data = (eval, evec)

    # rotation: BC
    sl1 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/BC2_' + str(long(1./h1)) + '_' + str(long(d1)) + '_' + str(n))
    bl1 = BULK_MATERIAL (sl1, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
    bd1 = BODY (sl1, 'FINITE_ELEMENT', COPY (mesh), bl1, form = 'BC')
    bd1.scheme = 'DEF_LIM'
    INITIAL_VELOCITY (bd1, (0, 0, 0), (1, 0, 0))
    RUN (sl1, sv, d1)

    # rotation: RO
    sl2 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/RO2_' + str(long(1./h1)) + '_' + str(long(d1)) + '_' + str(n))
    bl2 = BULK_MATERIAL (sl2, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
    bd2 = BODY (sl2, 'FINITE_ELEMENT', COPY (mesh), bl2, form = 'RO', modal = data)
    bd2.scheme = 'DEF_LIM'
    INITIAL_VELOCITY (bd2, (0, 0, 0), (1, 0, 0))
    RUN (sl2, sv, d1)

    # append with a tuple
    scaling_data.append ((str(n)+'x'+str(n)+'x'+str(10*n), sl1, sl2))

  if not VIEWER() and scaling_data [0][1].mode == 'READ':

    LABELS = []
    BCVALS = []
    ROVALS = []
    for data in scaling_data:
      th0 = HISTORY (data[1], ['TIMINT'], 0, d1)
      th1 = HISTORY (data[2], ['TIMINT'], 0, d1)

      tt0 = 0.0
      tt1 = 0.0
      for tt in th0[1]:
	tt0 += tt
      for tt in th1[1]:
	tt1 += tt

      LABELS.append (data[0])
      BCVALS.append (tt0/nstp)
      ROVALS.append (tt1/nstp)

    try:
      import matplotlib.pyplot as plt

      ind = np.arange (len(LABELS))
      wdt = 0.35
      plt.clf ()
      plt.title ('Rotating bar: scaling of avg. runtime / time step')
      plt.bar (ind, BCVALS, wdt, color = 'r', hatch = '//', label = 'BC')
      plt.bar (ind+wdt, ROVALS, wdt, color = 'y', hatch = '\\\\', label = 'RO')
      plt.xlabel ('Mesh size')
      plt.ylabel ('Runtime [s]')
      plt.xticks (ind+wdt, LABELS)
      plt.legend(loc = 'upper left')
      plt.savefig ('out/rotating-bar/rb_undamp_avgsca.eps')

    except ImportError:
      pass # no reaction

# single convergence rate test run
def runtest (formulation, E, damping, duration, step):
  solfec = SOLFEC ('DYNAMIC', step, 'out/rotating-bar/cnv_%s_%g_%g_%g_%g'%(formulation, E, damping, duration, step))
  bulk = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = E, poisson = PoissonRatio, density = MassDensity)
  bod = BODY (solfec, 'FINITE_ELEMENT', COPY (mesh), bulk, form = formulation, modal = data)
  bod.scheme = 'DEF_LIM'
  bod.damping = damping
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

# convergence rate tests
def convtest (formulation, E, damping, pow0, pow1, pow2, pow3):
   
  h0 = 1.0 / (2.0 ** pow0)
  d0 = 1.0 / (2.0 ** pow3)
  q0 = runtest (formulation, E, damping, d0, h0)

  dq = []
  hh = []
  for i in range (pow1, pow2):
    h =  1.0 / (2.0 ** i)
    q = runtest (formulation, E, damping, d0, h)
    dq.append (norm(diff(q0, q)))
    hh.append (h)

  print formulation, '=> convergence ratios:'
  for i in range (0, len(dq)-1):
    print dq[i] / dq[i+1]

  return (hh, dq)


# run conv. tests
def convergence_tests ():

  TL0 = convtest ('TL', 200E4, 0.000, 16, 7, 14, 4)
  BC0 = convtest ('BC', 200E4, 0.000, 16, 7, 14, 4)
  RO0 = convtest ('RO', 200E4, 0.000, 16, 7, 14, 4)
  TL1 = convtest ('TL', 200E9, 1.E-6, 16, 7, 14, 4)
  BC1 = convtest ('BC', 200E9, 1.E-6, 16, 7, 14, 4)
  RO1 = convtest ('RO', 200E9, 1.E-6, 16, 7, 14, 4)

  try:
    import matplotlib.pyplot as plt

    plt.clf ()
    plt.title ('Rotating bar: convergence rate (E = 200E4, $\eta=0$)')
    plt.loglog (TL0[0], TL0[1], label='TL', ls = '-.', lw=3)
    plt.loglog (BC0[0], BC0[1], label='BC')
    plt.loglog (RO0[0], RO0[1], label='RO', ls = '--', marker = 'o')
    plt.xlabel ('Time step $h$ [s]')
    plt.ylabel ('Solution distance $\|q_{ref} - q_h\|$ [m]')
    plt.legend(loc = 'lower right')
    plt.savefig ('out/rotating-bar/rb_undamp_convrate_E200E4.eps')

    plt.clf ()
    plt.title ('Rotating bar: convergence rate (E = 200E9, $\eta=1/10^6$)')
    plt.loglog (TL1[0], TL1[1], label='TL', ls = '-.', lw=3)
    plt.loglog (BC1[0], BC1[1], label='BC')
    plt.loglog (RO1[0], RO1[1], label='RO', ls = '--', marker = 'o')
    plt.xlabel ('Time step $h$ [s]')
    plt.ylabel ('Solution distance $\|q_{ref} - q_h\|$ [m]')
    plt.legend(loc = 'lower right')
    plt.savefig ('out/rotating-bar/rb_undamp_convrate_E200E9.eps')

  except ImportError:
    pass # no reaction

# run undapmped tests
rotation_test (1./64., 1, 0.0, 1, 1, 0)
rotation_test (1./256., 1, 0.0, 1, 1, 0)
rotation_test (1./64., 100, 0.0, 0, 0, 1)
rotation_test (1./64., 1000, 0.0, 0, 0, 1)

# compare rotations
rotation_comparison (2, 5, 12, 1, 200E9, 0.001)

# run damped test
rotation_test (1./64., 1, 0.01, 1, 1, 0)
rotation_test (1./64., 1, 0.05, 1, 1, 0)

# run runtime tests
undamped_rotation_runtimes (1./64., 1)

# run convergence tests
convergence_tests ()

'''
# quasistatic test   
h1 = 1./64.
d1 = 1.
sl1 = SOLFEC ('DYNAMIC', h1, 'out/rotating-bar/QS/RO')
bl1 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = 200E4, poisson = PoissonRatio, density = MassDensity)
bod = BODY (sl1, 'FINITE_ELEMENT', COPY (mesh), bl1, form = 'RO', modal = (eval, evec))
bod.scheme = 'DEF_LIM'
PRESSURE (bod, 0, 10)
PRESSURE (bod, 5, 10)
bod.damping = h1# O(h) damping
RUN (sl1, sv, d1)
'''
