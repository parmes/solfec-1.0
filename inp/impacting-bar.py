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

# base middle point
p0 = MASS_CENTER (mesh)

# obstacle mesh
obsm = HEX (nodes, 1, 1, 1, 0, [0, 1, 2, 3, 4, 5])
SCALE (obsm, (2, 2, 0.2))
TRANSLATE (obsm, (0, 0, -0.6))

# solver
sv = NEWTON_SOLVER ()

# compute eigen base
def eigenbase (h1, d1, E, v0, pow0, pow1, rest):
  # compute all eigenvalues and eigenvectors
  pt0 = 'out/impacting-bar/MK_%g_%g_%g_%g_%d_%d'%(h1, d1, E, v0, pow0, pow1)
  sl0 = SOLFEC ('DYNAMIC', 1E-3, pt0)
  bl0 = BULK_MATERIAL (sl0, model = 'KIRCHHOFF', young = E, poisson = PoissonRatio, density = MassDensity)
  bod = BODY (sl0, 'FINITE_ELEMENT', COPY (mesh), bl0)
  eval = [] # selected eigenvalue list
  evec = [] # selected eigenvector list (BODY command takes a tuple (eval, evec) argument for the RO formulation)
  vsel = (0,1,2,3,4,5,13,18,25,33,38)
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
  return (eval, evec)

# impact comparison
def impact_comparison (h1, d1, E, v0, pow0, pow1, rest):

  # compute eigen base
  data = eigenbase (h1, d1, E, v0, pow0, pow1, rest)

  toplot = []
  vdamp = []
  vout_BC = []
  vout_RO = []
  for p in range (pow0, pow1):

    damping = 1.0 / (2.0**p)

    vdamp.append (damping)

    # rotation: BC
    sl1 = SOLFEC ('DYNAMIC', h1, 'out/impacting-bar/BC_%g_%g_%g_%g_%d_%d_%g_%g'%(h1, d1, E, v0, pow0, pow1, damping, rest))
    SURFACE_MATERIAL (sl1, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = rest)
    bl1 = BULK_MATERIAL (sl1, model = 'KIRCHHOFF', young = E, poisson = PoissonRatio, density = MassDensity)
    bd1 = BODY (sl1, 'FINITE_ELEMENT', COPY (mesh), bl1, form = 'BC')
    bd1.scheme = 'DEF_LIM'
    bd1.damping = damping
    INITIAL_VELOCITY (bd1, (0, 0, v0), (0, 0, 0))
    BODY (sl1, 'OBSTACLE', COPY (obsm), bl1)
    RUN (sl1, sv, d1)
    if not VIEWER() and sl1.mode == 'READ':
      th1 = HISTORY (sl1, [(bd1, p0, 'VZ')], h1, d1)
      vout = 0.0
      for i in range (len(th1[1])-8, len(th1[1])):
	vout += th1[1][i]
      vout_BC.append (vout/8.0)



    # rotation: RO
    sl2 = SOLFEC ('DYNAMIC', h1, 'out/impacting-bar/RO_%g_%g_%g_%g_%d_%d_%g_%g'%(h1, d1, E, v0, pow0, pow1, damping, rest))
    SURFACE_MATERIAL (sl2, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = rest)
    bl2 = BULK_MATERIAL (sl2, model = 'KIRCHHOFF', young = E, poisson = PoissonRatio, density = MassDensity)
    bd2 = BODY (sl2, 'FINITE_ELEMENT', COPY (mesh), bl2, form = 'RO', modal = data)
    bd2.scheme = 'DEF_LIM'
    bd2.damping = damping
    INITIAL_VELOCITY (bd2, (0, 0, v0), (0, 0, 0))
    BODY (sl2, 'OBSTACLE', COPY (obsm), bl2)
    RUN (sl2, sv, d1)
    if not VIEWER() and sl2.mode == 'READ':
      th2 = HISTORY (sl2, [(bd2, p0, 'VZ')], h1, d1)
      vout = 0.0
      for i in range (len(th2[1])-8, len(th2[1])):
	vout += th2[1][i]
      vout_RO.append (vout/8.0)
    
    if not VIEWER() and sl1.mode == 'READ' and sl2.mode == 'READ':
      toplot.append (('$\eta=1/' + str(2**p) + '$', th1[0], th1[1], th2[1]))

  if not VIEWER() and sl1.mode == 'READ' and sl2.mode == 'READ' and 0: # comment out 0 to get these plots

    try:
      import matplotlib.pyplot as plt

      hstr = ' ($E = %g, h = 1/%d'%(E, 1/h1) + '$)'

      plt.clf ()
      plt.title ('Impacting bar: base point VZ: BC' + hstr)
      for dat in toplot:
        plt.plot (dat[1], dat[2], label=dat[0])
      plt.xlabel ('Time [s]')
      plt.ylabel ('VZ [m/s]')
      plt.legend(loc = 'best')
      plt.savefig ('out/impacting-bar/ib_vz_BC_h%d_d%d_E%g_v%g.eps'%(1/h1,1/d1,E,v0))

      plt.clf ()
      plt.title ('Impacting bar: base point VZ: BC' + hstr)
      for dat in toplot:
        plt.plot (dat[1], dat[3], label=dat[0])
      plt.xlabel ('Time [s]')
      plt.ylabel ('VZ [m/s]')
      plt.legend(loc = 'best')
      plt.savefig ('out/impacting-bar/ib_vz_RO_h%d_d%d_E%g_v%g.eps'%(1/h1,1/d1,E,v0))

    except ImportError:
      pass # no reaction

  return ('$h=\\frac{1}{%d}$'%(1/h1), vdamp, vout_BC, vout_RO)

#energy comparison
def enery_comparison (n1, h1, d1, E, v0, pow0, rest):

  if pow0 > 0.0: damping = 1.0 / (2.0**pow0)
  else: damping = 0.0

  # local mesh
  mesh = HEX (nodes, n1, n1, n1*10, 0, [0, 1, 2, 3, 4, 5])

  # rotation: BC
  sl1 = SOLFEC ('DYNAMIC', h1, 'out/impacting-bar/BC_%g_%g_%g_%g_%d_%d_%g_%g'%(h1, d1, E, v0, pow0, 0, damping, rest))
  SURFACE_MATERIAL (sl1, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = rest)
  bl1 = BULK_MATERIAL (sl1, model = 'KIRCHHOFF', young = E, poisson = PoissonRatio, density = MassDensity)
  bd1 = BODY (sl1, 'FINITE_ELEMENT', TRANSLATE(COPY (mesh), (0, 0, 0.005)), bl1, form = 'BC')
  bd1.scheme = 'DEF_LIM'
  bd1.damping = damping
  INITIAL_VELOCITY (bd1, (0, 0, v0), (0, 0, 0))
  BODY (sl1, 'OBSTACLE', COPY (obsm), bl1)
  RUN (sl1, sv, d1)

  if not VIEWER() and sl1.mode == 'READ':

    th1 = HISTORY (sl1, [(bd1, 'KINETIC'), (bd1, 'INTERNAL')], 0, d1)

    tot1 = []
    for (ek, ei) in zip (th1[1], th1[2]):
      tot1.append (ek + ei)

    try:
      import matplotlib.pyplot as plt

      if damping == 0.0:
	dmp = 0
        hstr = ' ($h = 1/%d, \eta=0, e=%g'%(1/h1, rest) + '$)'
      else:
	dmp = 1/damping
        hstr = ' ($h = 1/%d, \eta=1/%d, e=%g'%(1/h1, dmp, rest) + '$)'


      plt.clf ()
      plt.title ('Impacting %dx%dx%d bar energy: BC'%(n1,n1,n1*10) + hstr)
      plt.plot (th1[0], th1[1], label='KIN', ls='--', marker = 'o')
      plt.plot (th1[0], th1[2], label='INT')
      plt.plot (th1[0], tot1, label='TOT', lw=2)
      plt.xlabel ('Time [s]')
      plt.ylabel ('Energy [J]')
      plt.legend (loc = 'best')
      plt.savefig ('out/impacting-bar/ib_ene_BC_h%d_d%d_r%g.eps'%(1/h1,dmp,rest))

    except ImportError:
      pass # no reaction

  if n1 == 2:
    # compute eigen base
    data = eigenbase (h1, d1, E, v0, pow0, 0, rest)

    # rotation: RO
    sl2 = SOLFEC ('DYNAMIC', h1, 'out/impacting-bar/RO_%g_%g_%g_%g_%d_%d_%g_%g'%(h1, d1, E, v0, pow0, 0, damping, rest))
    SURFACE_MATERIAL (sl2, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = rest)
    bl2 = BULK_MATERIAL (sl2, model = 'KIRCHHOFF', young = E, poisson = PoissonRatio, density = MassDensity)
    bd2 = BODY (sl2, 'FINITE_ELEMENT', TRANSLATE(COPY (mesh), (0, 0, 0.005)), bl2, form = 'RO', modal = data)
    bd2.scheme = 'DEF_LIM'
    bd2.damping = damping
    INITIAL_VELOCITY (bd2, (0, 0, v0), (0, 0, 0))
    BODY (sl2, 'OBSTACLE', COPY (obsm), bl2)
    RUN (sl2, sv, d1)

    if not VIEWER() and sl2.mode == 'READ':

      th2 = HISTORY (sl2, [(bd2, 'KINETIC'), (bd2, 'INTERNAL')], 0, d1)

      tot2 = []
      for (ek, ei) in zip (th2[1], th2[2]):
	tot2.append (ek + ei)

      try:
	import matplotlib.pyplot as plt

	if damping == 0.0:
	  dmp = 0
	  hstr = ' ($h = 1/%d, \eta=0, e=%g'%(1/h1, rest) + '$)'
	else:
	  dmp = 1/damping
	  hstr = ' ($h = 1/%d, \eta=1/%d, e=%g'%(1/h1, dmp, rest) + '$)'

	plt.clf ()
	plt.title ('Impacting bar energy: RO' + hstr)
	plt.plot (th2[0], th2[1], label='KIN')
	plt.plot (th2[0], th2[2], label='INT')
	plt.plot (th2[0], tot2, label='TOT')
	plt.xlabel ('Time [s]')
	plt.ylabel ('Energy [J]')
	plt.legend (loc = 'best')
	plt.savefig ('out/impacting-bar/ib_ene_RO_h%d_d%d_r%g.eps'%(1/h1,dmp,rest))

      except ImportError:
	pass # no reaction

# impacting testes
e0 = 8
e1 = 17
dat0 = impact_comparison (1/1024., 1/32., 200E9, -1.0, e0, e1, 0)
# comment out below before using the viewer (FIXME: interaction with multiple nested RUN commands)
dat1 = impact_comparison (1/2048., 1/32., 200E9, -1.0, e0, e1, 0)
dat2 = impact_comparison (1/4096., 1/32., 200E9, -1.0, e0, e1, 0)
dat3 = impact_comparison (1/1024., 1/32., 200E9, -1.0, e0, e1, 1)
dat4 = impact_comparison (1/2048., 1/32., 200E9, -1.0, e0, e1, 1)
dat5 = impact_comparison (1/4096., 1/32., 200E9, -1.0, e0, e1, 1)

enery_comparison (2, 1/4096., 1/32., 200E9, -1.0, 0, 0)
enery_comparison (2, 1/4096., 1/32., 200E9, -1.0, 0, 1)
enery_comparison (2, 1/4096., 1/32., 200E9, -1.0, 16, 0)
enery_comparison (2, 1/4096., 1/32., 200E9, -1.0, 16, 1)
enery_comparison (2, 1/4096., 1/32., 200E9, -1.0, 10, 0)
enery_comparison (2, 1/4096., 1/32., 200E9, -1.0, 10, 1)

enery_comparison (2*2, 1/(2*4096.), 1/32., 200E9, -1.0, 0, 0)
enery_comparison (4*2, 1/(4*4096.), 1/32., 200E9, -1.0, 0, 0)
enery_comparison (2*2, 1/(2*4096.), 1/32., 200E9, -1.0, 0, 1)
enery_comparison (4*2, 1/(4*4096.), 1/32., 200E9, -1.0, 0, 1)

# avg. output velocity plots
if not VIEWER() and len(dat0[2]) > 0:

  try:
    import matplotlib.pyplot as plt

    plt.clf ()
    plt.title ('Impacting bar: avg. output velocity $u_z$, ($e=0$)')
    plt.plot (dat0[1], dat0[2], label='BC: '+dat0[0] + ', $e=0$')
    plt.plot (dat0[1], dat0[3], label='RO: '+dat0[0] + ', $e=0$', ls = '--', marker = 'o')
    plt.plot (dat1[1], dat1[2], label='BC: '+dat1[0] + ', $e=0$')
    plt.plot (dat1[1], dat1[3], label='RO: '+dat1[0] + ', $e=0$', ls = '--', marker = 'o')
    plt.plot (dat2[1], dat2[2], label='BC: '+dat2[0] + ', $e=0$')
    plt.plot (dat2[1], dat2[3], label='RO: '+dat2[0] + ', $e=0$', ls = '--', marker = 'o')
    xtic = []
    xlab = []
    for p in range (e0, min (e0+5, e1)):
      xtic.append (1./2.**p)
      xlab.append ('$\\frac{1}{%d}$'%2**p)
    plt.xticks (xtic, xlab)
    plt.xlabel ('Damping $\eta$')
    plt.ylabel ('Output velocity $u_z$ [m/s]')
    plt.legend(loc = 'best')
    plt.savefig ('out/impacting-bar/ib_e0_vz_out.eps')

    plt.clf ()
    plt.title ('Impacting bar: avg. output velocity $u_z$, ($e=1$)')
    plt.plot (dat3[1], dat3[2], label='BC: '+dat3[0] + ', $e=1$')
    plt.plot (dat3[1], dat3[3], label='RO: '+dat3[0] + ', $e=1$', ls = '--', marker = 'o')
    plt.plot (dat4[1], dat4[2], label='BC: '+dat4[0] + ', $e=1$')
    plt.plot (dat4[1], dat4[3], label='RO: '+dat4[0] + ', $e=1$', ls = '--', marker = 'o')
    plt.plot (dat5[1], dat5[2], label='BC: '+dat5[0] + ', $e=1$')
    plt.plot (dat5[1], dat5[3], label='RO: '+dat5[0] + ', $e=1$', ls = '--', marker = 'o')
    xtic = []
    xlab = []
    for p in range (e0, min (e0+5, e1)):
      xtic.append (1./2.**p)
      xlab.append ('$\\frac{1}{%d}$'%2**p)
    plt.xticks (xtic, xlab)
    plt.xlabel ('Damping $\eta$')
    plt.ylabel ('Output velocity $u_z$ [m/s]')
    plt.legend(loc = 'best')
    plt.savefig ('out/impacting-bar/ib_e1_vz_out.eps')

  except ImportError:
    pass # no reaction
