from scipy.linalg import eigh
from scipy.io import mmread
import shutil
import modred
import numpy

# model creation routine
def ro0_model (step, damping=0.0, femform='TL',
  robase=None, verbose='ON', kind='FINITE_ELEMENT',
  Young=200E4, overwrite=False, runduration=0.0):

  # output path
  if kind == 'FINITE_ELEMENT':
    if femform == 'TL': outpath = 'out/reduced-order0/ro0-fem-tl'
    elif femform == 'BC': outpath = 'out/reduced-order0/ro0-fem-bc'
    elif femform == 'BC-MODAL': outpath = 'out/reduced-order0/ro0-modal'
    elif femform == 'BC-RO': outpath = 'out/reduced-order0/ro0-reduced'
  elif kind == 'RIGID': outpath = 'out/reduced-order0/ro0-rigid'

  # remove previous output directory
  if overwrite:
    shutil.rmtree (outpath, True)

  # SOLFEC simulation object
  sol = SOLFEC ('DYNAMIC', step, outpath)
  sol.verbose = verbose

  # material
  mat = BULK_MATERIAL (sol, model = 'KIRCHHOFF',
    young = Young, poisson = 0.26, density = 7.8E3)

  # bar geometry
  nodes = [-0.05, -0.05, -0.5,
	    0.05, -0.05, -0.5,
	    0.05,  0.05, -0.5,
	   -0.05,  0.05, -0.5,
	   -0.05, -0.05,  0.5,
	    0.05, -0.05,  0.5,
	    0.05,  0.05,  0.5,
	   -0.05,  0.05,  0.5]
  mesh = HEX (nodes, 2, 2, 20, 0, [0]*6)

  # create body
  if kind == 'FINITE_ELEMENT':
    if robase != None:
      bod = BODY (sol, kind, mesh, mat, form=femform, base=robase)
    else: bod = BODY (sol, kind, mesh, mat, form=femform)
    bod.scheme = 'DEF_LIM'
    bod.damping = damping
  else: bod = BODY (sol, kind, mesh, mat)

  # initial angular velocity
  INITIAL_VELOCITY (bod, (0, 0, 0), (1, 0, 0))

  # run if requested
  if runduration > 0.0:
    RUN (sol, NEWTON_SOLVER(), runduration)

  # return simulation object
  return sol

# modal base calculation routine
def ro0_modal_base (use_scipy=False, verbose='OFF'):
  sol = ro0_model (1E-3, 0.0)
  bod = sol.bodies[0]
  eval = [] # selected eigenvalue list
  evec = [] # selected eigenvector list
  vsel = (0,1,2,3,4,5,13,18,25,33,38)
  if use_scipy:
    BODY_MM_EXPORT (bod, 'out/reduced-order0/M.mtx',
                         'out/reduced-order0/K.mtx')
    M = mmread ('out/reduced-order0/M.mtx').todense()
    K = mmread ('out/reduced-order0/K.mtx').todense()
    for j in range (0, K.shape[1]):
      for i in range (j+1, K.shape[0]):
	K [j, i] = K [i, j] # above diagonal = below diagonal
    x, y = eigh (K, M) # this produces y.T M y = 1 and y.T K y = x
    for j in vsel:
      eval.append (x[j].real)
      for z in y[:,j]:
	evec.append (z.real)
  else:
    data0 = MODAL_ANALYSIS (bod, 45, 'out/reduced-order0/modal',
                            1E-13, 1000, verbose)
    dofs = len(bod.velo)
    for j in vsel:
      eval.append (data0[0][j])
      for k in range (j*dofs,(j+1)*dofs):
	evec.append (data0[1][k])

  return (eval, evec)

# POD base calculation routine
def ro0_POD_base(rigid, deformable, num_modes=11, verbose=False):

  vecs = numpy.transpose(numpy.array(rigid+deformable))

  if verbose:
    svec = vecs.shape[0]
    nvec = vecs.shape[1]
    print 'Calculating', num_modes, 'POD modes for', nvec, 'input vectors of size', svec, '...'

  modes, eig_vals = modred.compute_POD_matrices_snaps_method(vecs, list(range(num_modes)))
  #modes, eig_vals = modred.compute_POD_matrices_direct_method(vecs, list(range(num_modes)))

  if verbose:
    print 'POD eigen values:'
    print eig_vals[0:num_modes]

  mod = numpy.transpose(modes).tolist()
  base = [x for vec in mod for x in vec]
  val = eig_vals.tolist()

  return (val[0:len(mod)], base)

# retrieve time history of time
def ro0_times (sol):
  dur = DURATION (sol)
  th = HISTORY (sol, ['STEP'], dur[0], dur[1])
  return th[0]

# retrieve elongation time history
def ro0_elongation (sol):
  p0 = (0.0, 0.0, -0.5)
  p1 = (0.0, 0.0, 0.5)
  bod = sol.bodies[0]
  dur = DURATION (sol)
  th = HISTORY (sol, [(bod, p0, 'DX'), (bod, p0, 'DY'), (bod, p0, 'DZ'), \
      (bod, p1, 'DX'), (bod, p1, 'DY'), (bod, p1, 'DZ')], dur[0], dur[1])
  lh = []
  for (dx0,dy0,dz0,dx1,dy1,dz1) in zip (th[1], th[2], th[3], th[4], th[5], th[6]):
    l = ((p0[0]+dx0-p1[0]-dx1)**2 + (p0[1]+dy0-p1[1]-dy1)**2 + (p0[2]+dz0-p1[2]-dz1)**2)**0.5
    lh.append (l-1.0)
  return lh
