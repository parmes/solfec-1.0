from scipy.linalg import eigh
from scipy.io import mmread
import shutil
import modred
import numpy

# define global path extension
pathext = ''
def ro0_path_extension(path_extension):
  global pathext
  pathext = path_extension

# model creation routine
def ro0_model (step, damping=0.0, femform='TL',
  robase=None, verbose='ON', kind='FINITE_ELEMENT',
  Young=200E4, overwrite=False, runduration=0.0):

  # output path
  if kind == 'FINITE_ELEMENT':
    if femform == 'TL': outpath = 'out/reduced-order0/ro0-fem-tl' + pathext
    elif femform == 'BC': outpath = 'out/reduced-order0/ro0-fem-bc' + pathext
    elif femform == 'BC-MODAL': outpath = 'out/reduced-order0/ro0-modal' + pathext
    elif femform == 'BC-RO': outpath = 'out/reduced-order0/ro0-reduced' + pathext
  elif kind == 'RIGID': outpath = 'out/reduced-order0/ro0-rigid' + pathext

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

# POD base calculation routine (keep rigid modes)
def ro0_POD_base_keep_rigid(rigid, deformable, num_modes=11, verbose=False):

  vecs = numpy.transpose(numpy.array(deformable))

  if verbose:
    svec = vecs.shape[0]
    nvec = vecs.shape[1]
    print 'Calculating', num_modes-6, 'POD modes for', nvec, 'input vectors of size', svec, '...'

  modes, eig_vals = modred.compute_POD_matrices_snaps_method(vecs, list(range(num_modes-6)))
  #modes, eig_vals = modred.compute_POD_matrices_direct_method(vecs, list(range(num_modes-6)))

  if verbose:
    print 'POD eigen values:'
    print eig_vals[0:num_modes-6]

  defo = numpy.transpose(modes).tolist()

  # re-orthogonalize deformable modes with
  # respect to the 6 rigid modes and themselves
  for i in range(0,len(defo)):
    for j in range(0, i): # defo(i) _|_ defo(j<i)
      x = defo[i]
      y = defo[j]
      dot = numpy.dot(x, y)
      z = numpy.array(x) - dot*numpy.array(y)
      defo[i] = z.tolist()
    for y in rigid: # defo(i) _|_ rigid
      x = defo[i]
      dot = numpy.dot(x, y)
      z = numpy.array(x) - dot*numpy.array(y)
      defo[i] = z.tolist()
    # normalize
    x = defo[i]
    invlen = 1.0/numpy.dot(x, x)**0.5
    z = numpy.array(x) * invlen
    defo[i] = z.tolist()

  '''
  for i in range(0,len(defo)):
    for j in range(0,len(defo)):
      x = defo[i]
      y = defo[j]
      dot = numpy.dot(x, y)
      print 'dot(%d,%d) = %g' % (i, j, dot)
    for y in rigid:
      x = defo[i]
      dot = numpy.dot(x, y)
      print 'dot(%d,rigid) = %g' % (i, dot)
   '''

  base = [x for vec in (rigid+defo) for x in vec]
  val = eig_vals.tolist()

  return ([0.]*6+[1.]*(num_modes-6), base)

# retrieve time history of time
def ro0_times (sol, progress='OFF'):
  dur = DURATION (sol)
  th = HISTORY (sol, ['STEP'], dur[0], dur[1], 1, progress)
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

# retrieve energy time history
def ro0_energy (sol):
  bod = sol.bodies[0]
  dur = DURATION (sol)
  th = HISTORY (sol, [(bod, 'KINETIC'), (bod, 'INTERNAL')], dur[0], dur[1], progress='ON')
  tot = []
  for (kinetic,internal) in zip(th[1],th[2]):
    tot.append (kinetic+internal)
  return tot

# calculate vector norm
def norm(q):
  sum = 0.0
  for x in q: sum += x**2
  return sum**0.5

# calculate vector difference
def diff(a, b):
  c = []
  for x, y in zip (a, b): c.append (x-y)
  return c

# convergence rate test
def ro0_convtest (femform, Young, damping, pow0, pow1, pow2, pow3, robase=None):

  h0 = 1.0 / (2.0 ** pow0)
  d0 = 1.0 / (2.0 ** pow3)
  print '\b%s --> h=%g ...' % (femform, h0),
  sol = ro0_model (h0, damping, femform, robase, '%',
                   'FINITE_ELEMENT', Young, 'ON', d0)
  q0 = sol.bodies[0].conf

  dq = []
  hh = []
  for i in range (pow1, pow2):
    h =  1.0 / (2.0 ** i)
    print '\b%s --> h=%g ...' % (femform, h),
    sol = ro0_model (h, damping, femform, robase, '%',
                     'FINITE_ELEMENT', Young, 'ON', d0)
    q = sol.bodies[0].conf
    dq.append (norm(diff(q0, q)))
    hh.append (h)

  print '\b%s --> convergence ratios:' % femform
  for i in range (0, len(dq)-1):
    print dq[i] / dq[i+1]

  return (hh, dq)
