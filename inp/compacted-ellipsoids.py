# compacted ellipsoids example

from random import randint
from random import random
from random import shuffle
from random import seed
from math import sqrt
import time

seed (1) # global seed for repeatibility of input / output

def is_string (s):
  try:
    str(s)
    return True
  except ValueError:
    return False

def is_number(s):
  try:
    float(s)
    return True
  except ValueError:
    return False

def READ_ELLS ():

  ell = []
  rmax = 0
  rmin = 1.0E10
  argv = NON_SOLFEC_ARGV()
  if argv != None and not is_string (argv [0]):  # file format : |---------------------------
							       # | N [...]
							       # | rx_1 ry_1 rz_1 [...]
							       # | [...]
							       # | rx_N ry_N rz_N [...]
    inp =  open (argv [0], 'r')
    if len (argv) > 1: coef = float (argv [1])
    else: coef = 1.0
    lin = inp.readline ()
    lst = lin.split ()
    num = int (lst [0])
    for i in range (0, num):
      lin = inp.readline ()
      lst = lin.split ()
      rx = float (lst [0]) * coef
      ry = float (lst [1]) * coef
      rz = float (lst [2]) * coef
      rmax = max (rmax, rx, ry, rz)
      rmin = min (rmin, rx, ry, rz)
      ell.append ((rx, ry, rz))
    return (ell, rmax, rmin) # output list of tuples of ellipsoid radii

  else: # generate random
    num = 1000
    if argv != None and is_number (argv [0]):
      num = int (argv [0])

    for i in range (0, num):
      r = 0.1 / float (num) ** (1.0/3.0)
      rx = 0.2*r + 0.7*r * random ()
      ry = 0.2*r + 0.7*r * random ()
      rz = 0.2*r + 0.7*r * random ()
      rmax = max (rmax, rx, ry, rz)
      rmin = min (rmin, rx, ry, rz)
      ell.append ((rx, ry, rz))
    return (ell, rmax, rmin) # output list of tuples of ellipsoid radii

def RND():
  return 0.25 * (0.5 - random())

def ONE_RUN (solver, mu, nt_delta, nt_epsilon, nt_maxiter):

  print 'Reading ...',
  (ell, rmax, rmin) = READ_ELLS ()
  print len (ell), 'ellipsoids (rmax = %g, rmin = %g).' % (rmax, rmin)

  shuffle (ell) # apply pseudo-random reordering (ellipsoids might be initially sorted according to volume)

  step = 1E-3 # arbitrary
  dist = rmin * 0.1 # motion per step

  if solver == 'GS':
    ending = 'mu-%g-gs-'%mu
  else:
    ending = 'mu-%g-nt-%g-%g-%d'%(mu,nt_delta, nt_epsilon, nt_maxiter)

  solfec = SOLFEC ('DYNAMIC', step, 'out/compacted-ellipsoids' + ending)
  bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E9, poisson = 0.3, density = 1E3)
  surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = mu, restitution = 0)

  n = 0
  dmax = 2 * rmax
  num = int (pow (float (len (ell)), 1./3.) + 1.0)
  for i in range (0, num):
    for j in range (0, num):
      for k in range (0, num):
	if n < len (ell):
	  x = i * dmax
	  y = j * dmax
	  z = k * dmax
	  radii = ell [n]
	  shp = ELLIP ((x, y, z), radii, 1, 1)
	  b = BODY (solfec, 'RIGID', shp, bulkmat) # ellipsoid bodies
	  INITIAL_VELOCITY (b, (RND(), RND(), RND()), (0, 0, 0))
	  n = n + 1

  # moving boundaries
  box = HULL ([0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1], 2, 2) # unit box

  bod = []
  pnt = []
  for d in range (0, 3):

    if d == 0:
      sca = (dmax, num * dmax, num * dmax)
      tlo = (-dmax-rmax, -rmax, -rmax)
      thi = (num * dmax - rmax, -rmax, -rmax)
      d1 = (0, 1, 0)
      d2 = (0, 0, 1)
    elif d == 1:
      sca = (num * dmax, dmax, num * dmax)
      tlo = (-rmax, -dmax-rmax, -rmax)
      thi = (-rmax, num * dmax - rmax, -rmax)
      d1 = (1, 0, 0)
      d2 = (0, 0, 1)
    else:
      sca = (num * dmax, num * dmax, dmax)
      tlo = (-rmax, -rmax, -dmax-rmax)
      thi = (-rmax, -rmax, num * dmax - rmax)
      d1 = (1, 0, 0)
      d2 = (0, 1, 0)

    dst = 5 * dmax

    lo = COPY (box)
    SCALE (lo, sca)
    TRANSLATE (lo, tlo)
    dir = [0, 0, 0]
    dir [d] = 1
    pt1 = MASS_CENTER (lo)
    pnq = TRANSLATE (pt1, (-dst*dir[0], -dst*dir[1], -dst*dir[2]))
    s = SPHERE (pnq, 2*dmax, 2, 2) # sphere to put anti-rotational constraints
    b = BODY (solfec, 'OBSTACLE', [lo, s], bulkmat)
    SET_VELOCITY (b, pt1, tuple (dir), dist / step)
    FIX_DIRECTION (b, pt1, d1)
    FIX_DIRECTION (b, pt1, d2)
    FIX_DIRECTION (b, pnq, d1)
    FIX_DIRECTION (b, pnq, d2)
    bod.append (b)
    pnt.append (pt1)

    hi = COPY (box)
    SCALE (hi, sca)
    TRANSLATE (hi, thi)
    dir = [0, 0, 0]
    dir [d] = -1
    pt2 = MASS_CENTER (hi)
    pnq = TRANSLATE (pt2, (-dst*dir[0], -dst*dir[1], -dst*dir[2]))
    s = SPHERE (pnq, 2*dmax, 2, 2) # sphere to put anti-rotational constraints
    b = BODY (solfec, 'OBSTACLE', [hi, s], bulkmat)
    SET_VELOCITY (b, pt2, tuple (dir), dist / step)
    FIX_DIRECTION (b, pt2, d1)
    FIX_DIRECTION (b, pt2, d2)
    FIX_DIRECTION (b, pnq, d1)
    FIX_DIRECTION (b, pnq, d2)
    bod.append (b)
    pnt.append (pt2)

  for i in range (0, len (bod)):
    for j in range (0, i):
      CONTACT_EXCLUDE_BODIES (bod [i], bod [j])

  dst = sqrt ((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2 + (pt2[2]-pt1[2])**2) - dmax
  print 'Initial distance is ', dst
  dif = dst - 0.1
  if dif > 0.0:
    stop = 0.5 * step * (dif / dist) # 0.5 since walls from both directions move
    print 'The simulation time is ', stop, ', that is %d time steps' % (stop / step)
  else:
    stop = 1.0

  if solver == 'GS':
    slv = GAUSS_SEIDEL_SOLVER (1, 1000, 1E-8)
  else:
    slv = NEWTON_SOLVER (1E-8, 500, maxmatvec = 1000, delta = nt_delta, epsilon = nt_epsilon, linmaxiter = nt_maxiter)

  IMBALANCE_TOLERANCE (solfec, 1.3, 0.5, 10)
  OUTPUT (solfec, stop / 20.0)

  t0 = time.time()
  RUN (solfec, slv, stop)
  elapsed = time.time() - t0
  if RANK() == 0 and solfec.mode == 'WRITE':   
      print "Analysis run time =", elapsed/3600.0, "hours"

  # post-processing
  if not VIEWER() and solfec.mode == 'READ':
    dur = DURATION (solfec)
    SEEK (solfec, dur [1])
    d = [1, 1, 1]
    for i in range (0, 3):
      p = pnt [2*i]
      q = pnt [2*i+1]
      a = DISPLACEMENT (bod [2*i], p)
      b = DISPLACEMENT (bod [2*i+1], q)
      p = TRANSLATE (p, a)
      q = TRANSLATE (q, b)
      d [i] = sqrt ((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2) - dmax
    l = (d[0]+d[1]+d[2])/3.0
    print 'Final distance is ', l

# main module

ONE_RUN ('GS', 0.0, 0, 0.25, 10)
