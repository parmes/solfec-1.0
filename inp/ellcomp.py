# ellipsoid compaction test

from random import randint
from random import random
from random import shuffle
from random import seed

seed (1) # global seed for repeatibility of input / output

def READ_ELLS ():

  ell = []
  rmax = 0
  argv = NON_SOLFEC_ARGV()
  if argv != None: # read from file => file format: |---------------------------
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
      ell.append ((rx, ry, rz))
    return (ell, rmax) # output list of tuples of ellipsoid radii

  else: # generate random
    for i in range (0, 1000):
      rx = random ()
      ry = random ()
      rz = random ()
      rmax = max (rmax, rx, ry, rz)
      ell.append ((rx, ry, rz))
    return (ell, rmax) # output list of tuples of ellipsoid radii


def COMPACTION_TEST (name, step):

  solfec = SOLFEC ('DYNAMIC', step, 'out/ellcomp_' + name)
  bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E9, poisson = 0.3, density = 1E3)
  surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0, restitution = 0)

  print 'Reading ...',
  (ell, rmax) = READ_ELLS ()
  print len (ell), 'ellipsoids (rmax = %g).' % (rmax)

  shuffle (ell) # apply pseudo-random reordering (ellipses might be initially sorted according to volume)

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
          BODY (solfec, 'RIGID', shp, bulkmat) # ellipsoid bodies
	  n = n + 1

  # force value
  fval = bulkmat.density * (dmax ** 3) * 1E4 # FIXME: smarten this up

  # moving boundaries
  box = HULL ([0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1], 2, 2) # unit box

  bod = []
  for d in range (0, 3):

    if d == 0:
      sca = (dmax, num * dmax, num * dmax)
      tlo = (-dmax, -rmax, -rmax)
      thi = (num * dmax - rmax, -rmax, -rmax)
      d1 = (0, 1, 0)
      d2 = (0, 0, 1)
    elif d == 1:
      sca = (num * dmax, dmax, num * dmax)
      tlo = (-rmax, -dmax, -rmax)
      thi = (-rmax, num * dmax - rmax, -rmax)
      d1 = (1, 0, 0)
      d2 = (0, 0, 1)
    else:
      sca = (num * dmax, num * dmax, dmax)
      tlo = (-rmax, -rmax, -dmax)
      thi = (-rmax, -rmax, num * dmax - rmax)
      d1 = (1, 0, 0)
      d2 = (0, 1, 0)

    dst = 5 * dmax

    lo = COPY (box)
    SCALE (lo, sca)
    TRANSLATE (lo, tlo)
    dir = [0, 0, 0]
    dir [d] = 1
    pnt = MASS_CENTER (lo)
    pnq = TRANSLATE (pnt, (-dst*dir[0], -dst*dir[1], -dst*dir[2]))
    s = SPHERE (pnq, 2*dmax, 2, 2) # sphere to put anti-rotational constraints
    b = BODY (solfec, 'RIGID', [lo, s], bulkmat)
    FORCE (b, 'SPATIAL', pnt, tuple (dir), fval)
    FIX_DIRECTION (b, pnq, d1)
    FIX_DIRECTION (b, pnq, d2)
    bod.append (b)

    hi = COPY (box)
    SCALE (hi, sca)
    TRANSLATE (hi, thi)
    dir = [0, 0, 0]
    dir [d] = -1
    pnt = MASS_CENTER (hi)
    pnq = TRANSLATE (pnt, (-dst*dir[0], -dst*dir[1], -dst*dir[2]))
    s = SPHERE (pnq, 2*dmax, 2, 2) # sphere to put anti-rotational constraints
    b = BODY (solfec, 'RIGID', [hi, s], bulkmat)
    FORCE (b, 'SPATIAL', pnt, tuple (dir), fval)
    FIX_DIRECTION (b, pnq, d1)
    FIX_DIRECTION (b, pnq, d2)
    bod.append (b)

  for i in range (0, len (bod)):
    for j in range (0, i):
      CONTACT_EXCLUDE_BODIES (bod [i], bod [j])

  return solfec

# main module

stop = 2.5
slv = GAUSS_SEIDEL_SOLVER (1E-3, 20)
sol = COMPACTION_TEST ('a', 1E-3)
RUN (sol, slv, stop)
