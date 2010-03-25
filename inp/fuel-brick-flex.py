# fuel brick with flexible suspension
import sys
sys.path.append ('inp/cores/inc')
from simple_core_base import *

def run_analysis (step, stop, scheme):
  solfec = SOLFEC ('DYNAMIC', step, 'out/fuel-brick-flex/' + scheme)
  SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3)
  graphite = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3)
  steel = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 200E9, poisson = 0.28, density = 8.0E3)
  GRAVITY (solfec, (0, 0, -10))

  shp = gcore_brick (0, 0, 0)
  a = BODY (solfec , 'RIGID', shp, graphite, label = 'A')

  vertices = [1, 0, 0,
	      1, 1, 0,
	      0, 1, 0,
	      0, 0, 0,
	      1, 0, 1,
	      1, 1, 1,
	      0, 1, 1,
	      0, 0, 1]

  outd = 0.4598
  height = 0.225
  hstep = 0.0098
  eps = 5E-3 * outd

  msh = HEX (vertices, 2, 2, 2, 1, [1, 1, 1, 1, 1, 1])
  SCALE (msh, (0.5, 0.1, 0.1))
  TRANSLATE (msh, (0.5 * outd - eps, -0.05, eps))
  p1 = msh.node (26)
  p2 = msh.node (24) 
  p3 = msh.node (8) 
  p4 = msh.node (6) 
  p5 = msh.node (20)
  p6 = msh.node (18) 
  p7 = msh.node (2) 
  p8 = msh.node (0) 
  b = BODY (solfec , 'FINITE_ELEMENT', msh, steel, label = 'B')
  b.scheme = scheme
  CONTACT_EXCLUDE_BODIES (solfec, a, b)
  PUT_RIGID_LINK (b, a, p1, p1)
  PUT_RIGID_LINK (b, a, p2, p2)
  PUT_RIGID_LINK (b, a, p3, p3)
  PUT_RIGID_LINK (b, a, p4, p4)
  FIX_POINT (b, p5)
  FIX_POINT (b, p6)
  FIX_POINT (b, p7)
  FIX_POINT (b, p8)

  z = (1.1) * (2*height - hstep)
  shp = gcore_brick (0, 0, z)
  BODY (solfec , 'RIGID', shp, graphite)

  solver = GAUSS_SEIDEL_SOLVER (1E-6, 100, failure = 'CONTINUE')
  RUN (solfec, solver, stop)
  return solfec

# main module
step = 1E-3
stop = 0.5

solfec1 = run_analysis (step, stop, 'DEF_EXP')
solfec2 = run_analysis (step, stop, 'DEF_LIM')
solfec3 = run_analysis (step, stop, 'DEF_LIM2')
solfec4 = run_analysis (step, stop, 'DEF_IMP')

if not VIEWER() and solfec1.mode == 'READ':

  list = [solfec1, solfec2, solfec3, solfec4]

  try:
    import matplotlib.pyplot as plt

    for solfec in list:
      a = BYLABEL (solfec, 'BODY', 'A')
      b = BYLABEL (solfec, 'BODY', 'B')
      th = HISTORY (solfec, [(a, (0, 0, 0), 'DZ')], 0, stop)
      plt.plot (th [0], th [1], label = b.scheme)

    plt.axis (xmin = 0, xmax = stop)
    plt.xlabel ('Time [s]')
    plt.ylabel ('Displacement [m]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/fuel-brick-flex/fuel-brick-flex-dz.eps')
  except ImportError:
    pass # no reaction

  for solfec in list:
    b = BYLABEL (solfec, 'BODY', 'B')
    timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL']
    dur = DURATION (solfec)
    th = HISTORY (solfec, timers, dur[0], dur[1])
    total = 0.0

    print 'SCHEME ' + b.scheme + ':'

    for i in range (0, len(timers)):
      sum = 0.0
      for tt in th [i+1]: sum += tt
      print timers [i], 'TIME:', sum
      total += sum

    print 'TOTAL TIME:', total

    print '--------------------------'
