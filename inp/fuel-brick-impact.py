# simple core model fuel brick drop
import sys
sys.path.append ('inp/cores/inc')
from simple_core_base import *

PI = 3.14159265358979323846 

def scene_base (material, solfec):

  vertices = [1, 0, 0,
              1, 1, 0,
              0, 1, 0,
              0, 0, 0,
              1, 0, 1,
              1, 1, 1,
              0, 1, 1,
              0, 0, 1]

  faces = [4, 0, 1, 5, 4, 3,
           4, 1, 2, 6, 5, 3,
           4, 2, 3, 7, 6, 3,
           4, 3, 0, 4, 7, 3,
           4, 0, 3, 2, 1, 3,
           4, 4, 5, 6, 7, 3]

  outd = 0.4598
  margin = 0.05
  thick = 0.1
  lx = outd + (margin + thick)
  ly = outd + (margin + thick)

  cvx = CONVEX (vertices, faces, 3)
  scl = (lx,  ly,  thick)
  vec = (-lx/2, -ly/2, -thick)
  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)

  BODY (solfec, 'OBSTACLE', cvx, material)

def scene_run (step, stop, kinem, formul, damp):

  solver = GAUSS_SEIDEL_SOLVER (1E-4, 1000)
  solfec = SOLFEC ('DYNAMIC', step, 'out/fuel-brick-impact/' + kinem + '_' + str(step) + '_' + str (stop) + '_' + formul + '_' + str (damp))
  bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E10, poisson = 0.2, density = 2E3)
  SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = 0.0)
  GRAVITY (solfec, (0, 0, -10))

  scene_base (bulkmat, solfec)
  shp = gcore_brick (0, 0, 1)
  msh = PIPE ((0, 0, 1), (0, 0, 0.45), 0.125, 0.13, 3, 8, 2, 0, [0, 0, 0, 0])
  bod = BODY (solfec , kinem, shp, bulkmat, mesh = msh, form = formul, label = 'FU')
  if kinem != 'RIGID': bod.scheme = 'DEF_LIM'
  #INITIAL_VELOCITY (bod, (0, 0, -1), (0, 0, 0))
  bod.damping = damp

  RUN (solfec, solver, stop)

  return (solfec, kinem)

# main module
#import rpdb2; rpdb2.start_embedded_debugger('a')

stop = 0.8
pairs = []
pairs.append (scene_run (0.0010, stop, 'FINITE_ELEMENT', 'BC', 0))
pairs.append (scene_run (0.0005, stop, 'FINITE_ELEMENT', 'BC', 0))
pairs.append (scene_run (0.0002, stop, 'FINITE_ELEMENT', 'BC', 0))
pairs.append (scene_run (0.0001, stop, 'FINITE_ELEMENT', 'BC', 0))

pairs.append (scene_run (0.0005, stop, 'FINITE_ELEMENT', 'BC', 1E-5))
pairs.append (scene_run (0.0005, stop, 'FINITE_ELEMENT', 'BC', 5E-5))
pairs.append (scene_run (0.0005, stop, 'FINITE_ELEMENT', 'BC', 1E-4))
pairs.append (scene_run (0.0005, stop, 'FINITE_ELEMENT', 'BC', 1E-3))

if not VIEWER() and pairs [0][0].mode == 'READ':

  try:
    import matplotlib.pyplot as plt

    d = []
    for i in range (0, 4):
      solfec = pairs [i][0]
      b = BYLABEL (solfec, 'BODY', 'FU')
      pnt = (0.1315, 0, 1.45)
      th = HISTORY (solfec, [(b, pnt, 'DZ')], 0.3, 0.8)
      d.append (th[0])
      d.append (th[1])

    plt.clf ()
    plt.title ( 'Various time steps, CR-FE')
    plt.plot (d[0], d[1], label = '$h=1^{-3}$')
    plt.plot (d[2], d[3], label = '$h=5^{-4}$')
    plt.plot (d[4], d[5], label = '$h=2^{-4}$')
    plt.plot (d[6], d[7], label = '$h=2^{-4}$')
    plt.xlabel ('Time [s]')
    plt.ylabel ('Displacement DZ [m]')
    plt.legend(loc = 'lower left')
    plt.savefig ('out/fuel-brick-impact/fb_dz_varstep.eps')
  except ImportError:
    pass # no reaction

  try:
    import matplotlib.pyplot as plt

    d = []
    for i in (1, 4, 5, 6, 7):
      solfec = pairs [i][0]
      b = BYLABEL (solfec, 'BODY', 'FU')
      pnt = (0.1315, 0, 1.45)
      th = HISTORY (solfec, [(b, pnt, 'DZ')], 0.3, 0.8)
      d.append (th[1])

    plt.clf ()
    plt.title ( 'Step 0.0005, CR-FE')
    plt.plot (th [0], d[0], label = '$\eta=0$')
    plt.plot (th [0], d[1], label = '$\eta=1^{-5}$')
    plt.plot (th [0], d[2], label = '$\eta=5^{-5}$')
    plt.plot (th [0], d[3], label = '$\eta=1^{-4}$')
    plt.plot (th [0], d[4], label = '$\eta=1^{-3}$')
    plt.xlabel ('Time [s]')
    plt.ylabel ('Displacement DZ [m]')
    plt.legend(loc = 'lower left')
    plt.savefig ('out/fuel-brick-impact/fb_dz_vareta.eps')
  except ImportError:
    pass # no reaction
