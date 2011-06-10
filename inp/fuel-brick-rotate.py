# simple core model fuel brick drop
import sys
sys.path.append ('inp/cores/inc')
from simple_core_base import *

PI = 3.14159265358979323846 

def scene_run (step, stop, kinem, formul, damp):

  solver = GAUSS_SEIDEL_SOLVER (1E-4, 1000)
  solfec = SOLFEC ('DYNAMIC', step, 'out/fuel-brick-rotate/' + kinem + '_' + str(step) + '_' + str (damp))
  bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E10, poisson = 0.2, density = 2E3)

  shp = gcore_brick (0, 0, 0)
  msh = PIPE ((0, 0, 0), (0, 0, 0.45), 0.125, 0.13, 3, 8, 2, 0, [0, 0, 0, 0])
  bod = BODY (solfec , kinem, shp, bulkmat, mesh = msh, form = formul, label = 'FU')
  if kinem != 'RIGID': bod.scheme = 'DEF_LIM'
  INITIAL_VELOCITY (bod, (0, 0, 0), (0, 0, 2*PI))
  bod.damping = damp

  RUN (solfec, solver, stop)

  return (solfec, kinem)

# main module
#import rpdb2; rpdb2.start_embedded_debugger('a')

pairs = []
pairs.append (scene_run (0.005, 1, 'FINITE_ELEMENT', 'BC', 0))
pairs.append (scene_run (0.005, 1, 'RIGID', 'BC', 0))
pairs.append (scene_run (0.001, 1, 'FINITE_ELEMENT', 'BC', 0))
pairs.append (scene_run (0.001, 1, 'RIGID', 'BC', 0))
pairs.append (scene_run (0.005, 1, 'FINITE_ELEMENT', 'BC', 1E-5))
pairs.append (scene_run (0.005, 1, 'FINITE_ELEMENT', 'BC', 1E-4))
pairs.append (scene_run (0.005, 1, 'FINITE_ELEMENT', 'BC', 1E-3))

if not VIEWER() and pairs [0][0].mode == 'READ':

  try:
    import matplotlib.pyplot as plt

    d = []
    for i in range (0, 2):
      solfec = pairs [i][0]
      b = BYLABEL (solfec, 'BODY', 'FU')
      p0 = (0.1315, 0, 0)
      p1 = (-0.1315, 0, 0)
      th = HISTORY (solfec, [(b, p0, 'DX'), (b, p0, 'DY'), (b, p1, 'DX'), (b, p1, 'DY')], 0, 1)

      d1 = []
      for j in range (0, len (th[0])):
	x1 = p0[0] + th [1][j]
	y1 = p0[1] + th [2][j]
	x2 = p1[0] + th [3][j]
	y2 = p1[1] + th [4][j]
	dx = x2 - x1
	dy = y2 - y1
	l1 = sqrt (dx*dx + dy*dy)
	d1.append (l1)

      d.append (d1)

    plt.title ('Step 0.005')
    plt.plot (th [0], d[0], label = 'CR-FE')
    plt.plot (th [0], d[1], label = 'RIG')
    plt.xlabel ('Time [s]')
    plt.ylabel ('Inner diamter [m]')
    plt.legend(loc = 'upper left')
    plt.savefig ('out/fuel-brick-rotate/fb_dd_step005_eta0.eps')
  except ImportError:
    pass # no reaction

  plt.clf ()

  try:
    import matplotlib.pyplot as plt

    for i in range (0, 2):
      solfec = pairs [i][0]
      b = BYLABEL (solfec, 'BODY', 'FU')
      p0 = (0.1315, 0, 0)
      th = HISTORY (solfec, [(b, p0, 'DY')], 0, 1)
      if i == 0: lab = 'DY: FE'
      elif i == 1: lab = 'DY: rigid'
      plt.plot (th [0], th [1], label = lab)

    plt.title ('Step 0.005')
    plt.xlabel ('Time [s]')
    plt.ylabel ('Displacement DY [m]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/fuel-brick-rotate/fb_dy_step005_eta0.eps')
  except ImportError:
    pass # no reaction

  plt.clf ()

  try:
    import matplotlib.pyplot as plt

    d = []
    for i in range (2, 4):
      solfec = pairs [i][0]
      b = BYLABEL (solfec, 'BODY', 'FU')
      p0 = (0.1315, 0, 0)
      p1 = (-0.1315, 0, 0)
      th = HISTORY (solfec, [(b, p0, 'DX'), (b, p0, 'DY'), (b, p1, 'DX'), (b, p1, 'DY')], 0, 1)

      d1 = []
      for j in range (0, len (th[0])):
	x1 = p0[0] + th [1][j]
	y1 = p0[1] + th [2][j]
	x2 = p1[0] + th [3][j]
	y2 = p1[1] + th [4][j]
	dx = x2 - x1
	dy = y2 - y1
	l1 = sqrt (dx*dx + dy*dy)
	d1.append (l1)

      d.append (d1)

    plt.title ('Step 0.001')
    plt.plot (th [0], d[0], label = 'CR-FE')
    plt.plot (th [0], d[1], label = 'RIG')
    plt.xlabel ('Time [s]')
    plt.ylabel ('Inner diamter [m]')
    plt.legend(loc = 'upper left')
    plt.savefig ('out/fuel-brick-rotate/fb_dd_step001_eta0.eps')
  except ImportError:
    pass # no reaction

  plt.clf ()

  try:
    import matplotlib.pyplot as plt

    for i in range (2, 4):
      solfec = pairs [i][0]
      b = BYLABEL (solfec, 'BODY', 'FU')
      p0 = (0.1315, 0, 0)
      th = HISTORY (solfec, [(b, p0, 'DY')], 0, 1)
      if i == 0: lab = 'DY: FE'
      elif i == 1: lab = 'DY: rigid'
      plt.plot (th [0], th [1], label = lab)

    plt.title ('Step 0.001')
    plt.xlabel ('Time [s]')
    plt.ylabel ('Displacement DY [m]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/fuel-brick-rotate/fb_dy_step001_eta0.eps')
  except ImportError:
    pass # no reaction

  plt.clf ()

  try:
    import matplotlib.pyplot as plt

    d = []
    for i in range (4, 7):
      solfec = pairs [i][0]
      b = BYLABEL (solfec, 'BODY', 'FU')
      p0 = (0.1315, 0, 0)
      p1 = (-0.1315, 0, 0)
      th = HISTORY (solfec, [(b, p0, 'DX'), (b, p0, 'DY'), (b, p1, 'DX'), (b, p1, 'DY')], 0, 1)

      d1 = []
      for j in range (0, len (th[0])):
	x1 = p0[0] + th [1][j]
	y1 = p0[1] + th [2][j]
	x2 = p1[0] + th [3][j]
	y2 = p1[1] + th [4][j]
	dx = x2 - x1
	dy = y2 - y1
	l1 = sqrt (dx*dx + dy*dy)
	d1.append (l1)

      d.append (d1)

    plt.title ('Step 0.005')
    plt.plot (th [0], d[0], label = 'CR-FE\eta=1E-4')
    plt.plot (th [0], d[1], label = 'CR-FE\eta=1E-3')
    plt.plot (th [0], d[2], label = 'CR-FE\eta=1E-2')
    plt.xlabel ('Time [s]')
    plt.ylabel ('Inner diamter [m]')
    plt.legend(loc = 'upper left')
    plt.savefig ('out/fuel-brick-rotate/fb_dd_step005_eta1.eps')
  except ImportError:
    pass # no reaction

  plt.clf ()

  try:
    import matplotlib.pyplot as plt

    for i in range (4, 7):
      solfec = pairs [i][0]
      b = BYLABEL (solfec, 'BODY', 'FU')
      p0 = (0.1315, 0, 0)
      th = HISTORY (solfec, [(b, p0, 'DY')], 0, 1)
      if i == 4: lab = 'CR-FE\eta=1E-4'
      elif i == 5: lab ='CR-FE\eta=1E-5'
      elif i == 6: lab ='CR-FE\eta=1E-6'
      plt.plot (th [0], th [1], label = lab)

    plt.title ('Step 0.001')
    plt.xlabel ('Time [s]')
    plt.ylabel ('Displacement DY [m]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/fuel-brick-rotate/fb_dy_step005_eta1.eps')
  except ImportError:
    pass # no reaction

  for pair in pairs:
    solfec = pair [0]
    kinem = pair [1]
    timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL']
    dur = DURATION (solfec)
    th = HISTORY (solfec, timers, dur[0], dur[1])
    total = 0.0

    for i in range (0, len (timers)):
      sum = 0.0
      for tt in th [i+1]: sum += tt
      print kinem, timers [i], 'TIME:', sum
      total += sum

    print kinem, 'TOTAL TIME:', total
