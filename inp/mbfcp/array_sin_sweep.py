# array of bricks
import sys
import commands
sys.path.append ('inp/mesh')
from abaqusread import *
from math import cos 

dbgplt = 0 # enable debug plots

PI = 3.14159265358979323846

def gen_acc_sin_sweep (step):

  f1 = 1.0
  f2 = 15.0
  Rlin = 0.3
  A = 1.0
  t2 = (f2 - f1) / Rlin
  t = 0.0
  acc = []
  th = []
  ah = []
  c1 = PI*(f2-f1)/t2
  c2 = 2*PI*f1
  step *= 0.25 # since acceleration will need to be integrated into
               # the velocity desner step will improve efficiency

  while t < t2:
    a = A * cos (c1*t*t + c2*t) # analytical integration into velocity form produces a too long expression
    acc.append (t)
    acc.append (a)
    if dbgplt:
      th.append (t)
      ah.append (a)
    t = t + step

  if dbgplt:
    try:
      import matplotlib.pyplot as plt
      plt.plot (th, vh)
      plt.xlabel ('Time')
      plt.ylabel ('Velocity')
      plt.savefig ('out/mbfcp/array_sin_sweep/acc.eps')
    except: pass

  return (t2, acc)

#import rpdb2; rpdb2.start_embedded_debugger ('a')

step = 1E-3

solfec = SOLFEC ('DYNAMIC', step, 'out/mbfcp/array_sin_sweep')

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.1)

if RANK () == 0:
  commands.getoutput ("bunzip2 inp/mesh/A1.inp.bz2") # only first CPU unpacks the input

BARRIER () # let all CPUs meet here

ABAQUS_READ ('inp/mesh/A1.inp', solfec)

if RANK () == 0:
  commands.getoutput ("bzip2 inp/mesh/A1.inp") # only first CPU pack the input

gen = gen_acc_sin_sweep (step)
stop = gen [0]
acc = TIME_SERIES (gen [1])

# boundary conditions
for b in solfec.bodies:
  c = b.center
  t = b.tensor
  if b.label == 'FBH':
    p = TRANSLATE (c, (0.1, -0.1, 0))
    FIX_DIRECTION (b, p, (0, 0, -1))
    p = TRANSLATE (c, (-0.1, -0.1, 0))
    FIX_DIRECTION (b, p, (0, 0, -1))
    p = TRANSLATE (c, (0, 0.15, 0))
    FIX_DIRECTION (b, p, (0, 0, -1))
  elif b.label == 'IBH':
    if abs (t [0]-t[4]) < 1E-6: # symmetrical case
      p = TRANSLATE (c, (0.08, -0.08, 0))
      FIX_DIRECTION (b, p, (0, 0, -1))
      p = TRANSLATE (c, (-0.08, -0.08, 0))
      FIX_DIRECTION (b, p, (0, 0, -1))
      p = TRANSLATE (c, (0, 0.1, 0))
      FIX_DIRECTION (b, p, (0, 0, -1))
    else: # half keys
      if c [0] < 0.0:
	p = TRANSLATE (c, (0.05, 0, 0))
	FIX_DIRECTION (b, p, (0, 0, -1))
	p = TRANSLATE (c, (-0.04, 0.08, 0))
	FIX_DIRECTION (b, p, (0, 0, -1))
	FIX_DIRECTION (b, p, (-1, 0, 0))
	p = TRANSLATE (c, (-0.04, -0.08, 0))
	FIX_DIRECTION (b, p, (0, 0, -1))
	FIX_DIRECTION (b, p, (-1, 0, 0))
      elif c [0] > 1.0:
	p = TRANSLATE (c, (-0.05, 0, 0))
	FIX_DIRECTION (b, p, (0, 0, -1))
	p = TRANSLATE (c, (0.04, 0.08, 0))
	FIX_DIRECTION (b, p, (0, 0, -1))
	FIX_DIRECTION (b, p, (1, 0, 0))
	p = TRANSLATE (c, (0.04, -0.08, 0))
	FIX_DIRECTION (b, p, (0, 0, -1))
	FIX_DIRECTION (b, p, (1, 0, 0))
	SET_ACCELERATION (b, p, (0, 1, 0), acc)
      elif c [1] < 0.0:
	p = TRANSLATE (c, (0, 0.05, 0))
	FIX_DIRECTION (b, p, (0, 0, -1))
	p = TRANSLATE (c, (0.08, -0.04, 0))
	FIX_DIRECTION (b, p, (0, 0, -1))
	p = TRANSLATE (c, (-0.08, -0.04, 0))
	FIX_DIRECTION (b, p, (0, 0, -1))
      else:
	p = TRANSLATE (c, (0, -0.05, 0))
	FIX_DIRECTION (b, p, (0, 0, -1))
	p = TRANSLATE (c, (0.08, 0.04, 0))
	FIX_DIRECTION (b, p, (0, 0, -1))
	p = TRANSLATE (c, (-0.08, 0.04, 0))
	FIX_DIRECTION (b, p, (0, 0, -1))
  elif b.label == 'LKH':
    if t [0] < 0.51 * t [4]: # key along x
      p = TRANSLATE (c, (0.0369, 0.016, 0))
      FIX_DIRECTION (b, p, (0, 0, -1))
      p = TRANSLATE (c, (-0.0369, 0.016, 0))
      FIX_DIRECTION (b, p, (0, 0, -1))
      p = TRANSLATE (c, (0, -0.01785, 0))
      FIX_DIRECTION (b, p, (0, 0, -1))
    elif t [4] < 0.51 * t [0]: # key along y
      p = TRANSLATE (c, (0.016, 0.0369, 0))
      FIX_DIRECTION (b, p, (0, 0, -1))
      p = TRANSLATE (c, (0.016, -0.0369, 0))
      FIX_DIRECTION (b, p, (0, 0, -1))
      p = TRANSLATE (c, (-0.01785, 0, 0))
      FIX_DIRECTION (b, p, (0, 0, -1))
    else: # rectangular keys
      p = TRANSLATE (c, (0.016, 0.016, 0))
      FIX_DIRECTION (b, p, (0, 0, -1))
      p = TRANSLATE (c, (-0.016, 0.016, 0))
      FIX_DIRECTION (b, p, (0, 0, -1))
      if c [0] < 0.0: FIX_DIRECTION (b, p, (-1, 0, 0))
      p = TRANSLATE (c, (0, -0.01785, 0))
      if c [0] < 0.0: FIX_DIRECTION (b, p, (-1, 0, 0))
      FIX_DIRECTION (b, p, (0, 0, -1))


# solver and run

slv  = GAUSS_SEIDEL_SOLVER (1E-3, 300, 1E-6)

if RANK() == 0 and solfec.mode == 'WRITE':
  print 'Running', stop, 'seconds or analysis with step', step, '...'

RUN (solfec, slv, stop)
