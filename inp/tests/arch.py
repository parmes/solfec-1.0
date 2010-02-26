# masonry arch test
from math import sin
from math import cos
from sys import stdout

step = 0.001
stepnum = 750
stop = step * stepnum
T = [] # plots
E1 = []
E2 = []

def masonry_arch_create (ratio, material, solfec):

  PI = 3.14159265358979323846 
  N = 27
  dalpha = PI / N
  radius = 10.0
  thickness = ratio * radius
  width = 5.0
  edge = radius + width

  dx = [0.05, 1.9, 0.05]

  surfaces = [0, 0, 0, 0, 0, 0]

  # create two obstacle base bodies
  alpha = -0.5*PI-dalpha 
  for i in range (2):

    xfarlo =  (radius + thickness/2.0) * sin (alpha)
    xnearlo = (radius - thickness/2.0) * sin (alpha)
    xfarup =  (radius + thickness/2.0) * sin (alpha + dalpha)
    xnearup = (radius - thickness/2.0) * sin (alpha + dalpha)
    zfarlo =  (radius + thickness/2.0) * cos (alpha)
    znearlo = (radius - thickness/2.0) * cos (alpha)
    zfarup =  (radius + thickness/2.0) * cos (alpha + dalpha)
    znearup = (radius - thickness/2.0) * cos (alpha + dalpha)

    nodes = [xfarlo, -width/2., zfarlo,
             xnearlo, -width/2., znearlo,
             xnearlo,  width/2., znearlo,
             xfarlo,  width/2., zfarlo,
             xfarup, -width/2.,  zfarup,
             xnearup, -width/2.,  znearup,
             xnearup,  width/2.,  znearup,
             xfarup,  width/2.,  zfarup]

    hex = HEX (nodes, 1, 1, 1, 0, surfaces)
    BODY (solfec, 'OBSTACLE', hex, material)
    alpha += (PI+dalpha)

  # create the remaining bricks
  alpha = -0.5*PI
  for i in range (N):
    xfarlo =  (radius + thickness/2.0) * sin (alpha)
    xnearlo = (radius - thickness/2.0) * sin (alpha)
    xfarup =  (radius + thickness/2.0) * sin (alpha + dalpha)
    xnearup = (radius - thickness/2.0) * sin (alpha + dalpha)
    zfarlo =  (radius + thickness/2.0) * cos (alpha)
    znearlo = (radius - thickness/2.0) * cos (alpha)
    zfarup =  (radius + thickness/2.0) * cos (alpha + dalpha)
    znearup = (radius - thickness/2.0) * cos (alpha + dalpha)

    nodes = [xfarlo, -width/2., zfarlo,
             xnearlo, -width/2., znearlo,
             xnearlo,  width/2., znearlo,
             xfarlo,  width/2., zfarlo,
             xfarup, -width/2.,  zfarup,
             xnearup, -width/2.,  znearup,
             xnearup,  width/2.,  znearup,
             xfarup,  width/2.,  zfarup]

    k = 1 + i % 4
    surfaces = [k, k, k, k, k, k]

    hex = HEX (nodes, 3, 2, 1, 0, surfaces, dx)
    BODY (solfec, 'RIGID', hex, material)
    alpha += dalpha


# first run progress
def progress_callback_one (solfec):
  print '\b\b\b\b\b%2d %%' % (50.0 * solfec.time / stop) , 
  stdout.flush ()
  T.append (solfec.time)
  E1.append (ENERGY (solfec) [0])
  return 1

# second run progress
def progress_callback_two (solfec):
  print '\b\b\b\b\b%2d %%' % (50.0 + 50.0 * solfec.time / stop) , 
  stdout.flush ()
  E2.append (ENERGY (solfec) [0])
  return 1

### main module ###
#import rpdb2; rpdb2.start_embedded_debugger('a')

gs = GAUSS_SEIDEL_SOLVER (1E-4, 10000)

# solfec one
solfec1 = SOLFEC ('DYNAMIC', step, 'out/tests/arch/one')
solfec1.verbose = 'OFF'
CONTACT_SPARSIFY (solfec1, 0.005)
SURFACE_MATERIAL (solfec1, model = 'SIGNORINI_COULOMB', friction = 0.5)
bulkmat1 = BULK_MATERIAL (solfec1, model = 'KIRCHHOFF', young = 1, poisson = 0, density = 1)
GRAVITY (solfec1, (0, 0, -9.81))
masonry_arch_create (0.1095, bulkmat1, solfec1)
CALLBACK (solfec1, step, solfec1, progress_callback_one)

# solfec two
solfec2 = SOLFEC ('DYNAMIC', step, 'out/tests/arch/two')
solfec2.verbose = 'OFF'
CONTACT_SPARSIFY (solfec2, 0.005)
SURFACE_MATERIAL (solfec2, model = 'SIGNORINI_COULOMB', friction = 0.5)
bulkmat2 = BULK_MATERIAL (solfec2, model = 'KIRCHHOFF', young = 1, poisson = 0, density = 1)
GRAVITY (solfec2, (0, 0, -9.81))
masonry_arch_create (0.1094, bulkmat2, solfec2)
CALLBACK (solfec2, step, solfec2, progress_callback_two)

print '    ' , 
RUN (solfec1, gs, stop)
RUN (solfec2, gs, stop)
print '\b\b\b\b\b\b\b' ,

if not VIEWER() and solfec1.mode == 'WRITE' and solfec2.mode == 'WRITE':

    if E1 [stepnum - 1] < 1E-4 and E2 [stepnum - 1] > 1E-3: print 'PASSED'
    else:
      print 'FAILED'
      print '(', 'Kinetic energy out of bounds: E1 = %g, E2 = %g' % (E1[stepnum-1], E2[stepnum-1]), ')'

    try:
      import matplotlib.pyplot as plt
      plt.clf ()
      plt.plot (T, E1, label='h/r = 0.1095')
      plt.plot (T, E2, label='h/r = 0.1094')
      plt.axis (xmin = 0, xmax = stop, ymin = -0.0005, ymax = 0.001)
      plt.xlabel ('Time [s]')
      plt.ylabel ('Kinetic energy [J]')
      plt.legend(loc = 'upper right')
      plt.savefig ('out/tests/arch/arch.eps')
    except ImportError:
      pass # no reaction
