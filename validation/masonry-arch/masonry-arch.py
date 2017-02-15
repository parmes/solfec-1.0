from math import sin, cos
from sys import stdout

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

def progress_callback (solfec, stop):
  print '\b\b\b\b\b%2d %%' % (50.0 + 50.0 * solfec.time / stop) , 
  stdout.flush ()
  return 1

def analysis (step, stop, frict, ratio):
  solfec = SOLFEC ('DYNAMIC', step, 'out/masonry-arch/fr-%g-ratio-%g' % (fr, ra))
  CALLBACK (solfec, step, (solfec, stop), progress_callback)
  if not VIEWER(): solfec.verbose = 'OFF'
  CONTACT_SPARSIFY (solfec, 0.005)
  SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = frict)
  bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1, poisson = 0, density = 1)
  GRAVITY (solfec, (0, 0, -9.81))
  masonry_arch_create (ratio, bulkmat, solfec)
  gs = GAUSS_SEIDEL_SOLVER (1E-4, 10000) # XXX
  if solfec.mode == 'WRITE': RUN (solfec, gs, stop)
  else: SEEK (solfec, stop)
  en = ENERGY(solfec)
  return en[0]

try:
  import numpy as np
except (ImportError, RuntimeError):
  import sys
  print "Unexpected error:", sys.exc_info()[1]
  print "Importing numpy has failed!"
  sys.exit (1)

# main module
GEOMETRIC_EPSILON (1E-4)
step = 0.001
stop = 0.1 # XXX
frict = [0.2, 0.25, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.5, 0.6, 0.7, 0.8]
ratio = [0.05, 0.07, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.23, 0.25]
X, Y = np.meshgrid (frict, ratio)
Z = np.zeros(X.shape)
n = len(frict)*len(ratio)
k = i = 0
for fr in frict:
  j = 0
  for ra in ratio:
    print 'fr-%.2f-ra-%.2f' % (fr, ra), '[' + '=' * (k/2) + ' ' * ((n/2)-(k/2)-1) + ']', '    ' , 
    Z[j][i] = analysis (step, stop, fr, ra)
    print
    j = j + 1
    k = k + 1
  i = i + 1

if not VIEWER():
  try:
    import matplotlib.pyplot as plt
    plt.clf()
    plt.figure()
    cp = plt.contour(X, Y, Z, levels = [0, 1E-3, 1E-1, 1])
    plt.clabel(cp, inline=True, fontsize=10)
    plt.axis('scaled')
    plt.title('Solfec - kinetic energy at %gs' % stop)
    plt.xlabel('Coefficient of friction')
    plt.ylabel('h/r')
    plt.savefig ('validation/masonry-arch/masonry-arch-energy.png', bbox_inches='tight', pad_inches=0)
    plt.clf()
    plt.figure()
    cp = plt.contour(X, Y, Z, levels = [0, 1E-3, 1E-1, 1], linewidths=3)
    plt.clabel(cp, inline=True, fontsize=10)
    plt.axis('scaled')
    plt.xticks([],[])
    plt.yticks([],[])
    plt.savefig ('validation/masonry-arch/masonry-arch-energy-trans.png', transparent=True, bbox_inches='tight', pad_inches=0)
  except (ImportError, RuntimeError):
    import sys
    print "Unexpected error:", sys.exc_info()[1]
    print "Plotting has failed!"
    pass
