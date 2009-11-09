# arch example

from math import sin
from math import cos

N = 27
PI = 3.14159265358979323846 

def masonry_arch_create (ratio, material, solfec):

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

### main module ###

step = 0.001

solfec = SOLFEC ('DYNAMIC', step, 'out/arch')

CONTACT_SPARSIFY (solfec, 0.005)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.4)

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1, poisson = 0, density = 1)

GRAVITY (solfec, (0, 0, -1), 9.81)

#import rpdb2; rpdb2.start_embedded_debugger('a')

masonry_arch_create (0.1070, bulkmat, solfec)

def gscallback (gs):
  print gs.error
  return 0

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'CALLBACK', callback = gscallback, diagsolver = 'PROJECTED_GRADIENT')

RUN (solfec, gs, 1000 * step)
