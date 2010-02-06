# masonry arch bridge example (MESH and RIGID)

from random import randint
from random import random
from random import seed
from math import sin
from math import cos

PI = 3.14159265358979323846 

def make_stone (x, y, z, r, kind, material, solfec):
  m = randint (8, 64)
  points = []
  for n in range (m):
    points.append (x + r * (1.0 - random()) * 2.0)
    points.append (y + r * (1.0 - random()) * 2.0)
    points.append (z + r * (1.0 - random()) * 2.0)

  hull = HULL (points, 2, 2)
  BODY (solfec, kind, hull, material)

def masonry_arch_create (base, radius, thickness, material, solfec, N, span):

  dalpha = PI / N
  width = 5.0
  edge = radius + width
  M = 4
  topup = radius * 0.3
  thick = thickness * 0.5

  # create two obstacle base bodies
  surfaces = [0, 0, 0, 0, 0, 0]
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
    TRANSLATE (hex, base)
    BODY (solfec, 'OBSTACLE', hex, material)
    alpha += (PI+dalpha)

  # create the remaining bricks
  surfaces = [1, 1, 1, 1, 1, 1]
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

    hex = HEX (nodes, 2, 3, 1, 1, surfaces)
    TRANSLATE (hex, base)
    BODY (solfec, 'RIGID', hex, material)

    # create side walls
    if i >= M and i < N - M:
      nodes = [xfarlo, -width/2., zfarlo,
	       xfarup, -width/2., zfarup,
	       xfarup, -width/2. + thick, zfarup,
	       xfarlo, -width/2. + thick, zfarlo,
	       xfarlo, -width/2., radius + topup,
	       xfarup, -width/2., radius + topup,
	       xfarup, -width/2. + thick, radius + topup,
	       xfarlo, -width/2. + thick, radius + topup]

      hex = HEX (nodes, 1, 2, 3, 1, surfaces)
      TRANSLATE (hex, base)
      BODY (solfec, 'RIGID', hex, material)

      nodes = [xfarlo, width/2. - thick, zfarlo,
	       xfarup, width/2. - thick, zfarup,
	       xfarup, width/2., zfarup,
	       xfarlo, width/2., zfarlo,
	       xfarlo, width/2. - thick, radius + topup,
	       xfarup, width/2. - thick, radius + topup,
	       xfarup, width/2., radius + topup,
	       xfarlo, width/2., radius + topup]

      hex = HEX (nodes, 1, 2, 3, 1, surfaces)
      TRANSLATE (hex, base)
      BODY (solfec, 'RIGID', hex, material)

    elif span < 1 and i > M and i < N-1:
      dxup = 2.0 * (radius + 0.5 * thickness - xfarup)
      dxlo = 2.0 * (radius + 0.5 * thickness - xfarlo)

      nodes = [xfarlo, -width/2., zfarlo,
	       xfarup, -width/2., zfarup,
	       xfarup, -width/2. + thick, zfarup,
	       xfarlo, -width/2. + thick, zfarlo,
	       xfarlo + dxlo, -width/2., zfarlo,
	       xfarup + dxup, -width/2., zfarup,
	       xfarup + dxup, -width/2. + thick, zfarup,
	       xfarlo + dxlo, -width/2. + thick, zfarlo]

      hex = HEX (nodes, 1, 2, 1, 1, surfaces)
      TRANSLATE (hex, base)
      BODY (solfec, 'RIGID', hex, material)

      nodes = [xfarlo, width/2. - thick, zfarlo,
	       xfarup, width/2. - thick, zfarup,
	       xfarup, width/2., zfarup,
	       xfarlo, width/2., zfarlo,
	       xfarlo + dxlo, width/2. - thick, zfarlo,
	       xfarup + dxup, width/2. - thick, zfarup,
	       xfarup + dxup, width/2., zfarup,
	       xfarlo + dxlo, width/2., zfarlo]

      hex = HEX (nodes, 1, 2, 1, 1, surfaces)
      TRANSLATE (hex, base)
      BODY (solfec, 'RIGID', hex, material)

    elif span == -1 and i < M and i > 0:

      x0 = -radius - 0.5 * thickness

      nodes = [xfarlo, -width/2., zfarlo,
	       xfarup, -width/2., zfarup,
	       xfarup, -width/2. + thick, zfarup,
	       xfarlo, -width/2. + thick, zfarlo,
	       x0, -width/2., zfarlo,
	       x0, -width/2., zfarup,
	       x0, -width/2. + thick, zfarup,
	       x0, -width/2. + thick, zfarlo]

      hex = HEX (nodes, 1, 2, 1, 1, surfaces)
      TRANSLATE (hex, base)
      BODY (solfec, 'RIGID', hex, material)

      nodes = [xfarlo, width/2. - thick, zfarlo,
	       xfarup, width/2. - thick, zfarup,
	       xfarup, width/2., zfarup,
	       xfarlo, width/2., zfarlo,
	       x0, width/2. - thick, zfarlo,
	       x0, width/2. - thick, zfarup,
	       x0, width/2., zfarup,
	       x0, width/2., zfarlo]

      hex = HEX (nodes, 1, 2, 1, 1, surfaces)
      TRANSLATE (hex, base)
      BODY (solfec, 'RIGID', hex, material)

    elif span == 1 and i >= N-M and i < N-1:

      x0 = radius + 0.5 * thickness

      nodes = [xfarlo, -width/2., zfarlo,
	       xfarup, -width/2., zfarup,
	       xfarup, -width/2. + thick, zfarup,
	       xfarlo, -width/2. + thick, zfarlo,
	       x0, -width/2., zfarlo,
	       x0, -width/2., zfarup,
	       x0, -width/2. + thick, zfarup,
	       x0, -width/2. + thick, zfarlo]

      hex = HEX (nodes, 1, 2, 1, 1, surfaces)
      TRANSLATE (hex, base)
      BODY (solfec, 'RIGID', hex, material)

      nodes = [xfarlo, width/2. - thick, zfarlo,
	       xfarup, width/2. - thick, zfarup,
	       xfarup, width/2., zfarup,
	       xfarlo, width/2., zfarlo,
	       x0, width/2. - thick, zfarlo,
	       x0, width/2. - thick, zfarup,
	       x0, width/2., zfarup,
	       x0, width/2., zfarlo]

      hex = HEX (nodes, 1, 2, 1, 1, surfaces)
      TRANSLATE (hex, base)
      BODY (solfec, 'RIGID', hex, material)

    alpha += dalpha

  # side rectangular fillings
  if span < 1:
    alpha = -0.5 * PI  + dalpha * (N - M - 1)
    xfarup =  (radius + thickness/2.0) * sin (alpha + dalpha)
    zfarup =  (radius + thickness/2.0) * cos (alpha + dalpha)

    dxup = 2.0 * (radius + 0.5 * thickness - xfarup)

    nodes = [xfarup, -width/2., zfarup,
	     xfarup + dxup, -width/2., zfarup,
	     xfarup + dxup, -width/2. + thick, zfarup,
	     xfarup, -width/2. + thick, zfarup,
             xfarup, -width/2., radius + topup,
	     xfarup + dxup, -width/2., radius + topup,
	     xfarup + dxup, -width/2. + thick, radius + topup,
	     xfarup, -width/2. + thick, radius + topup]

    hex = HEX (nodes, 1, 2, 1, 1, surfaces)
    TRANSLATE (hex, base)
    BODY (solfec, 'RIGID', hex, material)

    nodes = [xfarup, width/2. - thick, zfarup,
	     xfarup + dxup, width/2. - thick, zfarup,
	     xfarup + dxup, width/2., zfarup,
	     xfarup, width/2., zfarup,
             xfarup, width/2. - thick, radius + topup,
	     xfarup + dxup, width/2. - thick, radius + topup,
	     xfarup + dxup, width/2., radius + topup,
	     xfarup, width/2., radius + topup]

    hex = HEX (nodes, 1, 2, 1, 1, surfaces)
    TRANSLATE (hex, base)
    BODY (solfec, 'RIGID', hex, material)

  if span == -1 or span == 1:
    if span == -1: K = M - 1
    else: K = N - M - 1
    alpha = -0.5 * PI  + dalpha * K
    xfarup =  (radius + thickness/2.0) * sin (alpha + dalpha)
    zfarup =  (radius + thickness/2.0) * cos (alpha + dalpha)

    x0 = span * (radius + 0.5 * thickness)

    if span == 1:
      tmp = xfarup; xfarup = x0; x0 = tmp

    nodes = [x0, -width/2., zfarup,
	     xfarup, -width/2., zfarup,
	     xfarup, -width/2. + thick, zfarup,
	     x0, -width/2. + thick, zfarup,
             x0, -width/2., radius + topup,
	     xfarup, -width/2., radius + topup,
	     xfarup, -width/2. + thick, radius + topup,
	     x0, -width/2. + thick, radius + topup]

    hex = HEX (nodes, 1, 2, 1, 1, surfaces)
    TRANSLATE (hex, base)
    BODY (solfec, 'RIGID', hex, material)

    nodes = [x0, width/2. - thick, zfarup,
	     xfarup, width/2. - thick, zfarup,
	     xfarup, width/2., zfarup,
	     x0, width/2., zfarup,
             x0, width/2. - thick, radius + topup,
	     xfarup, width/2. - thick, radius + topup,
	     xfarup, width/2., radius + topup,
	     x0, width/2., radius + topup]

    hex = HEX (nodes, 1, 2, 1, 1, surfaces)
    TRANSLATE (hex, base)
    BODY (solfec, 'RIGID', hex, material)

    # end walls
    surfaces = [0, 0, 0, 0, 0, 0]
    if span == -1:
      nodes = [-radius - 0.5 * thickness, -width/2.0, 0.0,
	       -radius - 0.5 * thickness,  width/2.0, 0.0,
	       -radius - 0.5 * thickness - thick,  width/2.0, 0.0,
	       -radius - 0.5 * thickness - thick, -width/2.0, 0.0,
               -radius - 0.5 * thickness, -width/2.0, radius + topup,
	       -radius - 0.5 * thickness,  width/2.0, radius + topup,
	       -radius - 0.5 * thickness - thick,  width/2.0, radius + topup,
	       -radius - 0.5 * thickness - thick, -width/2.0, radius + topup]
    else:
      nodes = [radius + 0.5 * thickness + thick, -width/2.0, 0.0,
	       radius + 0.5 * thickness + thick,  width/2.0, 0.0,
	       radius + 0.5 * thickness,  width/2.0, 0.0,
	       radius + 0.5 * thickness, -width/2.0, 0.0,
               radius + 0.5 * thickness + thick, -width/2.0, radius + topup,
	       radius + 0.5 * thickness + thick,  width/2.0, radius + topup,
	       radius + 0.5 * thickness,  width/2.0, radius + topup,
	       radius + 0.5 * thickness, -width/2.0, radius + topup]

    hex = HEX (nodes, 1, 1, 1, 0, surfaces)
    TRANSLATE (hex, base)
    BODY (solfec, 'OBSTACLE', hex, material)

  NX = N
  NY = 6
  NZ = 10
  dx = (2.0 * radius + thickness) / NX
  dy = (width - 2.0 * thick) / NY
  dz = dy
  r  = min (dx, dy, dz) * 0.5

  # create gravel
  for i in range (0, NX):
    for j in range (0, NY):
      for k in range (0, NZ):
	x = -radius - 0.5 * thickness + dx * i + base [0]
	y = -0.5 * width + thick + dy * j + base [1]
	z = radius + topup + dz * k + base [2]
	make_stone (x, y, z, r, 'RIGID', material, solfec)

def masonry_bridge_create (material, solfec, N, M):

  ratio = 0.2
  radius = 10.0
  thickness = ratio * radius
  shift = 2.0 * radius + thickness

  for i in range (0, M):
    if i == 0: span = -1
    elif i == M-1: span = 1
    else: span = 0
    masonry_arch_create ((shift * i, 0, 0), radius, thickness, material, solfec, N, span)

### main module ###
#import rpdb2; rpdb2.start_embedded_debugger('a')

step = 0.001

seed (1)

solfec = SOLFEC ('DYNAMIC', step, 'out/bridge')

CONTACT_SPARSIFY (solfec, 0.005)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5)

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)

GRAVITY (solfec, (0, 0, -1), 9.81)

masonry_bridge_create (bulkmat, solfec, 27, 4)

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'CONTINUE', diagsolver = 'PROJECTED_GRADIENT')

OUTPUT (solfec, 50 * step)
RUN (solfec, gs, 10000 * step)

if not VIEWER() and solfec.mode == 'READ':

  timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL', 'PARBAL', 'GSINIT', 'GSRUN', 'GSCOM', 'GSMCOM']
  dur = DURATION (solfec)
  th = HISTORY (solfec, timers, dur[0], dur[1])
  total = 0.0

  for i in range (0, len(timers)):
    sum = 0.0
    for tt in th [i+1]: sum += tt
    print timers [i], 'TIME:', sum
    if i < 6: total += sum

  print 'TOTAL TIME:', total
