# Accuracy model with 2x fuel bricks, 2x interstitial bricks and 1x loose key

from solfec import *

from math import sin
from math import cos
from math import sqrt
from math import pi

PI = pi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def fuel_brick (
  surid = 0,
  inr = 0.12966,
  outd = 0.4455,
  keyinw = 0.03717,
  keyoutw1 = 0.04908,
  keyoutw3 = 0.04908,
  keyh = 0.03762,
  keylengthl = 0.3811,
  keylengths = 0.2749,
  height = 0.8893,
  basekeyoutr = 0.1867,
  baseringh = 0.0250):

  step = 2.0 * PI / 16.0
  angle = PI / 2.0 + step / 2.0

  faces = [4, 0, 1, 5, 4, surid,
         4, 1, 2, 6, 5, surid,
         4, 2, 3, 7, 6, surid,
         4, 3, 0, 4, 7, surid,
         4, 0, 3, 2, 1, surid,
         4, 4, 5, 6, 7, surid]

  zero = (0, 0, 0)
  zet = (0, 0, 1)
  toph =  outd * 0.5
  outr = toph / sin (angle)
  keyoutw2 = keyoutw1 + (keyoutw3 - keyoutw1) * (toph - inr) / (toph - inr - keyh)

  #define square blocks
  a = [keyoutw1/2.0, inr, 0,
       keyoutw2/2.0, toph-keyh, 0,
      -keyoutw2/2.0, toph-keyh, 0,
      -keyoutw1/2.0, inr, 0,
       keyoutw1/2.0, inr, height,
       keyoutw2/2.0, toph-keyh, height,
      -keyoutw2/2.0, toph-keyh, height,
      -keyoutw1/2.0, inr, height]

  list = []
  c1 = None
  c2 = None

  for i in range (8):

    cvy = CONVEX (a, faces, 0);
    ROTATE (cvy, zero, zet, i*45);

    if i == 1: c1 = cvy

    list.append (cvy)

  #define flat surfaces for key ways
  b = [keyoutw2/2.0, toph-keyh, 0,
       keyoutw3/2.0, toph, 0,
       keyinw/2.0, toph, 0,
       keyinw/2.0, toph-keyh, 0,
       keyoutw2/2.0, toph-keyh, height,
       keyoutw3/2.0, toph, height,
       keyinw/2.0, toph, height,
       keyinw/2.0, toph-keyh, height]

  c = [-keyinw/2.0, toph-keyh, 0,
       -keyinw/2.0, toph, 0,
       -keyoutw3/2.0, toph, 0,
       -keyoutw2/2.0, toph-keyh, 0,
       -keyinw/2.0, toph-keyh, height,
       -keyinw/2.0, toph, height,
       -keyoutw3/2.0, toph, height,
       -keyoutw2/2.0, toph-keyh, height]

  for i in range (8):
    cvx = CONVEX (b, faces, 0)
    cvy = CONVEX (c, faces, 0)

    ROTATE (cvx, zero, zet, i*45)
    ROTATE (cvy, zero, zet, i*45)

    if i == 1: c2 = cvx

    list.append (cvx)
    list.append (cvy)

  #define 6-sided bricks
  #define points of single brick
  d = [a [9], a [10], 0,
       c [6], c [7], 0,
       outr * cos (angle), outr * sin (angle), 0,
       outr * cos (angle+step), outr * sin (angle+step), 0,
       c2.vertex(1)[0], c2.vertex(1)[1], 0,
       c1.vertex(0)[0], c1.vertex(0)[1], 0,
       a [9], a [10], height,
       c [6], c [7], height,
       outr * cos (angle), outr * sin (angle), height,
       outr * cos (angle+step), outr * sin (angle+step), height,
       c2.vertex(1)[0], c2.vertex(1)[1], height,
       c1.vertex(0)[0], c1.vertex(0)[1], height]

  #define faces of shape
  g = [4, 0, 1, 7, 6, surid,
       4, 1, 2, 8, 7, surid,
       4, 2, 3, 9, 8, surid,
       4, 3, 4, 10, 9, surid,
       4, 4, 5, 11, 10, surid,
       4, 5, 0, 6, 11, surid,
       6, 0, 5, 4, 3, 2, 1, surid,
       6, 6, 7, 8, 9, 10, 11, surid]

  for i in range (8):

    cvy = CONVEX (d, g, 0)
    ROTATE (cvy, zero, zet, i*45);

    list.append (cvy)

  #define long key ways
  e = [keyinw/2.0, toph-keyh, 0,
       keyinw/2.0, toph, 0,
       -keyinw/2.0, toph, 0,
       -keyinw/2.0, toph-keyh, 0,
       keyinw/2.0, toph-keyh, height-keylengthl,
       keyinw/2.0, toph, height-keylengthl,
       -keyinw/2.0, toph, height-keylengthl,
       -keyinw/2.0, toph-keyh, height-keylengthl]

  #define short key ways
  h = [keyinw/2.0, toph-keyh, 0,
       keyinw/2.0, toph, 0,
       -keyinw/2.0, toph, 0,
       -keyinw/2.0, toph-keyh, 0,
       keyinw/2.0, toph-keyh, height-keylengths,
       keyinw/2.0, toph, height-keylengths,
       -keyinw/2.0, toph, height-keylengths,
       -keyinw/2.0, toph-keyh, height-keylengths]

  for i in range (8):
    if i % 2 != 0:
        cvy = CONVEX (e, faces, 0);
        ROTATE (cvy, zero, zet, i*45);
    elif i % 2 == 0:
        cvy= CONVEX (h, faces, 0);
        ROTATE (cvy, zero, zet, i*45);

    list.append (cvy)

  #define Base Ring
  nodes = [0, 0, 0,
	   1, 0, 0,
	   1, 1, 0,
	   0, 1, 0,
	   0, 0, 1,
	   1, 0, 1,
	   1, 1, 1,
	   0, 1, 1]

  msh = HEX (nodes, 1, 16, 1, 0, [surid, surid, surid, surid, surid, surid])
  SCALE (msh, (basekeyoutr-inr, 3.14159 * (360.0 / 180.0) * inr/cos(11.25*2*pi/360), baseringh))
  TRANSLATE (msh, (inr/cos(11.25*2*pi/360), 0, -baseringh))
  BEND (msh, (0, 0, 0), (0, 0, 1), 360)
  ROTATE (msh,(0,0,0),(0,0,1),11.25)
  list.append (msh)

  return list

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def loose_key(
    surid = 0,
    keyw = 0.03634,
    keyl = 0.0746,
    keyfh = 0.1158):

    vertices  = [1, 0, 0,
		 1, 1, 0,
		 0, 1, 0,
		 0, 0, 0,
		 1, 0, 1,
		 1, 1, 1,
		 0, 1, 1,
		 0, 0, 1]

    key = HEX (vertices, 8, 4, 2, 0, [surid, surid, surid, surid, surid, surid])

    SCALE (key, (keyw, keyl, keyfh))
    TRANSLATE (key, (-keyw/2.0, -keyl/2.0, 0))

    return key

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def interstitial_brick(
    surid = 0,
    iheight = 0.9217,
    ilengthil = 0.1112,
    ilengthis = 0.0461,
    ilengthol = 0.1861,
    ilengthos = 0.12424,
    keyawidth = 0.03596,
    keyadepth = 0.03141,
    keyaheight = 0.1518):

    intera = [ilengthis/2.0, ilengthil/2.0, 0,
              ilengthos/2.0, ilengthol/2.0, 0,
             -ilengthos/2.0, ilengthol/2.0, 0,
             -ilengthis/2.0, ilengthil/2.0, 0,
              ilengthis/2.0, ilengthil/2.0, iheight,
              ilengthos/2.0, ilengthol/2.0, iheight,
             -ilengthos/2.0, ilengthol/2.0, iheight,
             -ilengthis/2.0, ilengthil/2.0, iheight]

    faces = [4, 0, 1, 5, 4, surid,
           4, 1, 2, 6, 5, surid,
           4, 2, 3, 7, 6, surid,
           4, 3, 0, 4, 7, surid,
           4, 0, 3, 2, 1, surid,
           4, 4, 5, 6, 7, surid]

    list = []
    for i in range (4):
        cvy = CONVEX (intera, faces, 0);
        ROTATE (cvy, (0,0,0), (0,0,1), i*90);
        list.append (cvy)

    interb = [ilengthil/2.0, ilengthis/2.0, 0,
              ilengthol/2.0, ilengthos/2.0, 0,
              ilengthos/2.0, ilengthol/2.0, 0,
              ilengthis/2.0, ilengthil/2.0, 0,
              ilengthil/2.0, ilengthis/2.0, iheight,
              ilengthol/2.0, ilengthos/2.0, iheight,
              ilengthos/2.0, ilengthol/2.0, iheight,
              ilengthis/2.0, ilengthil/2.0, iheight]

    for i in range (4):
        cvy = CONVEX (interb, faces, 0);
        ROTATE (cvy, (0,0,0), (0,0,1), i*90);

        list.append (cvy)

    interkeya = [keyawidth/2.0, ilengthol/2.0, (iheight-keyaheight),
                 keyawidth/2.0, ((ilengthol/2.0)+keyadepth), (iheight-keyaheight),
                -keyawidth/2.0, ((ilengthol/2.0)+keyadepth), (iheight-keyaheight),
                -keyawidth/2.0, ilengthol/2.0, (iheight-keyaheight),
                 keyawidth/2.0, ilengthol/2.0, iheight,
                 keyawidth/2.0, ((ilengthol/2.0)+keyadepth), iheight,
                -keyawidth/2.0, ((ilengthol/2.0)+keyadepth), iheight,
                -keyawidth/2.0, ilengthol/2.0, iheight]

    for i in range (4):
	hex = HEX (interkeya, 6, 6, 2, 0, [surid, surid, surid, surid, surid, surid])
        ROTATE (hex, (0,0,0), (0,0,1), i*90);
        list.append (hex)

    # add small sphere on centreline of brick to provide a volume to attach a constraint to
    com = SPHERE((0,0,iheight/2),0.0001,0,surid)
    list.append(com)

    return list

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# main body

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# epsilon
GEOMETRIC_EPSILON (1E-6)

# geometry parameters
translation = 0.4589
half_trans = translation *0.5
brickheight = 0.9143
inr = 0.12966
icornery = 0.1097
icornerx = 0.0219
itransy = 0.1785
iheight = 0.9217

# control parameters
timestep = 1.0E-4
duration = 5          # simulation length, s

# setup analysis object
outpath = 'out/mbfcp/two_fuel_bricks'
solfec = SOLFEC ('DYNAMIC', timestep, outpath)

# materials
SURFACE_MATERIAL(solfec, model = 'SIGNORINI_COULOMB', friction = 0.1)
bulkmat = BULK_MATERIAL(solfec, label = 'GRAPHITE', model = 'KIRCHHOFF', young = 10.85E9, poisson = 0.21, density = 1810)

# Create geometry
fuels = []
inters =[]

for i in range(2):      # fuel bricks
    shp = fuel_brick(i)
    body_brick = BODY (solfec , 'RIGID', TRANSLATE(shp, (translation*i, 0, 0)), bulkmat, 'FUELBRICK%i' % i)
    fuels.append (body_brick)

for j in range (2):     # interstitial bricks
    shp_inter = interstitial_brick(2)
    ROTATE (shp_inter, (0,0,0), (0,0,1),45)
    body_inter = BODY (solfec, 'RIGID', TRANSLATE (shp_inter, (half_trans, half_trans-(translation*j), -itransy)), bulkmat, 'INTERBRICK%i' % j)
    inters.append (body_inter)

# loose key
shp_key = loose_key(4)
ROTATE (shp_key, (0,0,0), (0,0,1),90)
TRANSLATE (shp_key, (half_trans,0,0.6144))
loosekey = BODY (solfec , 'RIGID', shp_key, bulkmat, 'LOOSEKEY')

# Define Rigid links
PUT_RIGID_LINK (inters[0], inters[1], (half_trans, half_trans, (iheight/2)-itransy),
            (half_trans, -half_trans, (iheight/2)-itransy))


# Apply normal-to-plane boundary conditions
plane_normal = (0,0,1)
# LH  fuel brick
FIX_DIRECTION (fuels[0], (0, inr, 0), plane_normal)
FIX_DIRECTION (fuels[0], (0, -inr, 0), plane_normal)
FIX_DIRECTION (fuels[0], (-inr, 0, 0), plane_normal)
#RH fuel brick
FIX_DIRECTION (fuels[1], (translation, inr, 0), plane_normal)
FIX_DIRECTION (fuels[1], (translation, -inr, 0), plane_normal)
FIX_DIRECTION (fuels[1], (translation+inr, 0, 0), plane_normal)
# North interstitial brick
FIX_DIRECTION (inters[0], (half_trans, half_trans+icornery, 0), plane_normal)
FIX_DIRECTION (inters[0], (half_trans, half_trans-icornery, 0), plane_normal)
FIX_DIRECTION (inters[0], (half_trans-icornery, half_trans, 0), plane_normal)
# South interstitial brick
FIX_DIRECTION (inters[1], (half_trans, -half_trans+icornery, 0), plane_normal)
FIX_DIRECTION (inters[1], (half_trans, -half_trans-icornery, 0), plane_normal)
FIX_DIRECTION (inters[1], (half_trans-icornery, -half_trans, 0), plane_normal)
# Loose key
FIX_DIRECTION (loosekey, (half_trans, 0.03634*0.5, 0.6144), plane_normal)
FIX_DIRECTION (loosekey, (half_trans, -0.03634*0.5, 0.6144), plane_normal)
FIX_DIRECTION (loosekey, ((half_trans-(0.0746*0.2)), 0, 0.6144), plane_normal)

# Apply BCs normal to input motion to  LH fuel brick
FIX_DIRECTION (fuels[0], (0, inr, 0), (1,0,0))
FIX_DIRECTION (fuels[0], (0, -inr, 0), (1,0,0))

# Derive and appy input time history
numPeriod = 0.666
numDivision = 1000
numCycles = 7.5
mag = 0.25
sine_velocity = []
for i in range (int(numDivision*numCycles)):
    x = numPeriod*i/float(numDivision)
    y = mag*sin((i/float(numDivision))*2*pi)
    sine_velocity.append (x)
    sine_velocity.append (y)
motionHistory = TIME_SERIES (sine_velocity)
SET_VELOCITY (fuels[0], (0,inr,0.5), (0,1,0), motionHistory)

# Run analysis
OUTPUT(solfec, 1E-2)
gs = GAUSS_SEIDEL_SOLVER (1E-5, 100)
MBFCP_EXPORT (solfec, outpath + '/two_fuel_bricks.mbfcp')
RUN (solfec, gs, duration)

# Post-processing
def write_data (t, v, path):
  out = open (path, 'w')
  n = min (len(t), len(v))
  for i in range (n):
    out.write ('%g\t%g\n' % (t[i], v[i]))
  out.close ()

if not VIEWER() and solfec.mode == 'READ':

  dur = DURATION (solfec)
  th = HISTORY (solfec, [
	(fuels [1], fuels [1].center, 'DY'),
	(fuels [1], fuels [1].center, 'DX'),
        (loosekey, loosekey.center, 'DY'),
        (loosekey, loosekey.center, 'DX')
       ], dur[0], dur [1])

  write_data (th[0], th[1], outpath + '/Solfec_1E-4_fuel_2_dy.dat')
  write_data (th[0], th[3], outpath + '/Solfec_1E-4_loose_dy.dat')
