# two-body hybrid modeling example 
import matplotlib.pyplot as plt
import time, sys, os
sys.path.append(os.path.dirname(__file__))
from acc_sweep import *

step = 1E-5  # time step
stop = 5.0   # duration of the simulation
damp = 1E-5  # amount of stiffness proportional damping
lofq = 5     # low frequency for the sweep
hifq = 15    # high frequency for the sweep
amag = 10.0  # acceleration magnitude
nbodl = 2    # number of bodies
nele = 2     # number of elements per body (along x, y, z)
nmod = 6     # number of linear modes to produce
l = 0.1      # length of one body
w = 0.1      # widhth of one body
h = 0.1      # height of one body
gap = 0.002  # gap
ostep = 1E-3 # output step
wavg = 0.01  # energy averaging time window [t-wavg/2, t+wavg/2]
fstop = 7.0  # end frequency for averaging velocities
             # (hand tunded in order to cover the pre-drop area)

GEOMETRIC_EPSILON (1E-4*gap) # geometrical tolerance << gap

solfec = SOLFEC ('DYNAMIC', step, 'out/hybrid_modeling/two_cubes_nscd')
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.05, restitution = 0.0)
bulk = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E9, poisson = 0.25, density = 1E3)

nodes = [0, 0, 0,
	 w, 0, 0,
	 w, l, 0,
	 0, l, 0,
	 0, 0, h,
	 w, 0, h,
	 w, l, h,
	 0, l, h]

bodies = [] # empty list for FEM bodies
for j in range (0, nbodl):
  shape = HEX (nodes, nele, nele, nele, 1, [1]*6)
  TRANSLATE (shape, (0, j*(l+gap), 0))
  body = BODY (solfec, 'FINITE_ELEMENT', shape, bulk, form = 'BC')
  body.scheme = 'DEF_LIM' # implicit time integration
  body.damping = damp
  bodies.append (body) # append list

# base and side walls
s0 = HEX (nodes, 1, 1, 1, 1, [1]*6)
SCALE (s0, (1, ((nbodl+2)*l + (nbodl+1)*gap)/l, 0.1))
TRANSLATE (s0, (0, -gap-l, -0.1*h))
s1 = COPY (s0)
s2 = COPY (s0)
ROTATE (s1, (w/2,0,h/2), (0, 1, 0), +90)
ROTATE (s2, (w/2,0,h/2), (0, 1, 0), -90)
TRANSLATE (s1, (-gap, 0, 0))
TRANSLATE (s2, (+gap, 0, 0))

# moving obstacle
s3 = HEX (nodes, 1, 1, 1, 1, [1]*6)
TRANSLATE (s3, (0, -gap-l, 0))
s4 = HEX (nodes, 1, 1, 1, 1, [1]*6)
TRANSLATE (s4, (0, nbodl*(l+gap), 0))
body = BODY (solfec, 'OBSTACLE', [s0, s1, s2, s3, s4], bulk)
FIX_DIRECTION (body, (w,-(gap+l),0), (1, 0, 0))
FIX_DIRECTION (body, (w,-(gap+l),0), (0, 0, 1))
FIX_DIRECTION (body, (w,-gap,0), (1, 0, 0))
FIX_DIRECTION (body, (w,-gap,h), (1, 0, 0))
FIX_DIRECTION (body, (w,-gap,0), (0, 0, 1))

#enable gravity
GRAVITY (solfec, (0, 0, -10))

# apply velocity constraint corresponding to the acceleration sin sweep
if VIEWER():
  data = acc_sweep (step, stop, lofq, hifq, amag, None, None, None)
else:
  data = acc_sweep (step, stop, lofq, hifq, amag,
                  'out/hybrid_modeling/acc.png',
		  'out/hybrid_modeling/vel.png',
		  'out/hybrid_modeling/dsp.png')
SET_VELOCITY (body, (w/2.,-(gap+l)/2.,h/2.), (0, 1, 0), TIME_SERIES (data))

# create constraints solver
slv = NEWTON_SOLVER ()

# output results every 'ostep'
OUTPUT (solfec, ostep)

# run simulation
t0 = time.time()
RUN (solfec, slv, stop)
if solfec.mode == 'WRITE':
  print "Analysis run time:", (time.time() - t0)/60.0, "minutes"

# post-process results
if not VIEWER() and solfec.mode == 'READ':
  iavg = 1 + int (wavg / ostep) / 2
  tstop = 0.0
  data = []
  for b in bodies:
    data.append ((b, 'KINETIC'))
    data.append ((b, 'INTERNAL'))
    data.append ((b, b.center, 'VY'))
  th = HISTORY (solfec, data, 0, stop)
  n = len (th[0])
  for k in range (0, nbodl):
    fq = []
    ek = []
    ei = []
    vy = []
    for i in range (0, n):
      if i >= iavg and i < n-iavg-1:
	vek = 0.0
	vei = 0.0
	vvy = 0.0
	for j in range (i-iavg, i+iavg+1):
	  vek += th [3*k+1][j]
	  vei += th [3*k+2][j]
	  vvy += abs(th [3*k+3][j])

        f = lofq + (hifq-lofq)*(th[0][i]/stop)
	if f > fstop: tstop = th[0][i]
	fq.append (f)
	ek.append (vek/(2.0*iavg+1.0))
	ei.append (vei/(2.0*iavg+1.0))
	vy.append (vvy/(2.0*iavg+1.0))
       
    plt.clf ()
    plt.plot (fq, ek, lw = 2, label = 'kinetic')
    plt.plot (fq, ei, lw = 2, label = 'internal')
    plt.xlim ((lofq, hifq))
    plt.legend(loc = 'best')
    plt.xlabel ('Frequency $(Hz)$')
    plt.ylabel ('Energy $(J)$')
    plt.savefig ('out/hybrid_modeling/two_cubes_nscd/ene'+str(k)+'.png')

    plt.clf ()
    plt.plot (fq, vy, lw = 2)
    plt.xlim ((lofq, hifq))
    plt.xlabel ('Frequency $(Hz)$')
    plt.ylabel ('Velocity vy $(m/s)$')
    plt.savefig ('out/hybrid_modeling/two_cubes_nscd/vy'+str(k)+'.png')

    # averge pre-drop-off velocity for body k
    vavg = 0.0
    nvavg = 0.0
    for (f, v) in zip(fq, vy):
      if f < fstop:
        vavg += v
	nvavg += 1.0

    print 'Average pre-drop-off velocity for body', k, 'is', vavg/nvavg
  
  # average input impact velocity
  vavg = 0.0
  nvavg = 0.0
  SEEK (solfec, 0.0)
  while solfec.time < tstop:
    for con in solfec.constraints:
      if con.kind == 'CONTACT':
        vavg += con.V[2]
	nvavg += 1.0
    FORWARD (solfec, 1)

  print 'Avererage impact input velocity:', vavg/nvavg
