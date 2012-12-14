# ------------------ #
# sine sweep example #
# ------------------ #
import matplotlib.pyplot as plt
from math import sin, cos, pi

step = 1E-4  # time step
stop = 5.0   # duration of the simulation
damp = 1E-6  # amount of stiffness proportional damping
lofq = 5     # low frequency for the sweep
hifq = 10    # high frequency for the sweep
amag = 10.0  # acceleration magnitude
nbodl = 10   # number of bodies
nele = 10    # number of elements per body (along y)
nmod = 6     # number of linear modes to produce
l = 0.1      # length of one body
w = 0.1      # widhth of one body
h = 0.1      # height of one body
gap = 0.001  # gap

GEOMETRIC_EPSILON (1E-9) # tiny geometrical tolerance (<< gap)

solfec = SOLFEC ('DYNAMIC', step, 'out/sweep')
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = 0.0)
bulk = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3) # graphite

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
  shape = HEX (nodes, 1, nele, 1, 1, [1]*6)
  TRANSLATE (shape, (0, j*(l+gap), 0))
  body = BODY (solfec, 'FINITE_ELEMENT', shape, bulk, form = 'BC')
  body.scheme = 'DEF_LIM' # implicit time integration
  body.damping = damp
  FIX_DIRECTION (body, body.center, (1, 0, 0))
  FIX_DIRECTION (body, body.center, (0, 0, 1))
  FIX_DIRECTION (body, (w, j*(l+gap), 0), (1, 0, 0))
  FIX_DIRECTION (body, (w, j*(l+gap), 0), (0, 0, 1))
  FIX_DIRECTION (body, (w, j*(l+gap)+l, 0), (1, 0, 0))
  bodies.append (body) # append list

# compute low eigen modes of the first body and get corresponding eigenfrequencies
eigf = []
out = MODAL_ANALYSIS (bodies [0], 6+nmod, 'out/sweep/modes', abstol = 1E-14, verbose = 'ON')
for o in out [0][6:6+nmod]:
  eigf.append (o**0.5) # eigen frequency = sqrt (eigen value)
print 'Single body eigen frequencies are:', eigf

# moving obstacle
s1 = HEX (nodes, 1, 1, 1, 1, [1]*6)
TRANSLATE (s1, (0, -gap-l, 0))
s2 = HEX (nodes, 1, 1, 1, 1, [1]*6)
TRANSLATE (s2, (0, nbodl*(l+gap), 0))
body = BODY (solfec, 'RIGID', [s1, s2], bulk)
FIX_DIRECTION (body, (w,-(gap+l),0), (1, 0, 0))
FIX_DIRECTION (body, (w,-(gap+l),0), (0, 0, 1))
FIX_DIRECTION (body, (w,-gap,0), (1, 0, 0))
FIX_DIRECTION (body, (w,-gap,h), (1, 0, 0))
FIX_DIRECTION (body, (w,-gap,0), (0, 0, 1))

# acceleration sweep
t = 0.0
v = 0.0
vt = []
va = []
vv = []
while t < stop:
  x = t + step/2. # mid-step time
  a = amag * sin (2*pi*(lofq+(hifq-lofq)*x/stop)*x) # mid-step acceleration
  v = v + a * step # mid-step integration of dv / dt = a into v
  vt.append (t)
  va.append (a)
  vv.append (v)
  t += step

# velocity is now all positive, but we want it centered about zero
i = 0
while vv[i+1] > vv[i]: i += 1 # find first maximum
j = len(vv)-1
while vv[j-1] < vv[j]: j -= 1 # find last minimum (if any)
k = j
while vv[k-1] > vv[k]: k -= 1 # and previous maximum
j = k
while vv[j-1] < vv[j]: j -= 1 # and previous minimum
# y = Ax + B
A = (0.5*(vv[j]+vv[k])-0.5*vv[i])/(vt[j]-vt[i])
B = 0.5*vv[i] - A*vt[i]
# shift down be the line connecting middles of end extrema
for i in range (0, len (vv)): vv[i] -= A*vt[i]+B

# after velocity has been shifted down, produce displacement envelope
vd = []
d = 0.0
for v in vv:
 d = d + v * step  # integration of dd / dt = v
 vd.append (d)

if not VIEWER ():
  plt.clf ()
  plt.plot (vt, va)
  plt.xlabel ('time $(s)$')
  plt.ylabel ('acceleration $(m/s^2)$')
  plt.savefig ('out/sweep/acc.svg')
  plt.clf ()
  plt.plot (vt, vv)
  plt.xlabel ('time $(s)$')
  plt.ylabel ('velocity $(m/s)$')
  plt.savefig ('out/sweep/vel.svg')
  plt.clf ()
  plt.plot (vt, vd)
  plt.xlabel ('time $(s)$')
  plt.ylabel ('displacement $(m)$')
  plt.savefig ('out/sweep/dsp.svg')

# apply velocity constraint corresponding to the acceleration sin sweep
data = []
for (t, v) in zip (vt, vv): data += [t, v] # [t0,v0,t1,v1,...]
SET_VELOCITY (body, (w/2.,-(gap+l)/2.,h/2.), (0, 1, 0), TIME_SERIES (data))

# create constraints solver
slv = NEWTON_SOLVER ()

# run simulation
RUN (solfec, slv, stop)

# post-process results
if not VIEWER() and solfec.mode == 'READ':
  data = []
  for b in bodies:
    data.append ((b, 'KINETIC'))
    data.append ((b, 'INTERNAL'))
  th = HISTORY (solfec, data, 0, stop)
  fq = []
  for t in th[0]: fq.append (lofq + (hifq-lofq)*(t/stop))
  for i in range (0,nbodl):
    plt.clf ()
    plt.plot (fq, th[2*i+1], lw = 2, label = 'kinetic')
    plt.plot (fq, th[2*i+2], lw = 2, label = 'internal')
    plt.legend(loc = 'best')
    plt.xlabel ('Frequency $(Hz)$')
    plt.ylabel ('Energy $(J)$')
    plt.savefig ('out/sweep/ene'+str(i)+'.svg')
