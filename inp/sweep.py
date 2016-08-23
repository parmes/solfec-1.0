# ------------------ #
# sine sweep example #
# ------------------ #
import matplotlib.pyplot as plt
from math import sin, cos, pi

step = 1E-4  # time step
stop = 5.0   # duration of the simulation
damp = 1E-8  # amount of stiffness proportional damping
lofq = 5     # low frequency for the sweep
hifq = 15    # high frequency for the sweep
amag = 10.0  # acceleration magnitude
nbodl = 10   # number of bodies
nele = 10    # number of elements per body (along y)
nmod = 6     # number of linear modes to produce
l = 0.1      # length of one body
w = 0.1      # widhth of one body
h = 0.1      # height of one body
gap = 0.001  # gap
ostep = 1E-3 # output step
wavg = 0.01  # energy averaging time window [t-wavg/2, t+wavg/2]
fstop = 7.0  # end frequency for averaging velocities
             # (hand tunded in order to cover the pre-drop area)

GEOMETRIC_EPSILON (1E-9) # tiny geometrical tolerance (<< gap)

solfec = SOLFEC ('DYNAMIC', step, 'out/sweep')
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = 0.0)
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

t0 = 0.0
while t0 < stop:
  a0 = amag * sin (2*pi*(lofq+(hifq-lofq)*t0/stop)*t0)
  t0 += step
  a1 = amag * sin (2*pi*(lofq+(hifq-lofq)*t0/stop)*t0)
  if a1 < a0: break # find first acceleration maximum (integrated velocity will be more symmetrical about zero)

t = t0
v = 0.0
vt = []
va = []
vv = []
vf = []
while t < stop+t0:
  x = t + step/2. # mid-step time
  a = amag * sin (2*pi*(lofq+(hifq-lofq)*x/stop)*x) # mid-step acceleration
  v = v + a * step # mid-step integration of dv / dt = a into v
  vt.append (t-t0)
  va.append (a)
  vv.append (v)
  vf.append (lofq + (hifq-lofq)*(t/stop))
  t += step

for i in range (0, len(vv)):
  if vv[i+1] < vv[i]: break # find first velocity maximum

for j in range (0, len(vv)-i):
  vv[j] = vv[j+i] # shift velocity so it starts from the first maximum

while i > 0:
  vt.pop() # remove last i items from lists
  va.pop()
  vv.pop()
  vf.pop()
  i -= 1

# after velocity has been trimmed, produce displacement envelope
vd = []
d = 0.0
for v in vv:
 d = d + v * step  # integration of dd / dt = v
 vd.append (d)

# displacement has positive drift => find tangens of the positive drift angle 
i = len(vd)-1
while vd[i-1] > vd[i]: i -= 1 # first maximum
while vd[i-1] < vd[i]: i -= 1 # previous minimum
j = i
while vd[j-1] > vd[i]: j += 1 # previous maximum

# shift velocity down by the tangens of the drift angle
vshift = (vd[i]+vd[j]) / (vt[i]+vt[j])
for i in range (0, len(vv)): vv[i] -= vshift

# after velocity has been shifted down, produce displacement envelope
vd = []
d = 0.0
for v in vv:
 d = d + v * step  # integration of dd / dt = v
 vd.append (d)

if not VIEWER ():
  plt.clf ()
  plt.plot (vt, va)
  plt.xlim ((0, stop))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('acceleration $(m/s^2)$')
  plt.savefig ('out/sweep/acc.png')
  plt.clf ()
  plt.plot (vt, vv)
  plt.xlim ((0, stop))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('velocity $(m/s)$')
  plt.savefig ('out/sweep/vel.png')
  plt.clf ()
  plt.plot (vf, vd)
  plt.xlim ((lofq, hifq))
  plt.xlabel ('Frequency $(Hz)$')
  plt.ylabel ('displacement $(m)$')
  plt.savefig ('out/sweep/dsp.png')

# apply velocity constraint corresponding to the acceleration sin sweep
data = []
for (t, v) in zip (vt, vv): data += [t, v] # [t0,v0,t1,v1,...]
SET_VELOCITY (body, (w/2.,-(gap+l)/2.,h/2.), (0, 1, 0), TIME_SERIES (data))

# create constraints solver
slv = NEWTON_SOLVER ()

# output results every 'ostep'
OUTPUT (solfec, ostep)

# run simulation
RUN (solfec, slv, stop)

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
    plt.savefig ('out/sweep/ene'+str(k)+'.png')

    plt.clf ()
    plt.plot (fq, vy, lw = 2)
    plt.xlim ((lofq, hifq))
    plt.xlabel ('Frequency $(Hz)$')
    plt.ylabel ('Velocity vy $(m/s)$')
    plt.savefig ('out/sweep/vy'+str(k)+'.png')


    # output to file velocity for the 5th body
    if k == 4:
      fqout = open ('out/sweep/fq.txt', 'w')  
      vyout = open ('out/sweep/vy4.txt', 'w')  
      for item in fq: fqout.write (str(item)+'\n')
      for item in vy: vyout.write (str(item)+'\n')

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
