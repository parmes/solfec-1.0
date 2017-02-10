from math import sqrt, sin, cos

def expmap(O):
  W = (0, O[2], -O[1], -O[2], 0, O[0], O[1], -O[0], 0)
  WW = (W[0]*W[0]+W[3]*W[1]+W[6]*W[2],
        W[1]*W[0]+W[4]*W[1]+W[7]*W[2],
        W[2]*W[0]+W[5]*W[1]+W[8]*W[2],
        W[0]*W[3]+W[3]*W[4]+W[6]*W[5],
        W[1]*W[3]+W[4]*W[4]+W[7]*W[5],
        W[2]*W[3]+W[5]*W[4]+W[8]*W[5],
        W[0]*W[6]+W[3]*W[7]+W[6]*W[8],
        W[1]*W[6]+W[4]*W[7]+W[7]*W[8],
        W[2]*W[6]+W[5]*W[7]+W[8]*W[8])
  eye = (1, 0, 0, 0, 1, 0, 0, 0, 1)
  dp = O[0]*O[0]+O[1]*O[1]+O[2]*O[2]
  ln = sqrt(dp);
  if ln < 1E-15:
    R = eye
  else:
    R = tuple(eye[i] + (sin(ln)/ln)*W[i] + ((1-cos(ln))/dp)*WW[i] for i in range (0,9))
  return R

def torque(q, u, time, step):
  a = (q[6], q[7], q[8])
  b = (0, 0, 1)
  c = (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
  return (0, 0, 0, -20*c[0], -20*c[1], -20*c[2], 0, 0, 0)

scheme = ['RIG_POS', 'RIG_NEG', 'RIG_IMP']
step = [0.001, 0.003, 0.009]
stop = 10

T = list()
H = ((list(),list(),list()),
     (list(),list(),list()),
     (list(),list
     (),list()))

def callback(sol, bod, i, j):
  if i == 0 and j == 0: T.append(sol.time)
  ene = ENERGY(sol)
  L33 = bod.conf[8]
  H[i][j].append(ene[0]+20*L33)
  return 1

for i in range (0,len(scheme)):
  for j in range (0,len(step)):

    sol = SOLFEC ('DYNAMIC', step[j], 'out/heavy-symmetrical-top-%s-%.3f' % (scheme[i].lower(), step[j]))
    mat = BULK_MATERIAL (sol)

    sph = SPHERE ((0, 0, 0), 1, 1, 1) # a mock shape
    bod = BODY (sol, 'RIGID', sph, mat)
    bod.scheme = scheme[i]
    BODY_CHARS (bod, 1, 1, (0, 0, 0), (5, 0, 0, 0, 5, 0, 0, 0, 1))
    bod.conf = expmap((0.05, 0, 0)) + (0, 0, 0)
    bod.velo = (0, 0, 50, 0, 0, 0)
    FORCE (bod, 'SPATIAL', (0, 0, 0), (0, 0, 0), torque)

    CALLBACK (sol, 0.05, (sol, bod, i, j), callback)
    ns = NEWTON_SOLVER ()
    OUTPUT (sol, 0.05)
    RUN (sol, ns, stop)

L33 = expmap((0.05, 0, 0))[8]
H0 = 0.5*50**2 + 20*L33
print 'Initial R(3,3) is', L33
print 'Initial Hamiltonian value is', H0
print 'Plotting Hamiltonian time histories ...'

try:
  import matplotlib.pyplot as plt
  for i in range(0,len(scheme)):
    plt.clf()
    #ax = plt.gca()
    #ax.ticklabel_format(useOffset=False)
    for j in range(0,len(step)):
      plt.plot (T, H[i][j], label='h = %.3f'%step[j])
    plt.axis (xmin = 0, xmax = stop)
    plt.title (scheme[i])
    plt.xlabel ('Time [s]')
    plt.ylabel ('Hamiltonian [J]')
    plt.legend(loc = 'lower left')
    path = 'validation/heavy-symmetrical-top/' + scheme[i].lower() + '.png'
    print 'Saving', path
    plt.savefig (path)
except (ImportError, RuntimeError):
  import sys
  print "Unexpected error:", sys.exc_info()[1]
  print "Plotting has failed!"
  pass
