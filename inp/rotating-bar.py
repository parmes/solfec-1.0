# rotating bar example
# testing convergence rate of BC and RO formulations
# --------------------------------------------------

a = 0.1
nmesh = 2
lmesh = 8
nmod = 12

nodes = [-1.1*a, -a, -lmesh*a,
          1.3*a, -0.9*a, -lmesh*a,
          a,  0.7*a, -lmesh*a,
         -0.8*a,  1.2*a, -lmesh*a,
         -1.2*a, -0.8*a,  lmesh*a,
          0.7*a, -a,  lmesh*a,
          0.9*a,  1.3*a,  lmesh*a,
         -a,  1.1*a,  lmesh*a]

solfec_stack = []

def runtest (formulation, step, duration):
  mesh = HEX (nodes, nmesh, nmesh, lmesh*nmesh, 0, [0, 0, 0, 0, 0, 0])
  solfec = SOLFEC ('DYNAMIC', step, 'out/rotating-bar/' + formulation + '/' + str (step))
  solfec_stack.append (solfec)
  bulk = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 200E9, poisson = 0.26, density = 7.8E3) # A36 steel
  if formulation == 'RO':
    bod = BODY (solfec, 'FINITE_ELEMENT', COPY (mesh), bulk)
    data = MODAL_ANALYSIS (bod, nmod, 'out/rotating-bar/modal_data')
    DELETE (solfec, bod)
    bod = BODY (solfec, 'FINITE_ELEMENT', mesh, bulk, form = 'RO', modal = data)
  else:
    bod = BODY (solfec, 'FINITE_ELEMENT', mesh, bulk, form = formulation)
    bod.scheme = 'DEF_LIM'
  INITIAL_VELOCITY (bod, (0, 0, 0), (1, 0, 0))
  bod.damping = step # O(h) damping stabilizes the convergence rate (removes oscilations on the step level => investigate)
  sv = NEWTON_SOLVER ()
  RUN (solfec, sv, duration)
  if solfec.mode == 'READ':
    SEEK (solfec, duration)
  return bod.conf

def norm(q):
  sum = 0.0
  for x in q: sum += x**2
  return sum**0.5

def diff(a, b):
  c = []
  for x, y in zip (a, b): c.append (x-y)
  return c

h0 = 1.0 / (2.0 ** 16)
d0 = 1.0 / (2.0 ** 4)
q0 = runtest ('BC', h0, d0)

dq = []
for i in range (7, 14):
  h =  1.0 / (2.0 ** i)
  q = runtest ('BC', h, d0)
  dq.append (norm(diff(q0, q)))

print dq

if not VIEWER():
  for i in range (0, len(dq)-1):
    print dq[i] / dq[i+1]
