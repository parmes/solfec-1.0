# wire twisting example

step = 0.0001
stop = 0.2

nwir = 3

mx = 1
my = 1
mz = 50

nodes = [0, 0, 0,
         1, 0, 0,
         1, 1, 0,
         0, 1, 0,
         0, 0, 1,
         1, 0, 1,
         1, 1, 1,
         0, 1, 1]

solfec = SOLFEC ('DYNAMIC', step, 'out/wire')

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E6, poisson = 0.25, density = 1E3)

solver = GAUSS_SEIDEL_SOLVER (1E-3, 100)

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.2)

tms = TIME_SERIES ([0, 0.1, 10*step, 0.1, 11*step, 0.0, stop, 0.0])

mesh0 = HEX (nodes, 1, 1, 1, 0, [0, 0, 0, 0, 0, 0])
SCALE (mesh0, (nwir*0.01, nwir*0.01, 0.01))
b0 = BODY (solfec, 'RIGID', mesh0, bulkmat)
TORQUE (b0, 'SPATIAL', (0, 0, -1), tms)

mesh1 = HEX (nodes, 1, 1, 1, 0, [0, 0, 0, 0, 0, 0])
SCALE (mesh1, (nwir*0.01, nwir*0.01, 0.01))
TRANSLATE (mesh1, (0, 0, 0.11))
b1 = BODY (solfec, 'RIGID', mesh1, bulkmat)
TORQUE (b1, 'SPATIAL', (0, 0, 1), tms)

for i in range (1,nwir):
  for j in range (1,nwir):
    mesh2 = HEX (nodes, mx, my, mz, 1, [1, 1, 1, 1, 1, 1])
    SCALE (mesh2, (0.002, 0.002, 0.1))
    TRANSLATE (mesh2, (i*0.01, j*0.01, 0.01))

    nset0 = []
    nset1 = []
    for k in range (0, (mx+1)*(my+1)):
      nset0.append (mesh2.node (k))
      nset1.append (mesh2.node ((mx+1)*(my+1)*(mz+1)-1-k))

    b2 = BODY (solfec, 'FINITE_ELEMENT', mesh2, bulkmat)

    CONTACT_EXCLUDE_BODIES (b0, b2)
    CONTACT_EXCLUDE_BODIES (b1, b2)

    for nod in nset0:
      PUT_RIGID_LINK (b2, b0, nod, nod)
    for nod in nset1:
      PUT_RIGID_LINK (b2, b1, nod, nod)

RUN (solfec, solver, stop)
