# SOLFEC example --> articulated pendulum
# equivalent of parmec/examples/articulated_pendulum.py

nele=2 # number of elements
step= 0.0005
stop=10.

solfec = SOLFEC ('DYNAMIC', step, 'out/articulated_pendulum')

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF',
           young = 1E9, poisson = 0.25, density = 1E3)

bodies = []

for i in range (0, nele):
  nodes = [i,   0,   0,
	   i+1, 0,   0,
	   i+1, 0.1, 0,
	   i,   0.1, 0,
	   i,   0,   0.1,
	   i+1, 0,   0.1,
	   i+1, 0.1, 0.1,
	   i,   0.1, 0.1]

  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, 0]

  colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

  msh = MESH (nodes, elements, colors)

  bod = BODY (solfec, 'RIGID', msh, bulkmat)

  if i: 
    PUT_RIGID_LINK (bod, bodies[i-1], (i, 0.05, 0.05), (i, 0.05, 0.05))
    CONTACT_EXCLUDE_BODIES (bod, bodies[i-1])

  bodies.append(bod)

FIX_POINT (bod, (i+1, 0.05, 0.05))

GRAVITY (solfec, (0., 0., -10.))

ns = NEWTON_SOLVER()

RUN (solfec, ns, stop)

if not VIEWER() and solfec.mode == 'READ':

  print 'Generating plots ...',
  import sys
  sys.stdout.flush()

  th = HISTORY (solfec, [(bodies[0], (1, 0.05, 0.05), 'CX'),
                         (bodies[0], (1, 0.05, 0.05), 'CY'),
			 (bodies[0], (1, 0.05, 0.05), 'CZ'),
			 (bodies[1], (1, 0.05, 0.05), 'CX'),
			 (bodies[1], (1, 0.05, 0.05), 'CY'),
			 (bodies[1], (1, 0.05, 0.05), 'CZ')], 0, stop)
  try:
    import matplotlib.pyplot as plt

    dp = []
    for (x0,y0,z0,x1,y1,z1) in zip(th[1],th[2],th[3],th[4],th[5],th[6]):
      dp.append (((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)**0.5)
    plt.clf ()
    plt.plot (th[0], dp)
    plt.xlim ((0, th[0][-1]))
    plt.xlabel ('time (s)')
    plt.ylabel ('|p(body 0) - p(body 1)| (m)')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.title ('Joint point p motion difference (bodies: 0,1)')
    plt.savefig ('out/articulated_pendulum/articulated_pendulum_dp.png')

    print 'ok.'

  except:
    print 'failed. (check matplotlib)'
