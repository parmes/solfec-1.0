# sphere sliding inside of a pipe
from math import sqrt

step = 1E-3
duration = 10

solfec = SOLFEC ('DYNAMIC', step, 'out/sphere-in-pipe')
solver = GAUSS_SEIDEL_SOLVER (1E-6, 1000)
mat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0, restitution = 0)

pip = PIPE ((0, -0.5, 0), (0, 1, 0), 10, 1, 1, 512, 1, 0, [0, 0, 0, 0])
BODY (solfec, 'OBSTACLE', pip, mat)

sph = SPHERE ((-9.5, 0, 0), 0.5, 1, 1)
bod = BODY (solfec, 'RIGID', sph, mat)
INITIAL_VELOCITY (bod, (0, 0, -50), (0, 0, 0))

RUN (solfec, solver, duration)

if not VIEWER() and solfec.mode == 'READ':
  try:
    import matplotlib.pyplot as plt
     
    th = HISTORY (solfec, [(bod, bod.center, 'CX'), (bod, bod.center, 'CZ'), (solfec, 'KINETIC')], 0, duration)
    i = 0
    rh = []
    for t in th [0]:
      x = th [1][i]
      z = th [2][i]
      r = sqrt (x*x + z*z)
      rh.append (r)
      i = i + 1
 
    plt.plot (th [0], rh)
    plt.title ('Distance of sphere center from circle center')
    plt.xlabel ('Time')
    plt.ylabel ('Distance')
    plt.savefig ('out/sphere-in-pipe/sphere_distance_history.eps')

    plt.clf ()
    plt.plot (th [0], th [3])
    plt.title ('Kinetic energy')
    plt.xlabel ('Time')
    plt.ylabel ('Energy')
    plt.savefig ('out/sphere-in-pipe/sphere_energy.eps')
  except ImportError:
    pass
