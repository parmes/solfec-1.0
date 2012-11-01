# Advanced Manufacturing Techniques contact example
from math import pi
step = 0.001
stop = 6.5

sol = SOLFEC ('DYNAMIC', step, 'out/amt')

bulk = BULK_MATERIAL (sol, model = 'KIRCHHOFF', young = 1E9, poisson = 0.2, density = 100)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.1)

nodes = [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1]
msh = HEX (nodes, 30, 1, 2, 0, [0, 0, 0, 0, 0, 0])
SCALE (msh, (pi*4.0, 1, 1))
BEND (msh, (0, 0, 5.0), (0, -1, 0), 180)
bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk)
bod.selfcontact = 'ON'
bod.scheme = 'DEF_LIM'
bod.damping = 1E-3

msh = HEX (nodes, 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])
SCALE (msh, (pi*2.0, 1, 1))
TRANSLATE (msh, (0, 0, -1))
bod = BODY (sol, 'OBSTACLE', msh, bulk)

msh = HEX (nodes, 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])
SCALE (msh, (pi*2.0, 1, 1))
TRANSLATE (msh, (0, 0, 10))
p1 = msh.node (0)
p2 = msh.node (1)
p3 = msh.node (3)
p4 = msh.node (4)
bod = BODY (sol, 'RIGID', msh, bulk)
FIX_DIRECTION (bod, p1, (0, -1, 0))
FIX_DIRECTION (bod, p2, (0, -1, 0))
FIX_DIRECTION (bod, p4, (0, -1, 0))
FIX_DIRECTION (bod, p2, (1, 0, 0))
FIX_DIRECTION (bod, p3, (1, 0, 0))
FIX_DIRECTION (bod, p4, (-1, 0, 0))
con = SET_VELOCITY (bod, bod.center, (0, 0, 1), -1)
z0 = con.point [2]


sv = GAUSS_SEIDEL_SOLVER (1E-4, 200)

OUTPUT (sol, step)
RUN (sol, sv, stop)

if not VIEWER() and sol.mode == 'READ':

  t = 0.0
  SEEK (sol, t)
  d = []
  r = []
  while t <= stop:
    z1 = con.point [2]
    Rz = con.R [2]
    d.append (z0 - z1)
    r.append (-Rz)
    FORWARD (sol, 1)
    t += step

  try:
    import matplotlib.pyplot as plt
    plt.plot (d, r)
    plt.xlabel ('Displacement [m]')
    plt.ylabel ('Force [N]')
    plt.savefig ('out/amt/dr.eps')
  except ImportError:
    pass # no reaction
