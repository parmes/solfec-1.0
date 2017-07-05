c = 0.2
step = 1E-4
stop = 0.03
nsid = 2
GEOMETRIC_EPSILON (1E-6)

#sv = GAUSS_SEIDEL_SOLVER (1E-4, 1000, 1E-8)
sv = NEWTON_SOLVER ()
solfec = SOLFEC ('DYNAMIC', step, 'out/crack')
bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E7, poisson = 0.3, density = 500)
surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3)

table = HEX ([0, 0,  0,
	      c, 0,  0,
	      c, c,  0,
	      0, c,  0,
              0, 0, -0.1*c,
	      c, 0, -0.1*c,
	      c, c, -0.1*c,
	      0, c, -0.1*c], 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])

SCALE (table, (nsid, nsid, 1.0))
BODY (solfec, 'OBSTACLE', table, bulkmat)

for i in range (0, nsid):
  for j in range (0, nsid):
    shp = PIPE ((0, 0, 0), (0, 0, c), 0.2*c, 0.3*c, 4, 8, 2, 2, [2, 2, 2, 2])
    TRANSLATE (shp, (i*c+0.5*c, j*c+0.5*c, 0.5*c))
    bod = BODY (solfec, 'FINITE_ELEMENT', shp, bulkmat)
    SIMPLIFIED_CRACK (bod, bod.center, (1, 0, 0), (3, 3), 'TENSILE', ft=1E5)
    SIMPLIFIED_CRACK (bod, bod.center, (0, 1, 0), (3, 3), 'TENSILE', ft=1E5)
    SIMPLIFIED_CRACK (bod, bod.center, (0, 0, 1), (3, 3), 'TENSILE', ft=1E5)

GRAVITY (solfec, (0, 0, -1000))
OUTPUT (solfec, 2 * step)
RUN (solfec, sv, stop)

if not VIEWER() and solfec.mode == 'READ':
  try:
    import matplotlib.pyplot as plt

    dur = DURATION (solfec)
    th = HISTORY (solfec, [(solfec, 'KINETIC'), (solfec, 'INTERNAL'), (solfec, 'EXTERNAL'), 'BODS'], dur [0], dur [1], progress = 'ON')
    plt.plot (th [0], th [1], label='KIN')
    plt.plot (th [0], th [2], label='INT')
    plt.plot (th [0], th [3], label='EXT')
    tot = []
    for i in range(0, len (th[0])): tot.append (th[1][i] + th[2][i] - th[3][i])
    plt.plot (th [0], tot, label='TOT')
    plt.plot (th [0], th [4], label='NBD')
    plt.axis (xmin = dur [0], xmax = dur [1])
    plt.xlabel ('Time [s]')
    plt.ylabel ('Energy [J]')
    plt.legend(loc = 'upper right')
    plt.savefig ('out/crack/energy.eps')

  except ImportError:
    pass # no reaction
