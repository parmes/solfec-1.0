c = 0.1
d = 0.2 * c
step = 1E-5
stop = 1.0
velo = 1E-6
GEOMETRIC_EPSILON (1E-7)

sv = GAUSS_SEIDEL_SOLVER (1E-3, 1000, 1E-7)
#sv = NEWTON_SOLVER (maxiter = 500, theta = 0.1, locdyn = 'ON')

solfec = SOLFEC ('DYNAMIC', step, 'out/crack')

table = HEX ([-d, -d,  0,
	       c+d, -d,  0,
	       c+d,  c+d,  0,
	      -d,  c+d,  0,
              -d, -d, -d,
	       c+d, -d, -d,
	       c+d,  c+d, -d,
	      -d,  c+d, -d], 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E7, poisson = 0.3, density = 500)
surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3)

TRANSLATE (table, (0, 0, -0.5*c))
BODY (solfec, 'OBSTACLE', table, bulkmat)

box = HEX ([0, 0, 0,
	    c, 0, 0,
	    c, c, 0,
	    0, c, 0,
            0, 0, c,
	    c, 0, c,
	    c, c, c,
	    0, c, c], 3, 3, 3, 2, [2, 2, 2, 2, 2, 2])

box = TETRAHEDRALIZE (box, 'out/crack/tet1.dat', c**6, quality = 1.5)

bod = BODY (solfec, 'FINITE_ELEMENT', box, bulkmat)
#bod.nodecontact = 'ON'
p0 = TRANSLATE (bod.center, (c/2.0, 0, 0))
p1 = TRANSLATE (bod.center, (-c/2.0, 0, 0))
#SET_VELOCITY (bod, p0, (1, 0, 0), velo)
#SET_VELOCITY (bod, p1, (-1, 0, 0), velo)
#SIMPLIFIED_CRACK (bod, bod.center, (1, 0, 0), 1, 'TENSILE', 1E3, 10)
#SIMPLIFIED_CRACK (bod, bod.center, (0, 1, 0), 1, 'TENSILE', 1E3, 10)
#SIMPLIFIED_CRACK (bod, bod.center, (0, 0, 1), 1, 'TENSILE', 1E3, 10)

#FIXME: strangely preserved handing contact points afeter complete fragmentation
#FIXME: in the element to element contact mode
#FIXME: similarly in the node based contact mode

GRAVITY (solfec, (0, 0, -1000))
OUTPUT (solfec, 2 * step)
RUN (solfec, sv, stop)

if solfec.mode == 'READ':
  FORWARD (solfec, 1)
