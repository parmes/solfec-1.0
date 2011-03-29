c = 0.2
d = 0.3
step = 1E-5
stop = 1.0
velo = 1E-6
GEOMETRIC_EPSILON (1E-6)

sv = GAUSS_SEIDEL_SOLVER (1E-3, 1000, 1E-7)
#sv = NEWTON_SOLVER (maxiter = 500, theta = 0.1, locdyn = 'ON')
solfec = SOLFEC ('DYNAMIC', step, 'out/crack')
bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E7, poisson = 0.3, density = 500)
surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3)

table = HEX ([-d, -d,  0,
	       d, -d,  0,
	       d,  d,  0,
	      -d,  d,  0,
              -d, -d, -0.1*d,
	       d, -d, -0.1*d,
	       d,  d, -0.1*d,
	      -d,  d, -0.1*d], 1, 1, 1, 1, [1, 1, 1, 1, 1, 1])

BODY (solfec, 'OBSTACLE', table, bulkmat)

shp = HEX ([-c, -c, 0,
	     c, -c, 0,
	     c,  c, 0,
	    -c,  c, 0,
            -c, -c, 2*c,
	     c, -c, 2*c,
	     c,  c, 2*c,
	    -c,  c, 2*c], 3, 3, 3, 2, [2, 2, 2, 2, 2, 2])
shp= PIPE ((0, 0, 0), (0, 0, c), 0.2*c, 0.3*c, 3, 8, 2, 2, [2, 2, 2, 2])
#shp = TETRAHEDRALIZE (shp, 'out/crack/tet1.dat', c**6, quality = 1.5)
TRANSLATE (shp, (0, 0, 0.5*c))
#(a, b) = SPLIT (shp, (0, 0, 0.5*c), (0, 0, 1))
#shp = TETRAHEDRALIZE (a, 'out/crack/tet1.dat', c**6, quality = 1.5)
#shp = a

bod = BODY (solfec, 'FINITE_ELEMENT', shp, bulkmat)
#bod.nodecontact = 'ON' #FIXME: wrong normal directions after fragmentation
#p0 = TRANSLATE (bod.center, (c/2.0, 0, 0))
#p1 = TRANSLATE (bod.center, (-c/2.0, 0, 0))
#SET_VELOCITY (bod, p0, (1, 0, 0), velo)
#SET_VELOCITY (bod, p1, (-1, 0, 0), velo)
SIMPLIFIED_CRACK (bod, bod.center, (1, 0, 0), 1, 'TENSILE', 1E5, 10)
SIMPLIFIED_CRACK (bod, bod.center, (0, 1, 0), 1, 'TENSILE', 1E5, 10)
SIMPLIFIED_CRACK (bod, bod.center, (0, 0, 1), 1, 'TENSILE', 1E5, 10)

GRAVITY (solfec, (0, 0, -1000))
OUTPUT (solfec, 2 * step)
RUN (solfec, sv, stop)

if solfec.mode == 'READ':
  FORWARD (solfec, 1)
