# extended pseudo-rigid test - one dimensional bar

N = 5

def one_d_bar_create (n, material, solfec):

  a = 0.1
  b = 0.1
  c = 0.1
  sur = 1
  vol = 1
  shp = []

  for i in range (n):

    cvx = CONVEX (
	    [0, 0, 0,
	     a, 0, 0,
	     a, b, 0,
	     0, b, 0,
	     0, 0, c,
	     a, 0, c,
	     a, b, c,
	     0, b, c],
	    [4, 0, 3, 2, 1, sur,
	     4, 1, 2, 6, 5, sur,
	     4, 2, 3, 7, 6, sur,
	     4, 3, 0, 4, 7, sur,
	     4, 0, 1, 5, 4, sur,
	     4, 4, 5, 6, 7, sur], vol)

    TRANSLATE (cvx, (0, 0, i*c))

    shp.append (cvx);

  bod = BODY (solfec, 'EXTENDED_PSEUDO_RIGID', shp, material)

  FORCE (bod, 'SPATIAL', (.5*a, .5*a, 0), (0, 0, 1), 100)
  FORCE (bod, 'SPATIAL', (.5*a, .5*a, N*c), (0, 0, -1), 10)

  FIX_POINT (solfec, bod, (.5*a, .5*b, 0))
  #FIX_DIRECTION (solfec, bod, (a, 0, 0), (0, 0, 1))
  #FIX_DIRECTION (solfec, bod, (0, b, 0), (0, 0, 1))
  #FIX_DIRECTION (solfec, bod, (0, b, 0), (1, 0, 0))

### main module ###

step = 1E-3

solfec = SOLFEC ('DYNAMIC', step, 'out/epr0')

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E3, poisson = 0.25, density = 1E3)

GRAVITY (solfec, (0, 0, -1), 9.81)

#import rpdb2; rpdb2.start_embedded_debugger('a')

one_d_bar_create (N, bulkmat, solfec)

def gscallback (gs):
  print gs.error
  return 0

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'EXIT', diagsolver = 'PROJECTED_GRADIENT')

RUN (solfec, gs, 1E3 * step)
