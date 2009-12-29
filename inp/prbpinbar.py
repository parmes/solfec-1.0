# pinned bar example

def pinned_bar_create (material, solfec):

  nodes = [-0.05, -0.05, 0.0,
            0.05, -0.05, 0.0,
            0.05,  0.05, 0.0,
	   -0.05,  0.05, 0.0,
	   -0.05, -0.05, 1.0,
	    0.05, -0.05, 1.0,
	    0.05,  0.05, 1.0,
	   -0.05,  0.05, 1.0]

  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, 0]

  point = (0, 0, 0.75)

  vector = (0, 1, 0)

  fix1 = (0, -0.05, 0.75)

  fix2 = (0, 0.05, 0.75)

  msh = MESH (nodes, elements, 0)

  ROTATE (msh, point, vector, -30)

  bod = BODY (solfec, 'PSEUDO_RIGID', msh, material)
  bod.scheme = 'DEF_IMP'

  FIX_POINT (solfec, bod, fix1)
  FIX_POINT (solfec, bod, fix2)

  msh = MESH (nodes, elements, 0)

  vector = (-0.15, 0, -0.2)

  TRANSLATE (msh, vector)

  BODY (solfec, 'OBSTACLE', msh, material)

  return bod

### main module ###

step = 0.001
stop = 10.0

solfec = SOLFEC ('DYNAMIC', step, 'out/pinnedbar')

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.3)

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)

GRAVITY (solfec, (0, 0, -1), 9.8)

#import rpdb2; rpdb2.start_embedded_debugger('a')

bod = pinned_bar_create (bulkmat, solfec)

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'EXIT', diagsolver = 'PROJECTED_GRADIENT')

OUTPUT (solfec, 0.001)

RUN (solfec, gs, stop)

#if not VIEWER() and solfec.mode == 'READ':
  #kin = ENERGY_HISTORY (solfec, 'KINETIC', 0, stop)
  #pot = ENERGY_HISTORY (solfec, 'POTENTIAL', 0, stop)
  #ext = ENERGY_HISTORY (solfec, 'EXTWORK', 0, stop)
  #fri = ENERGY_HISTORY (solfec, 'FRICWORK', 0, stop)
