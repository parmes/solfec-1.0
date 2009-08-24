# a drum mixing example

def drum_create (material, solfec):
  box = HULL ([0, 0, 0,
               1, 0, 0,
	       1, 1, 0,
	       0, 1, 0,
               0, 0, 1,
               1, 0, 1,
	       1, 1, 1,
	       0, 1, 1], 0, 0)

  wallx1 = SCALE (COPY (box), (1, 0.1, 1))
  wallx2 = TRANSLATE (COPY (wallx1), (0, 0.9, 0))
  wally1 = SCALE (COPY (box), (0.1, 0.8, 1))
  wally2 = TRANSLATE (COPY (wally1), (0.9, 0, 0))
  TRANSLATE ([wally1, wally2], (0, 0.1, 0))
  wallz1 = SCALE (COPY (box), (1, 1, 0.1))
  wallz2 = TRANSLATE (COPY (wallz1), (0, 0, 1))
  TRANSLATE (wallz1, (0, 0, -0.1))

  shape = [wallx1, wallx2, wally1, wally2, wallz1, wallz2]
  bod = BODY (solfec, 'RIGID', shape, material)
  FIX_POINT (solfec, bod, (0.0, 0.5, 0.5))
  FIX_POINT (solfec, bod, (1.0, 0.5, 0.5))
  tms = TIME_SERIES ([0, 10, 100, 10])
  TORQUE (bod, 'SPATIAL', (1, 0, 0), tms)

def spheres_create (material, solfec):

  rng = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

  for x in rng:
    for y in rng:
      for z in rng:
	sph = SPHERE ((x, y, z), 0.025, 1, 1);
	BODY (solfec, 'RIGID', sph, material)

### main module ###

step = 0.001
skip = 10
dura = 2.0

solfec = SOLFEC ('DYNAMIC', step, 'out/drum')

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.4)

drumat = BULK_MATERIAL (solfec)

sphmat = BULK_MATERIAL (solfec, density = drumat.density / 10)

GRAVITY (solfec, (0, 0, -1), 9.8)

#import rpdb2; rpdb2.start_embedded_debugger('a')

drum_create (drumat, solfec)
spheres_create (sphmat, solfec)

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'EXIT')

def callback (sol):
  return 1

if not VIEWER(): CALLBACK (solfec, step * skip, solfec, callback)

OUTPUT (solfec, 0.02)

UNPHYSICAL_PENETRATION (solfec, 0.01)

RUN (solfec, gs, dura)
