# simple shock absorber example

step = 1E-3
stop = 1.0

solfec = SOLFEC ('DYNAMIC', step, 'out/dixdir')

bulk = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E9, poisson = 0.3, density = 1E3)

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.0, restitution = 0.0)

box = HEX ([0, 0, 0,
            1, 0, 0,
            1, 1, 0,
            0, 1, 0,
            0, 0, 1,
            1, 0, 1,
            1, 1, 1,
            0, 1, 1], 1, 1, 1, 1, [1]*6)

s = SCALE (COPY(box), (0.1, 0.1, 1))
s = TRANSLATE (s, (-0.05, -0.05, 0.0))
b1 = BODY (solfec, 'RIGID', s, bulk)

s = SCALE (COPY(box), (0.05, 0.05, 1))
s = TRANSLATE (s, (-0.025, -0.025, 0.5))
b2 = BODY (solfec, 'RIGID', s, bulk)

CONTACT_EXCLUDE_BODIES (b1, b2)

FIX_DIRECTION (b1, (0, 0, 0.5), (1, 0, 0), b2, (0, 0, 0.5))
FIX_DIRECTION (b1, (0, 0, 0.5), (0, 1, 0), b2, (0, 0, 0.5))
FIX_DIRECTION (b1, (0, 0, 1), (1, 0, 0), b2, (0, 0, 1))
FIX_DIRECTION (b1, (0, 0, 1), (0, 1, 0), b2, (0, 0, 1))

def springfunc (stroke, velocity):
  return -1000*stroke - 100*velocity

PUT_SPRING (b1, (0, 0, 1), b2, (0, 0, 0.5), springfunc, (-0.35, 0.35))

s = SCALE (COPY(box), (1, 1, 0.1))
s = TRANSLATE (s, (-0.5, -0.5, -1.1))
b3 = BODY (solfec, 'OBSTACLE', s, bulk)

GRAVITY (solfec, (0, 0, -10))

#slv = GAUSS_SEIDEL_SOLVER (1E-3, 200)
slv = NEWTON_SOLVER ()

RUN (solfec, slv, stop)
