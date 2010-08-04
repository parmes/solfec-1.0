# example of a motion restrained to a plane

step = 1E-3

stop = 0.5

plane_normal = (0, 0, 1)

solfec = SOLFEC ('DYNAMIC', step, 'out/inplane')

solver = GAUSS_SEIDEL_SOLVER (1E-5, 100)

material = BULK_MATERIAL (solfec)

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5, restitution = 0.25)

shape = HULL ([0, 0, 0,
               1, 0, 0,
	       1, 1, 0,
	       0, 1, 0,
	       0, 0, 1,
	       1, 0, 1,
	       1, 1, 1,
	       0, 1, 1], 1, 1)

oshp1 = TRANSLATE (SCALE (COPY (shape), (10, 5, 0.5)), (-0.25, -1, -0.75))
oshp2 = TRANSLATE (SCALE (COPY (shape), (0.25, 5, 0.25)), (9.5, -1, 0.25))
oshp3 = TRANSLATE (SCALE (COPY (shape), (9.5, 0.25, 0.25)), (0., -1, 0.25))
oshp4 = TRANSLATE (SCALE (COPY (shape), (9.5, 0.25, 0.25)), (0., 3.75, 0.25))
oshp5 = TRANSLATE (SCALE (COPY (shape), (0.25, 5, 0.25)), (-0.25, -1, 0.25))
BODY (solfec, 'OBSTACLE', [oshp1, oshp2, oshp3, oshp4, oshp5], material)

bod1 = BODY (solfec, 'RIGID', COPY (shape), material)
FIX_DIRECTION (bod1, (0, 0, 0), plane_normal)
FIX_DIRECTION (bod1, (1, 0, 0), plane_normal)
FIX_DIRECTION (bod1, (1, 1, 0), plane_normal)
INITIAL_VELOCITY (bod1, (2.5, 0, 0), (0, 0, 1))

bod2 = BODY (solfec, 'RIGID', TRANSLATE (COPY (shape), (0, 2, 0)), material)
FIX_DIRECTION (bod2, (0, 2, 0), plane_normal)
FIX_DIRECTION (bod2, (1, 2, 0), plane_normal)
FIX_DIRECTION (bod2, (1, 3, 0), plane_normal)
INITIAL_VELOCITY (bod2, (0.5, 0, 0), (0, 0, -1))

PUT_RIGID_LINK (bod1, bod2, (0.5, 0.5, 1), (0.5, 2.5, 1))

RUN (solfec, solver, stop)
