# tetrahedral mesh example
step = 1E-3
stop = 0.5
splitme = 0

solfec = SOLFEC ('DYNAMIC', step, 'out/stlpart')

solver = GAUSS_SEIDEL_SOLVER (1E-5, 100)

material = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 200E9, poisson = 0.25, density = 20E3)

SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5, restitution = 0.25)

if splitme:
  mesh = TETRAHEDRALIZE ('inp/mesh/hinge.stl', 'out/stlpart/mesh0.dat')
  (a, b) = SPLIT (mesh, (5, 5, 5), (1, 10, 10), (1, 2))
  a = TETRAHEDRALIZE (a, 'out/stlpart/mesh1.dat', quality = 1.5)
  b = TETRAHEDRALIZE (b, 'out/stlpart/mesh2.dat', quality = 1.5)
  TRANSLATE (b, (0, 0, 2))
  ROTATE (b, (0, 0, 0), (0, 0, 1), 45)
  BODY (solfec, 'FINITE_ELEMENT', a, material)
  BODY (solfec, 'FINITE_ELEMENT', b, material)
else:
  mesh = TETRAHEDRALIZE ('inp/mesh/hinge.stl', 'out/stlpart/mesh0.dat', quality = 1.5)
  BODY (solfec, 'FINITE_ELEMENT', mesh, material)

RUN (solfec, solver, stop)
