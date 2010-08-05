# cylinder section example

def CYLINDER_SECTION (center, radius, thick, height, angle, div_thick, div_height, div_angle, volid = 0, surfid = 0):

  nodes = [0, 0, 0,
	   1, 0, 0,
	   1, 1, 0,
	   0, 1, 0,
	   0, 0, 1,
	   1, 0, 1,
	   1, 1, 1,
	   0, 1, 1]

  msh = HEX (nodes, div_thick, div_angle, div_height, volid, [surfid, surfid, surfid, surfid, surfid, surfid])
  SCALE (msh, (thick, 3.14159 * (angle / 180.0) * radius, height))
  TRANSLATE (msh, (-thick-radius, 0, 0))
  BEND (msh, (0, 0, 0), (0, 0, -1), angle)

  return TRANSLATE (msh, center)


# main module

step = 0.001
stop = 1

sol = SOLFEC ('DYNAMIC', step, 'out/cycsec')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 15E9,
		      poisson = 0.3,
		      density = 2E3)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

CONTACT_SPARSIFY (sol, 0.005, 0.005)

shp = CYLINDER_SECTION ((0, 0, 0), 1, 1, 2, 90, 2, 4, 10)

bod = BODY (sol, 'RIGID', shp, bulk)

shp = CYLINDER_SECTION ((0, 0, -3), 1, 1, 0.5, 90, 2, 1, 10)

bod = BODY (sol, 'OBSTACLE', shp, bulk)

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000)

GRAVITY (sol, (0, 0, -10))

OUTPUT (sol, step)

RUN (sol, gs, stop)
