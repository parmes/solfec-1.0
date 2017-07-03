# splitting tests
import sys

sv = GAUSS_SEIDEL_SOLVER (1E-4, 1000)

sol = SOLFEC ('DYNAMIC', 0.001, 'out/splits')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E6,
		      poisson = 0.3,
		      density = 1E3)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

GRAVITY (sol, (0, 0, -10))

def numbered_pipe (base):
  shp = PIPE (base, (0, 0, 2), 0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])
  bod = BODY (sol, 'FINITE_ELEMENT', shp, bulk)
  msh = bod.mesh
  for i in range (0, msh.nnod):
    DISPLAY_POINT (bod, msh.node (i), '%d'%i)

def split_pipe (base, faces):
  shp = PIPE  (base, (0, 0, 2), 0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])
  out = MESH_SPLIT (shp, faceset = faces)
  for msh in out:
    BODY (sol, 'FINITE_ELEMENT', msh, bulk)

  shp = PIPE  (TRANSLATE (base, (0, 0, -1)), (0, 0, -1), 0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])
  BODY (sol, 'OBSTACLE', shp, bulk)

# numbered_pipe ((0, 0, 0))
f = [[0, 1, 17, 16], [2, 3, 19, 18], [4, 5, 21, 22],
     [6, 7, 23, 22], [8, 9, 25, 24], [10, 11, 27, 26],
     [12, 13, 29, 28], [14, 15, 31, 30]]

split_pipe ((0, 0, 0), [f[0]])
#split_pipe ((4, 0, 0), [f[1]])
#split_pipe ((8, 0, 0), [f[2]])
#split_pipe ((0, 4, 0), [f[3]])
#split_pipe ((4, 4, 0), [f[4]])
#split_pipe ((8, 4, 0), [f[5]])

RUN (sol, sv, 1.0)
