# MESH_SPLIT test
import sys

sv = GAUSS_SEIDEL_SOLVER (1E-4, 1000)

sol = SOLFEC ('DYNAMIC', 0.001, 'out/mesh_split')

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
  if out == None:
    print 'Warning --> spliting failed for face list:', faces
    sys.exit(1)
  else:
    for msh in out:
      BODY (sol, 'FINITE_ELEMENT', msh, bulk)

  shp = PIPE  (TRANSLATE (base, (0, 0, -1)), (0, 0, -1), 0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])
  BODY (sol, 'OBSTACLE', shp, bulk)

# numbered_pipe ((0, 0, 0))
f = [[0, 1, 17, 16], [2, 3, 19, 18], [4, 5, 21, 20],
     [6, 7, 23, 22], [8, 9, 25, 24], [10, 11, 27, 26],
     [12, 13, 29, 28], [14, 15, 31, 30]]

for i in range (0, 8):
  split_pipe ((0+4*(i%4), 0+4*(i>3), 0), [f[i]])

for i in range (0, 8):
  split_pipe ((0+4*(i%4), 8+4*(i>3), 0), [f[i], f[(i+1)%8]])

RUN (sol, sv, 1.0)
