# RT_SPLIT example: runtime mesh splitting from a callback routine
import os
d0 = os.path.dirname(os.path.realpath(__file__))

step = 1E-3 # time step
stop = 1.5 # duration

slv = NEWTON_SOLVER () # constraints solver
sol = SOLFEC ('DYNAMIC', step, 'out/runtime-split') # simulation object
bulk = BULK_MATERIAL (sol, model = 'KIRCHHOFF',
       young = 1E6, poisson = 0.3, density = 1E3) # elastic material
# contact material featuring non-penetration and Coulomb friction
SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)
GRAVITY (sol, (0, 0, -10)) # gravity acceleration

# pipe geometry routine
def pipe (point, name, tetgen = False):
  shp = PIPE (point, (0, 0, 2), 0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])
  if tetgen: shp = TETRAHEDRALIZE (shp, d0+'/'+name+'.mesh')
  ROTATE (shp, (0, 0, 1), (1, 0, 0), 90)
  bod = BODY (sol, 'FINITE_ELEMENT', shp, bulk, label = name)
  INITIAL_VELOCITY (bod, (0, 0, -1), (0, 0, 0))
  bod.selfcontact = 'ON' # detect self-contact
  bod.damping = 5E-3 # hand tuned
  if HERE(sol, bod):
    msh = bod.mesh
    for i in range (0, msh.nnod):
      # display points help to pick up
      # face node numbers in Soflec viewer
      DISPLAY_POINT (bod, msh.node (i), '%d'%i)

# create pipes
pipe ((0, 0, 0), 'pipe0')
pipe ((4, 0, 0), 'pipe1', True)

# create base disk
shp = PIPE  ((2, 0, -0.5), (0, 0, -0.1),
  0.05, 4.0, 1, 36, 1, 1, [1, 1, 1, 1, 1, 1])
BODY (sol, 'OBSTACLE', shp, bulk)

# pipe0 hex mesh faces picked visually in Solfec viewer
hf = [[0, 1, 17, 16], [2, 3, 19, 18], [4, 5, 21, 20],
      [6, 7, 23, 22], [8, 9, 25, 24], [10, 11, 27, 26],
      [12, 13, 29, 28], [14, 15, 31, 30]]

# three points defining a plane passing pipe1
p0 = (5.5, -1, 1) # picked visually in Solfec viewer
p1 = (5.5, 1, 1)
p2 = (4.5, 1, 1)
def sub(a, b): return (a[0]-b[0],a[1]-b[1],a[2]-b[2])
def prod(a, b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
n0 = prod(sub(p1, p0), sub(p2, p1)) # plane normal
def norm(a): return (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])**0.5
def scale(a, eps): return (a[0]*eps,a[1]*eps,a[2]*eps)
l0 = norm(n0)
n0 = scale(n0, 1.0/l0) # unit normal

# runtime callback
def callback (sol):
  bod = BYLABEL (sol, 'BODY', 'pipe0')
  # create first split (remains one body)
  if sol.time > stop*2.0/6.0 and bod != None and bod.mesh.nnod == 32:
    print 'Invoking 1st RT_SPLIT for pipe0'
    (bod2, lst1, lst2) = RT_SPLIT (bod, [hf[0]], 1, 2, 'pipe0/*')
    print 'bod2:', bod2
    print 'lst1:', lst1
    print 'lst2:', lst2
  # create second split (fragments in two bodies)
  bod = BYLABEL (sol, 'BODY', 'pipe0/*')
  if sol.time > stop*4.0/6.0 and bod != None and bod.mesh.nnod > 32:
    print 'Invoking 2nd RT_SPLIT for pipe0'
    ief = bod.mesh.inter_element_faces()
    (bod2, lst1, lst2) = RT_SPLIT (bod, [ief[1]], 1, 2,
                              'pipe0/1', 'pipe0/2')
    print 'bod2:', bod2
    print 'lst1:', lst1
    print 'lst2:', lst2
  # split second body in two fragemtns
  bod = BYLABEL (sol, 'BODY', 'pipe1')
  if sol.time > stop*2.0/6.0 and bod != None:
    faces = bod.mesh.inter_element_faces_on_plane (p0, n0)
    print faces
    print 'Invoking RT_SPLIT for pipe1'
    (bod2, lst1, lst2) = RT_SPLIT (bod, faces,
                          1, 2, 'pipe1/1', 'pipe1/2')
    print 'bod2:', bod2
    print 'lst1:', lst1
    print 'lst2:', lst2
  return 1

# set up callback routine
CALLBACK (sol, step, sol, callback)

# contact sparsification helps
# to avoid repeated contact points
CONTACT_SPARSIFY (sol, mindist=0.1)

# run simulation
RUN (sol, slv, stop)
