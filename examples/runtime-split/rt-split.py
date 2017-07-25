# RT_SPLIT example: runtime mesh splitting from a callback routine

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
def pipe (point, name):
  shp = PIPE (point, (0, 0, 2), 0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])
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
pipe ((4, 0, 0), 'pipe1')

# create base disk
shp = PIPE  ((2, 0, -0.5), (0, 0, -0.1),
  0.05, 4.0, 1, 36, 1, 1, [1, 1, 1, 1, 1, 1])
BODY (sol, 'OBSTACLE', shp, bulk)

# mesh faces picked visually in Solfec viewer
f = [[0, 1, 17, 16], [2, 3, 19, 18], [4, 5, 21, 20],
     [6, 7, 23, 22], [8, 9, 25, 24], [10, 11, 27, 26],
     [12, 13, 29, 28], [14, 15, 31, 30]]

# runtime callback
def callback (sol):
  bod = BYLABEL (sol, 'BODY', 'pipe0')
  # create first split (remains one body)
  if sol.time > stop*2.0/6.0 and bod != None and bod.mesh.nnod == 32:
    print 'Invoking 1st RT_SPLIT for pipe0'
    (bod2, lst1, lst2) = RT_SPLIT (bod, [f[0]], 1, 2, 'pipe0/*')
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
    print 'Invoking RT_SPLIT for pipe1'
    (bod2, lst1, lst2) = RT_SPLIT (bod, [f[1], f[5]],
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
