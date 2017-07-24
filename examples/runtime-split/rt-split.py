# Runtime mesh split example

step = 1E-3
stop = 1.0

slv = GAUSS_SEIDEL_SOLVER (1E-4, 1000)
sol = SOLFEC ('DYNAMIC', step, 'out/runtime-split')
bulk = BULK_MATERIAL (sol, model = 'KIRCHHOFF',
                      young = 1E6, poisson = 0.3, density = 1E3)
SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)
GRAVITY (sol, (0, 0, -10))

# pipe and base geometry module routine
def pipe_and_base (point, name):
  shp = PIPE (point, (0, 0, 2), 0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])
  bod = BODY (sol, 'FINITE_ELEMENT', shp, bulk, label = name)
  INITIAL_VELOCITY (bod, (0, 0, -1), (0, 0, 0))
  bod.selfcontact = 'ON'
  bod.scheme = 'DEF_LIM'
  bod.damping = 5E-3
  if HERE(sol, bod):
    msh = bod.mesh
    for i in range (0, msh.nnod):
      DISPLAY_POINT (bod, msh.node (i), '%d'%i)
      # display points can help to pick split surface faces

  shp = PIPE  (TRANSLATE (point, (0, 0, -0.25)), (0, 0, -1),
               0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])
  BODY (sol, 'OBSTACLE', shp, bulk)

# create several pipe and base modules
pipe_and_base ((0, 0, 0), 'pipe0')
pipe_and_base ((4, 0, 0), 'pipe1')

# mesh faces picked visually in Solfec viewer
f = [[0, 1, 17, 16], [2, 3, 19, 18], [4, 5, 21, 20],
     [6, 7, 23, 22], [8, 9, 25, 24], [10, 11, 27, 26],
     [12, 13, 29, 28], [14, 15, 31, 30]]

# dictionray mapping of old-to-new nodes for pipe0
o2n0 = {}

# runtime callback
def callback (sol):
  bod = BYLABEL (sol, 'BODY', 'pipe0')
  # create first split (remains one body)
  if sol.time > stop*2.0/6.0 and bod != None and len(o2n0) == 0:
    print 'Invoking 1st RT_SPLIT for pipe0'
    (bod2, lst1, lst2) = RT_SPLIT (bod, [f[0]], 1, 2)
    print 'bod2:', bod2
    print 'lst1:', lst1
    print 'lst2:', lst2
    for i in range (0, len(lst1)):
      o2n0[lst1[i]] = i # don't mind overwrites:
                        # only the unique part is used
  # create second split (fragments in two bodies)
  if sol.time > stop*4.0/6.0 and bod != None:
    print 'Invoking 2nd RT_SPLIT for pipe0'
    f1 = []
    for n in f[3]: f1.append (o2n0[n])
    (bod2, lst1, lst2) = RT_SPLIT (bod, [f1], 1, 2, 'pipe0/1', 'pipe0/2')
    print 'bod2:', bod2
    print 'lst1:', lst1
    print 'lst2:', lst2

  # split second body in two fragemtns
  bod = BYLABEL (sol, 'BODY', 'pipe1')
  if sol.time > stop*3.0/6.0 and bod != None:
    print 'Invoking RT_SPLIT for pipe1'
    (bod2, lst1, lst2) = RT_SPLIT (bod, [f[2], f[6]], 1, 2, 'pipe1/1', 'pipe1/2')
    print 'bod2:', bod2
    print 'lst1:', lst1
    print 'lst2:', lst2

  return 1

CALLBACK (sol, step, sol, callback)
RUN (sol, slv, stop)
