# RT_SPLIT test
step = 1E-3
tet = True
single = True

sv = GAUSS_SEIDEL_SOLVER (1E-4, 1000)

sol = SOLFEC ('DYNAMIC', step, 'out/rt_split')

bulk = BULK_MATERIAL (sol,
                      model = 'KIRCHHOFF',
		      young = 1E6,
		      poisson = 0.3,
		      density = 1E3)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.3)

GRAVITY (sol, (0, 0, -10))

def pipe_and_base (point):
  shp = PIPE (point, (0, 0, 2), 0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])
  if tet:
    msh = TETRAHEDRALIZE (shp, 'out/pipe.mesh')
  else: msh = shp
  print 'inter-element-faces:', msh.inter_element_faces()
  bod = BODY (sol, 'FINITE_ELEMENT', msh, bulk, label = 'pipe')
  #if single: bod.selfcontact = 'ON'
  if HERE(sol, bod):
    msh = bod.mesh
    for i in range (0, msh.nnod):
      DISPLAY_POINT (bod, msh.node (i), '%d'%i)

  shp = PIPE  (TRANSLATE (point, (0, 0, -1)), (0, 0, -1), 0.5, 1, 1, 8, 1, 1, [1, 1, 1, 1, 1, 1])
  BODY (sol, 'OBSTACLE', shp, bulk)

pipe_and_base ((0, 0, 0))

hf = [[0, 1, 17, 16], [2, 3, 19, 18], [4, 5, 21, 20],
      [6, 7, 23, 22], [8, 9, 25, 24], [10, 11, 27, 26],
      [12, 13, 29, 28], [14, 15, 31, 30]]

tf = [[2, 3, 5], [2, 5, 7], [17, 19, 16], [19, 18, 16]]

def callback (sol):
  bod = BYLABEL (sol, 'BODY', 'pipe')
  if sol.time > 0.5 and bod != None:
    print 'Invoking RT_SPLIT at time 0.5s'
    if tet and single: (bod2, lst1, lst2) = RT_SPLIT (bod, [tf[0], tf[1]], 1, 2, 'frag1', 'frag2')
    elif tet: (bod2, lst1, lst2) = RT_SPLIT (bod, tf, 1, 2, 'frag1', 'frag2')
    elif single: (bod2, lst1, lst2) = RT_SPLIT (bod, [hf[0]], 1, 2, 'frag1', 'frag2')
    else: (bod2, lst1, lst2) = RT_SPLIT (bod, [hf[0], hf[1]], 1, 2, 'frag1', 'frag2')
    print 'bod2:', bod2
    print 'lst1:', lst1
    print 'lst2:', lst2
  return 1

CALLBACK (sol, step, sol, callback)

RUN (sol, sv, 1.0)

if sol.mode == 'READ':
  dur = DURATION (sol)
  pnt = (1.06066, 1.06066, 0)
  th = HISTORY (sol, ('frag1',pnt, 'DZ'), dur[0], dur[1])
  print th
