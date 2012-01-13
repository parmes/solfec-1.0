# reduced order model test

a = 0.1
n = 1
l = 4
m = 24
step = 0.01

nodes = [-a, -a, -l*a,
          a, -a, -l*a,
          a,  a, -l*a,
         -a,  a, -l*a,
         -a, -a,  l*a,
          a, -a,  l*a,
          a,  a,  l*a,
         -a,  a,  l*a]

mesh = HEX (nodes, n, n, l*n, 0, [0, 1, 2, 3, 4, 5])
nod = mesh.node (0)

print 'Computing ', m, 'eigenpairs of a ', 3 * mesh.nnod, 'system ...'

solfec = SOLFEC ('DYNAMIC', step, 'out/ro')

bulk = BULK_MATERIAL (solfec,
                      model = 'KIRCHHOFF',
		      young = 15E9,
		      poisson = 0.2,
		      density = 2E3)

bod = BODY (solfec, 'FINITE_ELEMENT', mesh, bulk, form = 'RO')
bod.scheme = 'DEF_LIM'

MODAL_ANALYSIS (bod, m)

ns = NEWTON_SOLVER ()

#GRAVITY (solfec, (0, 0, -10))

INITIAL_VELOCITY (bod, (0, 0, 0), (1, 0, 0))

RUN (solfec, ns, 10)
