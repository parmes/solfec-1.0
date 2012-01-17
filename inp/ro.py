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

print 'Computing ', m, 'eigenpairs of a ', 3 * mesh.nnod, 'system ...'

solfec = SOLFEC ('DYNAMIC', step, 'out/ro')

bulk = BULK_MATERIAL (solfec,
                      model = 'KIRCHHOFF',
		      young = 15E9,
		      poisson = 0.2,
		      density = 2E3)

bod = BODY (solfec, 'FINITE_ELEMENT', COPY (mesh), bulk)
data = MODAL_ANALYSIS (bod, m, 'out/ro/modal')
DELETE (solfec, bod)

for i in range (0, 5):
  bod = BODY (solfec, 'FINITE_ELEMENT', TRANSLATE (COPY (mesh), (0, i*2*a, 0)), bulk, form = 'RO', modal = data)
  INITIAL_VELOCITY (bod, (0, 0, 0), (1, 0, 0))

#sv = NEWTON_SOLVER (delta = 1E-4)
sv = GAUSS_SEIDEL_SOLVER (1E-4, 1000)

#GRAVITY (solfec, (0, 0, -10))

RUN (solfec, sv, 10)
