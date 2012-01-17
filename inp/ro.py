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

bod = BODY (solfec, 'FINITE_ELEMENT', COPY (mesh), bulk, form = 'RO')
bod.scheme = 'DEF_LIM'
MODAL_ANALYSIS (bod, m)
INITIAL_VELOCITY (bod, (0, 0, 0), (1, 0, 0))

for i in range (1, 5): # run -np 4 to get malloc errors
  if 0:
    bod = BODY (solfec, 'FINITE_ELEMENT', TRANSLATE (COPY (mesh), (0, i, 0)), bulk, form = 'RO')
    bod.scheme = 'DEF_LIM'
    MODAL_ANALYSIS (bod, m)
    INITIAL_VELOCITY (bod, (0, 0, 0), (1, 0, 0))
  else:
    b = CLONE (bod, (0, i*2*a, 0))
    INITIAL_VELOCITY (b, (0, 0, 0), (1, 0, 0))

#sv = NEWTON_SOLVER ()
sv = GAUSS_SEIDEL_SOLVER (1E-4, 1000)

#GRAVITY (solfec, (0, 0, -10))

RUN (solfec, sv, 10)
