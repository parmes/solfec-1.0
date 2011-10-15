# simple modal analysis

a = 0.1
n = 10
l = 1
m = 50

nodes = [-a, -a, -l*a,
          a, -a, -l*a,
          a,  a, -l*a,
         -a,  a, -l*a,
         -a, -a,  l*a,
          a, -a,  l*a,
          a,  a,  l*a,
         -a,  a,  l*a]

mesh = HEX (nodes, n, n, l*n, 0, [0, 1, 2, 3, 4, 5])

solfec = SOLFEC ('DYNAMIC', 1, 'out/modal')

bulk = BULK_MATERIAL (solfec,
                      model = 'KIRCHHOFF',
		      young = 15E9,
		      poisson = 0.2,
		      density = 2E3)

bod = BODY (solfec, 'FINITE_ELEMENT', mesh, bulk)

print 'Computing ', m, 'eigenpairs of a ', 3*(n+1)**3, 'system ...'

out = MODAL_ANALYSIS (bod, m)

print 'Eigenvalues:', out [0]
