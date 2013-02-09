# simple modal analysis

a = 0.1
n = 5
l = 5
m = 20

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

solfec = SOLFEC ('DYNAMIC', 1, 'out/modal')

bulk = BULK_MATERIAL (solfec,
                      model = 'KIRCHHOFF',
		      young = 15E9,
		      poisson = 0.2,
		      density = 2E3)

bod = BODY (solfec, 'FINITE_ELEMENT', mesh, bulk)

out = MODAL_ANALYSIS (bod, m, 'out/modal/data', verbose = 'ON')

print 'Eigenvalues:', out [0]
