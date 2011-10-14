# simple modal analysis

a = 0.1
n = 1

nodes = [-a, -a, -a,
          a, -a, -a,
          a,  a, -a,
         -a,  a, -a,
         -a, -a,  a,
          a, -a,  a,
          a,  a,  a,
         -a,  a,  a]

mesh = HEX (nodes, n, n, n, 0, [0, 1, 2, 3, 4, 5])

solfec = SOLFEC ('DYNAMIC', 1, 'out/modal')

bulk = BULK_MATERIAL (solfec,
                      model = 'KIRCHHOFF',
		      young = 1,
		      poisson = 0,
		      density = 1)

bod = BODY (solfec, 'FINITE_ELEMENT', mesh, bulk)

out = MODAL_ANALYSIS (bod, 3*(n+1)**3)
