step = 1E-4
stop = 5.0
nbod = 10

solfec = SOLFEC ('DYNAMIC', step, 'out/put_spring_mpi')
#solfec.verbose = 'OFF'

bulk = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E9, poisson = 0.3, density = 1E3)

volid = 0

def springfunc (stroke, velocity):
  return -1E6*stroke -1E3*velocity

prevbod = None

for i in range (0, nbod):
  nodes = [i, i, i,
	   i+1, i, i,
	   i+1, i+1, i,
	   i, i+1, i,
	   i, i, i+1,
	   i+1, i, i+1,
	   i+1, i+1, i+1,
	   i, i+1, i+1]

  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, volid]

  colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

  mesh = MESH (nodes, elements, colors)

  body = BODY (solfec, 'RIGID', mesh, bulk)

  if i: 
    PUT_SPRING (body, (i, i, i), prevbod, (i, i, i), springfunc, (-0.05, 0.05))
    CONTACT_EXCLUDE_BODIES (body, prevbod)

  prevbod = body

FIX_POINT (body, (i+1, i+1, i+1))

GRAVITY (solfec, (0, 0, -10))

#slv = NEWTON_SOLVER () # linearised calculation of spring values
slv = GAUSS_SEIDEL_SOLVER (1E-3, 100) # direct calculation of spring values

OUTPUT (solfec, 0.01)

RUN (solfec, slv, stop)
