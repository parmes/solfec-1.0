# pinned bar example

def pinned_bar_create (material, solfec):

  nodes = [-0.05, -0.05, 0.0,
            0.05, -0.05, 0.0,
            0.05,  0.05, 0.0,
	   -0.05,  0.05, 0.0,
	   -0.05, -0.05, 1.0,
	    0.05, -0.05, 1.0,
	    0.05,  0.05, 1.0,
	   -0.05,  0.05, 1.0]

  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, 0]

  point = (0, 0, 0.75)

  vector = (0, 1, 0)

  fix1 = (0, -0.05, 0.75)

  fix2 = (0, 0.05, 0.75)

  msh = MESH (nodes, elements, 0)

  ROTATE (msh, point, vector, -30)

  bod = BODY (solfec, 'RIGID', msh, material)

  FIX_POINT (bod, fix1)
  FIX_POINT (bod, fix2)

  return bod

# callback stop function
def callback (bod):
  if bod.conf == None: return 1 # in parallel bod.conf might be None for a not-owner processor
  elif bod.conf [9] <= 0:
    exact = 2.121
    error = abs (bod.velo [1] - exact) / exact;
    if (error < 0.01): print 'PASSED'
    else:
      print 'FAILED'
      print '(','Angular velocity y [mass center x = 0] was %.3f' % bod.velo [1], 'while the reference value is %.3f' % exact, ')'
    return 0
  else: return 1

# main module
# import rpdb2; rpdb2.start_embedded_debugger('a')

step = 0.001
stop = 0.5

solfec = SOLFEC ('DYNAMIC', step, 'out/tests/pinned-bar')
solfec.verbose = 'OFF'

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 200e9, poisson = 0.3, density = 1)

GRAVITY (solfec, (0, 0, -9.8))

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'EXIT')

bod = pinned_bar_create (bulkmat, solfec)

if not VIEWER(): CALLBACK (solfec, step, bod, callback)

if not VIEWER() and solfec.mode == 'READ': print '\nPrevious test results exist. Please "make del" and rerun tests'
else: RUN (solfec, gs, stop)
