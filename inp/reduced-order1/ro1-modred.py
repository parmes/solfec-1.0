import numpy
import modred
import pickle
import time

try:
  mod = pickle.load(open('out/reduced-order1/mod.pickle', 'rb'))
  val = pickle.load(open('out/reduced-order1/val.pickle', 'rb'))
except:
  t0 = time.time()
  num_modes = 10
  try:
    rig = pickle.load(open('out/reduced-order1/rig.pickle', 'rb'))
    dsp = pickle.load(open('out/reduced-order1/dsp.pickle', 'rb'))
  except:
    print 'File out/reduced-order1/modal.h5 not found',
    print '--> run ro0-fem-tl.py example first!'
    import sys
    sys.exit(0)
  vecs = numpy.transpose(numpy.array(rig+dsp))
  svec = vecs.shape[0]
  nvec = vecs.shape[1]
  print 'Solving for', nvec, 'input vectors of size', svec, '...'
  modes, eig_vals = modred.compute_POD_matrices_snaps_method(vecs, list(range(num_modes)))
  #modes, eig_vals = modred.compute_POD_matrices_direct_method(vecs, list(range(num_modes)))

  print 'Eig vals:'
  print eig_vals

  mod = numpy.transpose(modes).tolist()
  val = eig_vals.tolist()
  pickle.dump(mod, open('out/reduced-order1/mod.pickle', 'wb'))
  pickle.dump(val, open('out/reduced-order1/val.pickle', 'wb'))

  t1 = time.time()
  print '\bTotal runtime: %.3f seconds' % (t1-t0)
