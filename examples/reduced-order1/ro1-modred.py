import numpy
import modred
import pickle
import time
import os, sys
dirpath = os.path.dirname(os.path.realpath(__file__))
if modred.__version__[0] != '1':
  print 'modred version 1 is needed!'
  sys.exit(1)

try:
  podbase = pickle.load(open('out/reduced-order1/podbase.pickle', 'rb'))
except:
  t0 = time.time()
  num_modes = 18 # number of POD modes to use
  try:
    rig = pickle.load(open('out/reduced-order1/rig.pickle', 'rb'))
    dsp = pickle.load(open('out/reduced-order1/dsp.pickle', 'rb'))
  except:
    print 'File(s) out/reduced-order1/[rig,dsp].pickle not found',
    print '--> running ro0-fem-tl.py ...',
    global percentage
    percentage = True
    execfile (dirpath + '/ro1-fem-tl.py')
    try:
      rig = pickle.load(open('out/reduced-order1/rig.pickle', 'rb'))
      dsp = pickle.load(open('out/reduced-order1/dsp.pickle', 'rb'))
    except:
      print 'Running ro0-fem-tl.py has failed --> report a bug!' 
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
  basevec = [x for vec in mod for x in vec]
  podbase = (val[0:len(mod)], basevec)
  pickle.dump(podbase, open('out/reduced-order1/podbase.pickle', 'wb'))

  t1 = time.time()
  print '\bPOD runtime: %.3f seconds' % (t1-t0)
