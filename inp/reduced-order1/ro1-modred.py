import numpy
import modred
import pickle
import time

t0 = time.time()
num_modes = 10
rig = pickle.load(open('out/reduced-order1/rig.pickle', 'rb'))
dsp = pickle.load(open('out/reduced-order1/dsp.pickle', 'rb'))
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
