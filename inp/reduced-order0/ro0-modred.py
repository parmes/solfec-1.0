import numpy
import modred
import pickle

num_modes = 10
rig = pickle.load(open('out/reduced-order0/rig.pickle', 'rb'))
dsp = pickle.load(open('out/reduced-order0/dsp.pickle', 'rb'))
vecs = numpy.transpose(numpy.array(rig+dsp))
print 'Input vectors shape:', vecs.shape
print 'Calculating POD...'
modes, eig_vals = modred.compute_POD_matrices_snaps_method(vecs, list(range(num_modes)))
#modes, eig_vals = modred.compute_POD_matrices_direct_method(vecs, list(range(num_modes)))

print 'Modes:'
print modes
print 'Eig vals:'
print eig_vals

mod = numpy.transpose(modes).tolist()
val = eig_vals.tolist()
pickle.dump(mod, open('out/reduced-order0/mod.pickle', 'wb'))
pickle.dump(val, open('out/reduced-order0/val.pickle', 'wb'))
