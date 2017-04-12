# import model preparation routines
import os
dirpath = os.path.dirname(os.path.realpath(__file__))
execfile (dirpath + '/ro0-prep.py')

# calculate POD base
sol = ro0_model (1.0/256, 0.0, 'TL', verbose='OFF', overwrite=True)
rig = RIGID_DISPLACEMENTS (sol.bodies[0])
defo = COROTATED_DISPLACEMENTS (sol, sol.bodies[0])
RUN (sol, NEWTON_SOLVER(), 1.0)
pod_base = ro0_POD_base (rig, defo)

# generate modal base
modal_base = ro0_modal_base ()

# Save POD and modal bases
import pickle
pickle.dump(pod_base, open('out/reduced-order0/pod_base.pickle', 'wb')) 
pickle.dump(modal_base, open('out/reduced-order0/modal_base.pickle', 'wb')) 
