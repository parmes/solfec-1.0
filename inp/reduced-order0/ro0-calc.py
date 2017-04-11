# import model preparation routines
import os
dirpath = os.path.dirname(os.path.realpath(__file__))
execfile (dirpath + '/ro0-prep.py')

# solver to pass to RUN
slv = NEWTON_SOLVER ()

if VIEWER():
  # create five basic simulation kinds
  # for viewing purposes --> TODO
  sol = ro0_model (1E-3, 0.0)
  RUN (sol, slv, 1.0)

else:
  # run initial 'TL' simulation to sample displacements
  sol = ro0_model (1E-3, 0.0)
  if sol.mode == 'WRITE':
    rig = RIGID_DISPLACEMENTS (sol.bodies[0])
    defo = COROTATED_DISPLACEMENTS (sol, sol.bodies[0])
  RUN (sol, slv, 1.0)
  ro_base = ro0_POD_base (rig, defo, verbose=True)
  # FIXME --> OUTPUT interval not specified and it seems
  # like deformational displacements are sampled just once
