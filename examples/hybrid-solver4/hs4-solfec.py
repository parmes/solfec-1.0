import os
d0 = os.path.dirname(os.path.realpath(__file__))
execfile (d0 + '/hs4-globals.py')

percent = False # percentage progress
sweep = False # acceleration sweep
argv = NON_SOLFEC_ARGV()
if argv <> None and '-percent' in argv: percent = True
if argv <> None and '-sweep' in argv: sweep = True

if sweep:
  lofq = lofq_sweep
  hifq = hifq_sweep

nodes0 = [dO0, 0, 0,
          dO0+d01, 0, 0,
          dO0+d01, d12, 0,
	  dO0, d12, 0,
	  dO0, d12-d34, 0,
	  dO0+d45, d12-d34, 0,
	  dO0+d45, d12-d34-d47, 0,
	  dO0, d12-d34-d47, 0,
          dO0, 0, dOz,
          dO0+d01, 0, dOz,
	  dO0+d01, d12, dOz,
	  dO0, d12, dOz,
	  dO0, d12-d34, dOz,
	  dO0+d45, d12-d34, dOz,
	  dO0+d45, d12-d34-d47, dOz,
	  dO0, d12-d34-d47, dOz]


elements0 = [8, 0, 1, 6, 7, 8, 9, 14, 15, 0,
             8, 6, 1, 2, 5, 14, 9, 10, 13, 0,
	     8, 4, 5, 2, 3, 12, 13, 10, 11, 0]


mesh0 = MESH (nodes0, elements0, 0)

nodes1 = [0, d34+gap, 0,
          dO0+d45-gap, d34+gap, 0,
	  dO0+d45-gap, d34+d47-gap, 0,
	  0, d34+d47-gap, 0,
          0, d34+gap, dOz,
          dO0+d45-gap, d34+gap, dOz,
	  dO0+d45-gap, d34+d47-gap, dOz,
	  0, d34+d47-gap, dOz]

mesh1 = HEX(nodes1, 1, 1, 2, 0, [0]*6)

sol = SOLFEC ('DYNAMIC', step, 'out/hs4-solfec' + ('-sweep' if sweep else ''))
if percent: sol.verbose = '%'

mat = BULK_MATERIAL (sol, model = 'KIRCHHOFF',
    young = 1E9, poisson = 0.25, density = 1E3)

bod0 = BODY (sol, 'FINITE_ELEMENT', mesh0, mat)
bod0.scheme = 'DEF_LIM'
bod0.damping = step

bod1 = BODY (sol, 'RIGID', mesh1, mat)

FIX_DIRECTION (bod0, tuple(nodes0[1*3:1*3+3]), (1, 0, 0))
FIX_DIRECTION (bod0, tuple(nodes0[1*3:1*3+3]), (0, 0, 1))
FIX_DIRECTION (bod0, tuple(nodes0[9*3:9*3+3]), (1, 0, 0))
FIX_DIRECTION (bod0, tuple(nodes0[9*3:9*3+3]), (0, 0, 1))

(vt, vd, vv, va) = acc_sweep (step, stop, lofq, hifq, amag)
tsv = [None]*(len(vt)+len(vd))
tsv[::2] = vt
tsv[1::2] = vv
tsv = TIME_SERIES (tsv)

SET_VELOCITY (bod0, tuple(nodes0[1*3:1*3+3]), (0, 1, 0), tsv)
SET_VELOCITY (bod0, tuple(nodes0[9*3:9*3+3]), (0, 1, 0), tsv)

FIX_POINT (bod1, tuple(nodes1[0*3:0*3+3]))
FIX_POINT (bod1, tuple(nodes1[4*3:4*3+3]))

OUTPUT (sol, 1E-2)

RUN (sol, NEWTON_SOLVER(), stop)

if sol.mode == 'READ' and not VIEWER():
  try:
    import matplotlib.pyplot as plt
    dur = DURATION (sol)
    th = HISTORY (sol, [(bod0, tuple(nodes0[0:3]), 'DY')], dur [0], dur [1])
    plt.plot (th[0], th[1], label='baseline')
    plt.xlabel ('Time [s]')
    plt.ylabel ('DY of node 0 [m]')
    plt.legend(loc = 'upper right')
    plt.title('Baseline: Solfec only')
    plt.savefig (sol.outpath+'/node0dy.png')
    import pickle
    pickle.dump (th[0], open(sol.outpath+'/node0t.pickle', 'wb'))
    pickle.dump (th[1], open(sol.outpath+'/node0dy.pickle', 'wb'))
  except ImportError:
    pass # no reaction
