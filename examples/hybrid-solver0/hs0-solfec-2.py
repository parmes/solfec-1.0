# as in hs0-solfec-1.py
step = 0.01
sol = SOLFEC ('DYNAMIC', step, 'out/hybrid-solver0')
GRAVITY (sol, (0, 0, -10))
mat = BULK_MATERIAL (sol, model = 'KIRCHHOFF', young = 1, poisson = 0.25, density = 1)
SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.1)
nodes = [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1]
msh = HEX (nodes, 1, 1, 1, 0, [0, 1, 2, 3, 4, 5])
bod1 = BODY (sol, 'RIGID', msh, mat) # boundary bodies are rigid
msh = HEX (nodes, 2, 2, 2, 0, [0, 1, 2, 3, 4, 5])
TRANSLATE (msh, (0, 0, 1.1))
bod2 = BODY (sol, 'RIGID', msh, mat)
ns = NEWTON_SOLVER ()

# simulation end time
stop = 10.0

# parmec's output files are written to the same location as the input path;
# for that to be the solfec's output directory, we copy parmec's input file there
from shutil import copyfile
copyfile('examples/hybrid-solver0/hs0-parmec.py', 'out/hybrid-solver0/hs0-parmec.py')

# nubering of bodies in Parmec starts from 0 hence below we
# use dictionary {0 : bod1.id} as the parmec2solfec mapping
hs = HYBRID_SOLVER ('out/hybrid-solver0/hs0-parmec.py', step, {0 : 1}, ns)

import parmec as parmec # parmec module becomes available after creation of the HYBRID_SOLVER
tms0 = parmec.TSERIES ([0, 0.03, 5, 0.03, 5+step, 0.1, stop, 0.1])
hs.parmec_interval = (tms0, 0.03) # variable output file interval and constant output history interval

# parmec time history plots are collected at runtime
t_parmec = parmec.HISTORY ('TIME')
dz_parmec = parmec.HISTORY ('DZ', 0)

# solfec time history plots at runtime can be collected via a callback
t_solfec = []
dz_solfec = []
def plot_callback(sol, bod):
  t_solfec.append(sol.time)
  disp = DISPLACEMENT(bod, bod.center)
  dz_solfec.append(disp[2])
  return 1
CALLBACK (sol, 0.03, (sol, bod2), plot_callback)

import solfec as solfec # we need to be specific when using the OUTPUT command
solfec.OUTPUT (sol, 0.03) # since 'OUTPUT' in Solfec collides with 'OUTPUT' in Parmec

# run simulation
RUN (sol, hs, stop)

# plot time histories
try:
  import matplotlib.pyplot as plt
  print 'Plotting time histories...'
  plt.clf ()
  plt.plot (t_parmec, dz_parmec, label = 'lower body (parmec)')
  plt.plot (t_solfec, dz_solfec, label = 'upper body (solfec)')
  plt.legend (loc = 'upper right')
  plt.xlim ((0, stop))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('dz $(m)$')
  plt.savefig ('out/hybrid-solver0/hs0-dz.png')
except:
  print 'Plotting using matplotlib has failed.'

# XDMF export
if sol.mode == 'WRITE' and not VIEWER():
  print 'Run one more times to export XDMF files!'
elif sol.mode == 'READ' and not VIEWER():
  XDMF_EXPORT (sol, (0.0, stop), 'out/hybrid-solver0/hs0-solfec', [3, 4, 5])
