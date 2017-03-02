from shutil import copyfile

M = 5 # must be same as in hs2-paremc.py
N = 3 # must be same as in hs2-paremc.py
gap = 0.001 # must be same as in hs2-paremc.py
step = 5E-4
stop = 5 # must be <= stop in hs2-parmec.py

sol = SOLFEC ('DYNAMIC', step, 'out/hybrid-solver2')

mat = BULK_MATERIAL (sol, model = 'KIRCHHOFF', young = 1E6, poisson = 0.25, density = 100)

SURFACE_MATERIAL (sol, model = 'SIGNORINI_COULOMB', friction = 0.1)

nodes = [0.0, 0.0, 0.0,
         0.1, 0.0, 0.0,
	 0.1, 0.1, 0.0,
	 0.0, 0.1, 0.0,
	 0.0, 0.0, 0.1,
	 0.1, 0.0, 0.1,
	 0.1, 0.1, 0.1,
	 0.0, 0.1, 0.1]

iparmec = 0
parmec2solfec = {}
isolfec = []
for i in range (0,M+N+M):
  for j in range (0,M+N+M):
    if i >= M and j >= M and i < M+N and j < M+N: # inner Solfec bodies
      msh = HEX (nodes, 2, 2, 2, 0, [0, 1, 2, 3, 4, 5])
      TRANSLATE (msh, (i*(0.1+gap), j*(0.1+gap), 0))
      p1 = msh.node(0)
      p2 = msh.node(2)
      p3 = msh.node(8)
      bod = BODY (sol, 'FINITE_ELEMENT', msh, mat)
      bod.scheme = 'DEF_LIM' # semi-implicit time integration
      bod.damping = 1E-4 # damping out free vibrations
      isolfec.append(bod.id)
      FIX_DIRECTION (bod, p1, (0, 0, 1))
      FIX_DIRECTION (bod, p2, (0, 0, 1))
      FIX_DIRECTION (bod, p3, (0, 0, 1))
    else:
      if (i in [M-1,M+N] and j in range(M,M+N)) or \
	 (j in [M-1,M+N] and i in range(M,M+N)) or \
	 (i in [M-1,M+N] and j in [M-1,M+N]): # Solfec-Parmec boundary
	msh = HEX (nodes, 1, 1, 1, 0, [0, 1, 2, 3, 4, 5])
	TRANSLATE (msh, (i*(0.1+gap), j*(0.1+gap), 0))
	bod = BODY (sol, 'RIGID', msh, mat) # boundary bodies are rigid
	parmec2solfec[iparmec] = bod.id
      iparmec = iparmec + 1 # next parmec body

ns = NEWTON_SOLVER ()

# parmec's output files are written to the same location as the input path
# for that to be the solfec's output directory, we copy parmec's input file there
copyfile('examples/hybrid-solver2/hs2-parmec.py', 'out/hybrid-solver2/hs2-parmec.py')

hs = HYBRID_SOLVER ('out/hybrid-solver2/hs2-parmec.py', 1E-4, parmec2solfec, ns, 0.03)

import solfec as solfec # we need to be specific when using the OUTPUT command
solfec.OUTPUT (sol, 0.03) # since 'OUTPUT' in Solfec collides with 'OUTPUT' in Parmec

RUN (sol, hs, stop)

if sol.mode == 'READ' and not VIEWER():
  XDMF_EXPORT (sol, (0.0, stop), 'out/hybrid-solver2/hs2-solfec', isolfec)
