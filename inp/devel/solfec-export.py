# Set up domino toppling example
solfec = SOLFEC ('DYNAMIC', 1E-3, 'out/domino')
GRAVITY (solfec, (0, 0, -9.81))
mat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF',
    young = 15E9, poisson = 0.3, density = 1.8E3)
SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB',
                   friction = 0.5, restitution = 0.25)
cube = HEX ([0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1,
            1, 0, 1, 1, 1, 1, 0, 1, 1], 2, 2, 2, 1, [1]*6)
BODY (solfec, 'OBSTACLE', SCALE(COPY(cube), (1, 1, 0.1)), mat)
for i in range (0, 4):
  piece = SCALE (COPY(cube), (0.2, 0.05, 0.4))
  TRANSLATE (piece, (0.4, i*0.2, 0.1))
  BODY (solfec, 'RIGID', piece, mat, label = 'Domino' + str(i+1))
ball = BODY (solfec, 'RIGID', SPHERE ((0.5, -0.5, 0.4), 0.1, 3, 3), mat)
INITIAL_VELOCITY (ball, (0, 3, 0), (0, 0, 0))

# Export initial state or run simulation
argv = NON_SOLFEC_ARGV()
if argv != None and '--geom0' in argv:
  SOLFEC_EXPORT (solfec, 0.0, 'out/sxptest0')
  solfec.cleanup = 'ON'
else: RUN (solfec, NEWTON_SOLVER(), 1.0)

# Export results
if solfec.mode == 'READ' and not VIEWER():
  # export simulation state at t = 0.5
  SOLFEC_EXPORT (solfec, 0.5, 'out/sxptest1')
  # export entire simulation:
  SOLFEC_EXPORT (solfec, (0.0, 1.0), 'out/sxptest2')
  # 101 time instants:
  times = [0.01*i for i in range(0, 101)]
  # export all bodies at 101 times:
  SOLFEC_EXPORT (solfec, times, 'out/sxptest3')
  # export a subset of bodies at 101 times:
  SOLFEC_EXPORT (solfec, times, 'out/sxptest4',
    subset = [(0.4, 0, 0, 0.6, 0.05, 0.2), 'Domino2', ball])
  # export Domino2 and Domino3 at 101 times:
  SOLFEC_EXPORT (solfec, times, 'out/sxptest5', subset = 'Domino[2,3]')
