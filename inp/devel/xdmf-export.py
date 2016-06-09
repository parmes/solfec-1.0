# XDMF export example; domino toppling

step = 1E-3

solfec = SOLFEC ('DYNAMIC', step, 'out/domino')

GRAVITY (solfec, (0, 0, -9.81))

table = HULL ([0, 0, 0,
               0, 1, 0,
	       1, 1, 0,
	       1, 0, 0,
               0, 0, -0.1,
               0, 1, -0.1,
	       1, 1, -0.1,
	       1, 0, -0.1], 1, 1)

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5, restitution = 0.25)

BODY (solfec, 'OBSTACLE', table, bulkmat)

hex = HEX ([0, 0, 0,
	    1, 0, 0,
	    1, 1, 0,
	    0, 1, 0,
	    0, 0, 1,
	    1, 0, 1,
	    1, 1, 1,
	    0, 1, 1], 2, 2, 2, 2, [2, 2, 2, 2, 2, 2])

SCALE (hex, (0.2, 0.05, 0.4))

TRANSLATE (hex, (0.4, 0, 0))

for i in range (0, 4):
  shp = COPY (hex)
  TRANSLATE (shp, (0, i * 0.2, 0))
  b = BODY (solfec, 'RIGID', shp, bulkmat, label = 'Domino' + str(i+1))

shp = SPHERE ((0.5, -0.5, 0.3), 0.1, 3, 3)

ball = BODY (solfec, 'RIGID', shp, bulkmat)

INITIAL_VELOCITY (ball, (0, 3, 0), (0, 0, 0))

sv = NEWTON_SOLVER ()

OUTPUT (solfec, step)

RUN (solfec, sv, 1.0)

if solfec.mode == 'READ':
  times = [0.01*i for i in range(0, 101)]
  XDMF_EXPORT (solfec, 0.0, 'out/xmftest0') # initial state
  XDMF_EXPORT (solfec, (0.0, 1.0), 'out/xmftest1') # entire simulation
  XDMF_EXPORT (solfec, times, 'out/xmftest2') # all bodies at 101 times
  XDMF_EXPORT (solfec, times, 'out/xmftest3', subset = [(0.4, 0, 0, 0.6, 0.05, 0.04), 'Domino2', ball]) # a subset of bodies at 101 times
  XDMF_EXPORT (solfec, times, 'out/xmftest4', subset = 'Domino[2,3]', attributes = ['VELO', 'RELV']) # Domino2 and Domino3 with two attributes at 101 times
