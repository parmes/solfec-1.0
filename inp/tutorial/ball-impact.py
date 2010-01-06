# Ball impacting a plate

# domain creation
def ball_impact (step, stop, spring_value, dashpot_value, output):
  w = 2
  l = 2
  h = 1
  floor_vid = 1
  floor_sid = 1

  floor = HULL (
	 [-w/2, -l/2, -h,
	   w/2, -l/2, -h,
	   w/2,  l/2, -h,
	  -w/2,  l/2, -h,
	  -w/2, -l/2,  0,
	   w/2, -l/2,  0,
	   w/2,  l/2,  0,
	  -w/2,  l/2,  0], floor_vid, floor_sid)

  solfec = SOLFEC ('DYNAMIC', step, output)

  bulk = BULK_MATERIAL (solfec,
			model = 'KIRCHHOFF',
			young = 15E9,
			poisson = 0.3,
			density = 2E3)

  BODY (solfec, 'OBSTACLE', floor, bulk)

  sphere = SPHERE ((0, 0, 1.0), 1, 1, 1)
  body = BODY (solfec, 'RIGID', sphere, bulk)
  INITIAL_VELOCITY (body, (0, 0, -5), (0, 0, 0))

  SURFACE_MATERIAL (solfec, model = 'SPRING_DASHPOT', friction = 0.0, spring = spring_value, dashpot = dashpot_value)

  GRAVITY (solfec, (0, 0, -1), 10)

  xs = EXPLICIT_SOLVER ()

  RUN (solfec, xs, stop)

  return solfec

# main module
step = 1E-3
stop = 2.0
spring = 1E+9

sol1 = ball_impact (step, stop, spring, 0E0, 'out/tutorial/sphere-impact-1')
sol2 = ball_impact (step, stop, spring, 1E6, 'out/tutorial/sphere-impact-2')
sol3 = ball_impact (step, stop, spring, 1E7, 'out/tutorial/sphere-impact-3')

# plotting
if not VIEWER() and sol1.mode == 'READ':
  import matplotlib.pyplot as plt
  th = HISTORY (sol1, (sol1, 'KINETIC'), 0, stop)
  plt.plot (th [0], th [1], lw = 2, label='kin (0)')
  th = HISTORY (sol2, (sol2, 'KINETIC'), 0, stop)
  plt.plot (th [0], th [1], lw = 2, label='kin (1E6)')
  th = HISTORY (sol3, (sol3, 'KINETIC'), 0, stop)
  plt.plot (th [0], th [1], lw = 2, label='kin (1E7)')
  plt.axis (xmin = 0, xmax = 2, ymin=-10000, ymax=110000)
  plt.legend(loc = 'upper right')
  plt.savefig ('doc/figures/ball-impact.eps')
