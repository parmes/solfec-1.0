# display point test

solfec = SOLFEC ('DYNAMIC', 1E-3, 'out/deve/display-point')

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)

hex = HEX ([0, 0, 0,
	    1, 0, 0,
	    1, 1, 0,
	    0, 1, 0,
	    0, 0, 1,
	    1, 0, 1,
	    1, 1, 1,
	    0, 1, 1], 2, 2, 2, 2, [2, 2, 2, 2, 2, 2])

b = BODY (solfec, 'RIGID', hex, bulkmat)
DISPLAY_POINT (b, (0, 0, 0), '(0, 0, 0)')
DISPLAY_POINT (b, (1, 0, 0), '(1, 0, 0)')
DISPLAY_POINT (b, (1, 1, 0), '(1, 1, 0)')
DISPLAY_POINT (b, (0, 1, 0), '(0, 1, 0)')

print b.display_points

gs = GAUSS_SEIDEL_SOLVER (1E-3, 100)

RUN (solfec, gs, 1.0)
