# a block moved by a set velocity signal;
# the idea is to test whether there are issues with large input signals (>50MB files);
# run:
# python generate.py
# in this directory to generate 'series.txt' (90MB) before trying to use this example
# ----------
# 22/02/2015


step = 1E-3

solfec = SOLFEC ('DYNAMIC', step, 'out/devel/long-series')

GRAVITY (solfec, (0, 0, -9.81))

bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)

surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5, restitution = 0.25)

hex = HEX ([0, 0, 0,
	    1, 0, 0,
	    1, 1, 0,
	    0, 1, 0,
	    0, 0, 1,
	    1, 0, 1,
	    1, 1, 1,
	    0, 1, 1], 2, 2, 2, 2, [2, 2, 2, 2, 2, 2])

bod = BODY (solfec, 'RIGID', hex, bulkmat)

SET_VELOCITY(bod, (0, 0, 0), (0, 0, 1), 0)
SET_VELOCITY (bod, (1, 0, 0), (0, 0, 1), 0)
SET_VELOCITY (bod, (1, 1, 0), (0, 0, 1), 0)
SET_VELOCITY (bod, (0, 0, 0), (0, 1, 0), 0)
SET_VELOCITY (bod, (1, 0, 0), (0, 1, 0), 0)
ts = TIME_SERIES ('inp/devel/long-input-motion/series.txt');
SET_VELOCITY (bod, (1, 1, 0), (1, 0, 0), ts)

gs = GAUSS_SEIDEL_SOLVER (1E-3, 100)

OUTPUT (solfec, step)

RUN (solfec, gs, 1.0)
