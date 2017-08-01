# partially cached TIME_SERIES ('data.txt', cache = N) example
import os
d0 = os.path.dirname(os.path.realpath(__file__))

# create SOLFEC object and bulk material
solfec = SOLFEC ('DYNAMIC', 1E-3, 'out/cached-ts')
mat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF',
      young = 1E9, poisson = 0.3, density = 1E3)

# create 10 rigid spheres following along the y-direction
# the time series based displacement history; note that
# in this case we use a unique time series for each body;
# should the number of bodies and the size of the time series
# be large - this would easily lead to using up all memory;
# using partially cached time series allows to avoid this
# issue; in our case we have 100 data points per series and
# we set the partial cache size to 10; this means that only
# 10 points are stored in memory, per series, at any time;
for i in range (1, 11):
  tms = TIME_SERIES (d0+'/ts-data-%d.txt' % i, cache = 10)
  sph = SPHERE ((i, 0, 0), 0.4, 1, 1)
  bod = BODY (solfec, 'RIGID', sph, mat)
  SET_DISPLACEMENT (bod, (i, 0, 0), (0, 1, 0), tms)

# create constraints solver and run simulation
slv = NEWTON_SOLVER()
RUN (solfec, slv, 1.0)
