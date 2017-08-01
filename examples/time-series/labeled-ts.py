# labeled TIME_SERIES ('data.txt', label = 'string') example
import os
d0 = os.path.dirname(os.path.realpath(__file__))

# create SOLFEC object and bulk material
solfec = SOLFEC ('DYNAMIC', 1E-3, 'out/labeled-ts')
mat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF',
      young = 1E9, poisson = 0.3, density = 1E3)

# create labeled time series object
tms = TIME_SERIES (d0+'/ts-data-1.txt', label = 'data-1')

# create 10 rigid spheres following along the y-direction
# the time series based displacement history; note that
# in this case, because the time series object was labeled,
# only a single copy of the time series object will be used
# internally - saving memory; should 'tms' be unlabeled 10
# separate data sets would be used;
for i in range (0, 10):
  sph = SPHERE ((i, 0, 0), 0.4, 1, 1)
  bod = BODY (solfec, 'RIGID', sph, mat)
  SET_DISPLACEMENT (bod, (i, 0, 0), (0, 1, 0), tms)

# create constraints solver and run simulation
slv = NEWTON_SOLVER()
RUN (solfec, slv, 1.0)
