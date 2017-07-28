# partially cached TIME_SERIES ('data.txt', cache = N) example

step = 1E-3
solfec = SOLFEC ('DYNAMIC', step, 'out/cached-time-series')
mat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 1E9, poisson = 0.3, density = 1E3)
sph = SPHERE ((0, 0, 0), 1, 1, 1)
bod = BODY (solfec, 'RIGID', sph, mat)
tms = TIME_SERIES ('data.txt', cache = 1000)
print tms.times
#SET_VELOCITY (bod, (0, 0, 0), (1, 0, 0), tms)
#slv = NEWTON_SOLVER()
#RUN (solfec, slv, 1.0)
