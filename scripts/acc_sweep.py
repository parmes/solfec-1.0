# acceleration sine sweep signal generation
from math import sin, cos, pi

# acc_sweep: generate constant magnitude acceleration sine sweep signal
# ---------------------------------------------------------------------
# step - signal time step
# stop - duration of the signal
# lofq - low frequency for the sweep
# hifq - high frequency for the sweep
# amag - acceleration magnitude
# -------------------------------------------------------------------------------
# returned: (vt, vd, vv, va), where
#         - vt is a list of time instants
#         - vd is a list of displacement values
#         - vv is a list of velocity values
#         - va is a list of acceleration values, at those time instants
# -------------------------------------------------------------------------------
def acc_sweep (step, stop, lofq, hifq, amag):

  t = 0.0
  v = 0.0
  va = []
  vv = []
  vf = []
  extend = 0.0
  while t < stop+extend:
    x = t + step/2. # mid-step time
    a = amag * sin (2*pi*(lofq+(hifq-lofq)*x/stop)*x) # mid-step acceleration
    v = v + a * step # mid-step integration of dv / dt = a into v
    va.append (a)
    vv.append (v)
    vf.append (lofq + (hifq-lofq)*(t/stop))
    if extend == 0.0 and len(vv) > 2 and vv[-1] < vv[-2]: extend = t
    t += step

  # find stabilized velocity level
  # by avaraging the last 5 minima and maxima
  imax = 0.0
  vmax = 0.0
  imin = 0.0
  vmin = 0.0
  i = len(vv)-2
  while i > 0:
    if vv[i-1] < vv[i] and vv[i] > vv[i+1]:
      imax = imax + 1.0
      vmax = vmax + vv[i]
    if vv[i-1] > vv[i] and vv[i] < vv[i+1]:
      imin = imin + 1.0
      vmin = vmin + vv[i]
    if imax == 5.0 and imin == 5.0: break
    i = i - 1
  vlevel = 0.1*(vmax+vmin)

  # find when this level is crossed from the start
  i = 0
  while vv[i] < vlevel: i = i + 1

  # trim histories to this moment
  while i > 0:
    va.pop(0)
    vv.pop(0)
    vf.pop(0)
    i -= 1

  # now produce displacement and time history
  vt = []
  vd = []
  d = 0.0
  t = 0.0
  for v in vv:
   vt.append (t)
   vd.append (d)
   t = t + step
   d = d + v * step  # integration of dd / dt = v

  # displacement has positive drift => find tangens of the positive drift angle 
  i = len(vd)-1
  while vd[i-1] > vd[i]: i -= 1 # first maximum
  while vd[i-1] < vd[i]: i -= 1 # previous minimum
  j = i
  while vd[j-1] > vd[i]: j += 1 # previous maximum

  # shift velocity down by the tangens of the drift angle
  vshift = (vd[i]+vd[j]) / (vt[i]+vt[j])
  for i in range (0, len(vv)): vv[i] -= vshift

  # after velocity has been shifted down, produce displacement envelope
  vd = []
  d = 0.0
  for v in vv:
   d = d + v * step  # integration of dd / dt = v
   vd.append (d)

  return (vt, vd, vv, va)
