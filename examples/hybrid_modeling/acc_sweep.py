# acceleration sine sweep signal generation
import matplotlib.pyplot as plt
from math import sin, cos, pi

# acc_sweep: generate constant magnitude acceleration sine sweep signal
# ---------------------------------------------------------------------
# step - signal time step
# stop - duration of the signal
# lofq - low frequency for the sweep
# hifq - high frequency for the sweep
# amag - acceleration magnitude
# acc_plot - path to acceleration signal plot
# vel_plot - path to velocity signal plot
# dsp_plot - path to displacement signal plot
# -------------------------------------------------------------------------------
# returned - [t0,v0,t1,v1,...] linear spline of the corresponding velocity signal
# -------------------------------------------------------------------------------
def acc_sweep (step, stop, lofq, hifq, amag, acc_plot = None, vel_plot = None, dsp_plot = None):

  t0 = 0.0
  while t0 < stop:
    a0 = amag * sin (2*pi*(lofq+(hifq-lofq)*t0/stop)*t0)
    t0 += step
    a1 = amag * sin (2*pi*(lofq+(hifq-lofq)*t0/stop)*t0)
    if a1 < a0: break # find first acceleration maximum
                      # (integrated velocity will be more symmetrical about zero)
  t = t0
  v = 0.0
  vt = []
  va = []
  vv = []
  vf = []
  while t < stop+t0:
    x = t + step/2. # mid-step time
    a = amag * sin (2*pi*(lofq+(hifq-lofq)*x/stop)*x) # mid-step acceleration
    v = v + a * step # mid-step integration of dv / dt = a into v
    vt.append (t-t0)
    va.append (a)
    vv.append (v)
    vf.append (lofq + (hifq-lofq)*(t/stop))
    t += step

  for i in range (0, len(vv)):
    if vv[i+1] < vv[i]: break # find first velocity maximum

  for j in range (0, len(vv)-i):
    vv[j] = vv[j+i] # shift velocity so it starts from the first maximum

  while i > 0:
    vt.pop() # remove last i items from lists
    va.pop()
    vv.pop()
    vf.pop()
    i -= 1

# after velocity has been trimmed, produce displacement envelope
  vd = []
  d = 0.0
  for v in vv:
   d = d + v * step  # integration of dd / dt = v
   vd.append (d)

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

  if acc_plot != None:
    plt.clf ()
    plt.plot (vt, va)
    plt.xlim ((0, stop))
    plt.xlabel ('time $(s)$')
    plt.ylabel ('acceleration $(m/s^2)$')
    plt.savefig (acc_plot)

  if vel_plot != None:
    plt.clf ()
    plt.plot (vt, vv)
    plt.xlim ((0, stop))
    plt.xlabel ('time $(s)$')
    plt.ylabel ('velocity $(m/s)$')
    plt.savefig (vel_plot)

  if dsp_plot != None:
    plt.clf ()
    plt.plot (vf, vd)
    plt.xlim ((lofq, hifq))
    plt.xlabel ('Frequency $(Hz)$')
    plt.ylabel ('displacement $(m)$')
    plt.savefig (dsp_plot)

  # generate return linear spline list
  data = []
  for (t, v) in zip (vt, vv): data += [t, v] # [t0,v0,t1,v1,...]
  return data
