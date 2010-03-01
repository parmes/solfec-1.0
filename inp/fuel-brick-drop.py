# simple core model fuel brick drop
import sys
sys.path.append ('inp/cores/inc')
from simple_core_base import *

def scene_base (material, solfec):

  vertices = [1, 0, 0,
              1, 1, 0,
              0, 1, 0,
              0, 0, 0,
              1, 0, 1,
              1, 1, 1,
              0, 1, 1,
              0, 0, 1]

  faces = [4, 0, 1, 5, 4, 3,
           4, 1, 2, 6, 5, 3,
           4, 2, 3, 7, 6, 3,
           4, 3, 0, 4, 7, 3,
           4, 0, 3, 2, 1, 3,
           4, 4, 5, 6, 7, 3]

  outd = 0.4598
  margin = 0.05
  thick = 0.1
  lx = outd + (margin + thick)
  ly = outd + (margin + thick)
  shape = []

  cvx = CONVEX (vertices, faces, 3)
  scl = (lx,  ly,  thick)
  vec = (-lx/2, -ly/2, -thick)
  SCALE (cvx, scl)
  TRANSLATE (cvx, vec)
  shape.append (cvx)

  BODY (solfec, 'OBSTACLE', shape, material)

def scene_generate (material, solfec, kinematics):

  scene_base (material, solfec)

  shp = gcore_brick (0, 0, 0)

  vertices = [0, 0, 0,
              1, 0, 0,
              1, 2, 0,
              0, 2, 0,
              0, 0, 3,
              1, 0, 3,
              1, 2, 3,
              0, 2, 3]

  TRANSLATE (shp, (0, 0, 0.11))
  ROTATE (shp, (0, 0, 0.11), (1, 1, 1), 35)
  msh = ROUGH_HEX (shp, 2, 2, 2)
  bod = BODY (solfec , kinematics, shp, material, mesh=msh)
  INITIAL_VELOCITY (bod, (0, 0, -10), (0, 0, 0))
  if (kinematics == 'FINITE_ELEMENT'): bod.scheme = 'DEF_IMP'

def scene_run (solver, kinematics):

  step = 1E-2
  stop = 3.0

  solfec = SOLFEC ('DYNAMIC', step, 'out/fuel-brick-drop/' + kinematics)
  SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.7)
  bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3)
  GRAVITY (solfec, (0, 0, -10))
  scene_generate (bulkmat, solfec, kinematics)
  OUTPUT (solfec, 0.02)
  RUN (solfec, solver, stop)
  return (solfec, kinematics)

# main module
#import rpdb2; rpdb2.start_embedded_debugger('a')

gs = GAUSS_SEIDEL_SOLVER (1E-3, 1000, failure = 'CONTINUE', diagsolver = 'PROJECTED_GRADIENT')

pairs = []
pairs.append (scene_run (gs, 'FINITE_ELEMENT'))
pairs.append (scene_run (gs, 'PSEUDO_RIGID'))
pairs.append (scene_run (gs, 'RIGID'))

if not VIEWER() and pairs [0][0].mode == 'READ':

  for pair in pairs:
    solfec = pair [0]
    kinematics = pair [1]
    timers = ['TIMINT', 'CONUPD', 'CONDET', 'LOCDYN', 'CONSOL']
    dur = DURATION (solfec)
    th = HISTORY (solfec, timers, dur[0], dur[1])
    total = 0.0

    for i in range (0, len (timers)):
      sum = 0.0
      for tt in th [i+1]: sum += tt
      print kinematics, timers [i], 'TIME:', sum
      total += sum

    print kinematics, 'TOTAL TIME:', total
