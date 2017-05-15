# import model library routines
import os, sys
dirpath = os.path.dirname(os.path.realpath(__file__))
execfile (dirpath + '/ro0-lib.py')

# find path to executable
def where(program):
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, program)):
      return path
  return None

if VIEWER():
  # create 'TL','BC','BC-RO', 'BC-MODAL'
  # simulations for viewing purposes
  ro0_model (1.0/256, 0.0, 'TL')
  ro0_model (1.0/256, 0.0, 'BC')

  # locate sofec executable
  path = where('solfec')
  if path == None:
    print 'ERROR: solfec executable not found in PATH!'
    print '       Download and compile solfec; add solfec directory to PATH variable;'
    sys.exit(1)

  # run solfec externally to generate POD and modal bases
  from subprocess import call
  call (['solfec', dirpath + '/ro0-gen-bases.py'])

  # load POD and modal bases
  import pickle
  pod_base = pickle.load(open('out/reduced-order0/pod_base.pickle', 'rb'))
  modal_base = pickle.load(open('out/reduced-order0/modal_base.pickle', 'rb'))

  # create reduced models
  ro0_model (1.0/256, 0.0, 'BC-RO', pod_base)
  ro0_model (1.0/256, 0.0, 'BC-MODAL', modal_base)
