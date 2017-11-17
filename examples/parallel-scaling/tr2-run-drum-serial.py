# auxiliary serial multi-running script
import subprocess
import os
dirpath = os.path.dirname(os.path.realpath(__file__))

argv = [dirpath+'/rotating-drum.py -kifo RG -nsdl 0.0 -stop 1.0',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 1E-14 -stop 1.0',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 1E-13 -stop 1.0',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 1E-12 -stop 1.0',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 1E-11 -stop 1.0',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 1E-10 -stop 1.0',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 1E-9 -stop 1.0',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 1E-8 -stop 1.0',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 1E-7 -stop 1.0',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 1E-6 -stop 1.0',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 1E-5 -stop 1.0']

for args in argv:
  command = 'solfec ' + args
  print 'Running:', command
  process = subprocess.Popen(command, shell=True)
  process.wait()
