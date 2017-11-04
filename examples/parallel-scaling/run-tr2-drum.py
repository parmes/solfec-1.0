# auxiliary serial multi-running script
import subprocess
import os
dirpath = os.path.dirname(os.path.realpath(__file__))

'''
argv = [dirpath+'/rotating-drum.py -kifo RG -nsdl 0.0',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 1E-8',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 2E-8',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 4E-8',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 8E-8',
        dirpath+'/rotating-drum.py -kifo RG -nsdl 1E-7',
        dirpath+'/rotating-drum.py -kifo PR -nsdl 0.0',
        dirpath+'/rotating-drum.py -kifo PR -nsdl 1E-8',
        dirpath+'/rotating-drum.py -kifo PR -nsdl 2E-8',
        dirpath+'/rotating-drum.py -kifo PR -nsdl 4E-8',
        dirpath+'/rotating-drum.py -kifo PR -nsdl 8E-8',
        dirpath+'/rotating-drum.py -kifo PR -nsdl 1E-7']
'''

argv = [dirpath+'/rotating-drum.py -kifo RG -nsep 0.75',
        dirpath+'/rotating-drum.py -kifo RG -nsep 0.5',
        dirpath+'/rotating-drum.py -kifo RG -nsep 0.25',
        dirpath+'/rotating-drum.py -kifo RG -nsep 0.2',
        dirpath+'/rotating-drum.py -kifo RG -nsep 0.1',
        dirpath+'/rotating-drum.py -kifo RG -nsep 0.05']

for args in argv:
  command = 'solfec ' + args
  print 'Running:', command
  process = subprocess.Popen(command, shell=True)
  process.wait()
