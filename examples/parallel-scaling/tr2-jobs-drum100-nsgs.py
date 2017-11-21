from shutil import copy
import fileinput
import subprocess
import os, sys

class Job:
  def __init__ (self, name, command):
    self.name = name
    self.command = command

jobs = []
jobs.append (Job('dr0_e_ns', 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS'))
jobs.append (Job('dr0_e_gs', 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS'))
jobs.append (Job('dr1_e_ns', 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS'))
jobs.append (Job('dr1_e_gs', 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS'))
jobs.append (Job('dr0_s_ns', 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS -sphs'))
jobs.append (Job('dr0_s_gs', 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS -sphs'))
jobs.append (Job('dr1_s_ns', 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS -sphs'))
jobs.append (Job('dr1_s_gs', 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS -sphs'))

for job in jobs: # schedule jobs
  print '***'
  print '*** scheduling: %s' % job.name
  print '***'
  copy ('examples/parallel-scaling/run.sh.serial', 'run.sh')
  for line in fileinput.input('run.sh', inplace=True):
    if 'solfec input-file.py' in line:
      print job.command
    else: print line,

  '''
  if '--post' in sys.argv:
    print job.command
    os.system('read -s -n 1 -p "Press any key to continue..."')
    print
  else:
    print 'sbatch -J %s run.sh' % job.name
    os.system('read -s -n 1 -p "Press any key to continue..."')
    print
  '''

  if '--post' in sys.argv: process = subprocess.Popen(job.command, shell=True)
  else: process = subprocess.Popen('sbatch -J %s run.sh' % job.name, shell=True)
  process.wait()
