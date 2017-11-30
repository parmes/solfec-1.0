from shutil import copy
import fileinput
import subprocess
import os, sys

class Job:
  def __init__ (self, name, command):
    self.name = name
    self.command = command

dens = [10., 50., 100., 500., 1000., 5000., 10000.]
nsdl = [0.0, 0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5]
jobs = []

for i in range(0,len(dens)):
  for j in range(0, len(nsdl)):
    jobs.append (Job('0dns_%d_%d' % (i, j), 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo RG -dens %g -lmxi 20 -leps 0.25 -nsdl %g -rldl avgWii -solv NS' % (dens[i], nsdl[j])))

for i in range(0,len(dens)):
  for j in range(0, len(nsdl)):
    jobs.append (Job('1dns_%d_%d' % (i, j), 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo PR -dens %g -lmxi 10 -leps 0.025 -nsdl %g -rldl avgWii -solv NS' % (dens[i], nsdl[j])))

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
