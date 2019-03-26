from shutil import copy
import fileinput
import subprocess
import os, sys

class Job:
  def __init__ (self, name, command):
    self.name = name
    self.command = command

leps = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 0.25]
lmxi = [3, 5, 10, 15, 20, 25, 30, 40, 50, 70, 100]
nsdl = [0.0, 0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5]
kifo = ['BC', 'RG']
lepsC = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01]
lmxiC = [  5,    5,   5,    5,     5,    5,  10,   10,  10,   10,    10,   10,  15,   15,  15,   15,    15,   15,  20,   20,  20,   20,    20,   20]
jobs = []

for k in range (0, len(kifo)):
  for i in range(0, len(leps)):
    for j in range(0, len(nsdl)):
      jobs.append (Job('0a%d_%d_%d' % (i,j,k), 'solfec examples/parallel-scaling/array-of-cubes.py -subd reldelta -outi 0.3 -kifo %s -lmxi 1000 -leps %g -nsdl %g -rldl avgWii -prfx A' % (kifo[k], leps[i], nsdl[j])))

for k in range (0, len(kifo)):
  for i in range(0, len(lmxi)):
    for j in range(0, len(nsdl)):
      jobs.append (Job('0b%d_%d_%d' % (i,j,k), 'solfec examples/parallel-scaling/array-of-cubes.py -subd reldelta -outi 0.3 -kifo %s -lmxi %d -leps 1E-3 -nsdl %g -rldl avgWii -prfx B' % (kifo[k], lmxi[i], nsdl[j])))

for k in range (0, len(kifo)):
  for i in range(0, len(lepsC)):
    for j in range(0, len(nsdl)):
      jobs.append (Job('0c%d_%d_%d' % (i,j,k), 'solfec examples/parallel-scaling/array-of-cubes.py -subd reldelta -outi 0.3 -kifo %s -lmxi %d -leps %g -nsdl %g -rldl avgWii -prfx C' % (kifo[k], lmxiC[i], lepsC[i], nsdl[j])))

if '--post' in sys.argv: # delete recreated file
  cmd = 'rm -f out/array-of-cubes/reldelta/TIMINGS'
  print cmd
  process = subprocess.Popen(cmd, shell=True)
  process.wait()

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

if '--post' in sys.argv: # copy stats to renamed files
  cmd = 'cp out/array-of-cubes/reldelta/ITERS tr2-cubes25-reldelta-iters'
  print cmd
  process = subprocess.Popen(cmd, shell=True)
  process.wait()
  cmd = 'cp out/array-of-cubes/reldelta/RUNTIMES tr2-cubes25-reldelta-runtimes'
  print cmd
  process = subprocess.Popen(cmd, shell=True)
  process.wait()
  cmd = 'cp out/array-of-cubes/reldelta/TIMINGS tr2-cubes25-reldelta-timings'
  print cmd
  process = subprocess.Popen(cmd, shell=True)
  process.wait()
