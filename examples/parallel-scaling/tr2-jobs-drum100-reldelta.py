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
rldl = ['minWii', 'avgWii', 'maxWii']
lepsC = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01]
lmxiC = [  5,    5,   5,    5,     5,    5,  10,   10,  10,   10,    10,   10,  15,   15,  15,   15,    15,   15,  20,   20,  20,   20,    20,   20]
jobs = []

for i in range(0,len(leps)):
  for j in range(0, len(nsdl)):
    for k in range (0, len(rldl)):
      jobs.append (Job('0a%d_%d_%d' % (i,j,k), 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo RG -lmxi 1000 -leps %g -nsdl %g -rldl %s -prfx A' % (leps[i], nsdl[j], rldl[k])))

for i in range(0,len(lmxi)):
  for j in range(0, len(nsdl)):
    for k in range (0, len(rldl)):
      jobs.append (Job('0b%d_%d_%d' % (i,j,k), 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo RG -lmxi %d -leps 1E-3 -nsdl %g -rldl %s -prfx B' % (lmxi[i], nsdl[j], rldl[k])))

for i in range(0,len(lepsC)):
  for j in range(0, len(nsdl)):
    for k in range (0, len(rldl)):
      jobs.append (Job('0c%d_%d_%d' % (i,j,k), 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo RG -lmxi %d -leps %g -nsdl %g -rldl %s -prfx C' % (lmxiC[i], lepsC[i], nsdl[j], rldl[k])))

for i in range(0,len(leps)):
  for j in range(0, len(nsdl)):
    for k in range (0, len(rldl)):
      jobs.append (Job('1a%d_%d_%d' % (i,j,k), 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo PR -lmxi 1000 -leps %g -nsdl %g -rldl %s -prfx A' % (leps[i], nsdl[j], rldl[k])))

for i in range(0,len(lmxi)):
  for j in range(0, len(nsdl)):
    for k in range (0, len(rldl)):
      jobs.append (Job('1b%d_%d_%d' % (i,j,k), 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo PR -lmxi %d -leps 1E-3 -nsdl %g -rldl %s -prfx B' % (lmxi[i], nsdl[j], rldl[k])))

for i in range(0,len(lepsC)):
  for j in range(0, len(nsdl)):
    for k in range (0, len(rldl)):
      jobs.append (Job('1c%d_%d_%d' % (i,j,k), 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo PR -lmxi %d -leps %g -nsdl %g -rldl %s -prfx C' % (lmxiC[i], lepsC[i], nsdl[j], rldl[k])))

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
