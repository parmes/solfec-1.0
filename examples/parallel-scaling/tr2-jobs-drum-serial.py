from shutil import copy
import fileinput
import subprocess
import os

class Job:
  def __init__ (self, name, command):
    self.name = name
    self.command = command

leps = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 0.25]
lmxi = [3, 5, 10, 15, 20, 25, 30, 40, 50, 70, 100]
nsdl = [0.0, 1E-14, 1E-13, 1E-12, 1E-11, 1E-10, 1E-9, 1E-8, 1E-7, 1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 0.1]
jobs = []

for i in range(0,len(leps)):
  for j in range(0, len(nsdl)):
    jobs.append (Job('dr_a%d_%d' % (i,j), 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo RG -lmxi 1000 -leps %g -nsdl %g -prfx A' % (leps[i], nsdl[j])))

for i in range(0,len(lmxi)):
  for j in range(0, len(nsdl)):
    jobs.append (Job('dr_b%d_%d' % (i,j), 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo RG -lmxi %d -leps 1E-3 -nsdl %g -prfx B' % (lmxi[i], nsdl[j])))

lepsC = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01]
lmxiC = [  5,    5,   5,    5,     5,    5,  10,   10,  10,   10,    10,   10,  15,   15,  15,   15,    15,   15,  20,   20,  20,   20,    20,   20]

for i in range(0,len(lepsC)):
  for j in range(0, len(nsdl)):
    jobs.append (Job('dr_c%d_%d' % (i,j), 'solfec examples/parallel-scaling/rotating-drum.py -outi 0.3 -kifo RG -lmxi %d -leps %g -nsdl %g -prfx C' % (lmxiC[i], lepsC[i], nsdl[j])))

for job in jobs: # schedule jobs
  print '***'
  print '*** scheduling: %s' % job.name
  print '***'
  copy ('examples/parallel-scaling/run.sh.serial', 'run.sh')
  for line in fileinput.input('run.sh', inplace=True):
    if 'solfec input-file.py' in line:
      print job.command
    else: print line,

  #print 'sbatch -J %s run.sh' % job.name
  #os.system('read -s -n 1 -p "Press any key to continue..."')
  #print
  process = subprocess.Popen('sbatch -J %s run.sh' % job.name, shell=True)
  process.wait()
