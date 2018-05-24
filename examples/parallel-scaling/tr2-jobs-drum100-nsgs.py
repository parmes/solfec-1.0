from shutil import copy
import fileinput
import subprocess
import os, sys

class Job:
  def __init__ (self, name, command):
    self.name = name
    self.command = command

jobs = []
# Gauss-Seidel jobs
jobs.append (Job('dr0_e_gs', 'solfec examples/parallel-scaling/rotating-drum.py -subd nsgs -outi 0.3 -kifo RG -solv GS'))
jobs.append (Job('dr1_e_gs', 'solfec examples/parallel-scaling/rotating-drum.py -subd nsgs -outi 0.3 -kifo PR -solv GS'))
jobs.append (Job('dr0_s_gs', 'solfec examples/parallel-scaling/rotating-drum.py -subd nsgs -outi 0.3 -kifo RG -solv GS -sphs'))
jobs.append (Job('dr1_s_gs', 'solfec examples/parallel-scaling/rotating-drum.py -subd nsgs -outi 0.3 -kifo PR -solv GS -sphs'))
# best timing Newton absolute delta jobs (updated after running all absdelta jobs and TR2 shortest runtime script)
jobs.append (Job('dr0_s_ns', 'solfec examples/parallel-scaling/rotating-drum.py -subd absdelta -outi 0.3 -kifo RG -lmxi  15 -leps 0.25 -nsdl 0.0001 -sphs'))
jobs.append (Job('dr1_s_ns', 'solfec examples/parallel-scaling/rotating-drum.py -subd absdelta -outi 0.3 -kifo PR -lmxi  20 -leps 0.25 -nsdl 1E-13 -sphs'))
# best timing Newton relative delta jobs (updated after running all reldelta jobs and TR2 shortest runtime script)
jobs.append (Job('dr2_s_ns', 'solfec examples/parallel-scaling/rotating-drum.py -subd reldelta -outi 0.3 -kifo RG -lmxi 10 -leps 0.25 -nsdl 0.01 -rldl avgWii -sphs'))
jobs.append (Job('dr3_s_ns', 'solfec examples/parallel-scaling/rotating-drum.py -subd reldelta -outi 0.3 -kifo PR -lmxi 5 -leps 0.01 -nsdl 0.01 -rldl avgWii -sphs'))

if '--post' in sys.argv: # delete recreated file
  cmd = 'rm -f out/rotating-drum/nsgs/TIMINGS'
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
  cmd = 'cp out/rotating-drum/nsgs/ITERS tr2-dru100-nsgs-iters'
  print cmd
  process = subprocess.Popen(cmd, shell=True)
  process.wait()
  cmd = 'cp out/rotating-drum/nsgs/RUNTIMES tr2-dru100-nsgs-runtimes'
  print cmd
  process = subprocess.Popen(cmd, shell=True)
  process.wait()
  cmd = 'cp out/rotating-drum/nsgs/TIMINGS tr2-dru100-nsgs-timings'
  print cmd
  process = subprocess.Popen(cmd, shell=True)
  process.wait()
