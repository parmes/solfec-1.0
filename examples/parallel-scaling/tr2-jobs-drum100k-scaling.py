from shutil import copy
import fileinput
import subprocess
import os, sys

BASE = '-subd scaling -npar 100000 -outi 0.3 -step 1E-4'
NSRG = '-kifo RG -lmxi 10 -leps 0.25 -nsdl 0.01 -rldl avgWii -solv NS'
NSPR = '-kifo PR -lmxi 5 -leps 0.01 -nsdl 0.01 -rldl avgWii -solv NS'

class Job:
  def __init__ (self, name, command, nodes, ranks, taskspernode=24):
    self.modules = ['icc', 'ifort', 'impi']
    self.name = name
    self.command = command
    self.nodes = nodes
    self.ranks = ranks
    self.taskspernode = taskspernode

jobs = []

jobs.append (Job('dr0ns_1', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s %s' % (BASE, NSRG), 1, 24))
jobs.append (Job('dr0ns_2', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s %s' % (BASE, NSRG), 2, 48))
jobs.append (Job('dr0ns_4', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s %s' % (BASE, NSRG), 4, 96))
jobs.append (Job('dr0ns_8', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s %s' % (BASE, NSRG), 8, 192))
jobs.append (Job('dr0ns_16', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s %s' % (BASE, NSRG), 16, 384))

jobs.append (Job('dr0gs_1', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s -kifo RG -solv GS' % BASE, 1, 24))
jobs.append (Job('dr0gs_2', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s -kifo RG -solv GS' % BASE, 2, 48))
jobs.append (Job('dr0gs_4', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s -kifo RG -solv GS' % BASE, 4, 96))
jobs.append (Job('dr0gs_8', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s -kifo RG -solv GS' % BASE, 8, 192))
jobs.append (Job('dr0gs_16', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s -kifo RG -solv GS' % BASE, 16, 384))

jobs.append (Job('dr1ns_1', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s %s' % (BASE, NSPR), 1, 24))
jobs.append (Job('dr1ns_2', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s %s' % (BASE, NSPR), 2, 48))
jobs.append (Job('dr1ns_4', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s %s' % (BASE, NSPR), 4, 96))
jobs.append (Job('dr1ns_8', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s %s' % (BASE, NSPR), 8, 192))
jobs.append (Job('dr1ns_16', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s %s' % (BASE, NSPR), 16, 384))

jobs.append (Job('dr1gs_1', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s -kifo PR -solv GS' % BASE, 1, 24))
jobs.append (Job('dr1gs_2', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s -kifo PR -solv GS' % BASE, 2, 48))
jobs.append (Job('dr1gs_4', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s -kifo PR -solv GS' % BASE, 4, 96))
jobs.append (Job('dr1gs_8', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s -kifo PR -solv GS' % BASE, 8, 192))
jobs.append (Job('dr1gs_16', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py %s -kifo PR -solv GS' % BASE, 16, 384))

if '--post' in sys.argv: # delete recreated file
  cmd = 'rm -f out/rotating-drum/scaling/TIMINGS'
  print cmd
  process = subprocess.Popen(cmd, shell=True)
  process.wait()

for job in jobs: # schedule jobs
  print '***'
  print '*** scheduling: %s' % job.name
  print '***'
  copy ('examples/parallel-scaling/run.sh.parallel', 'run.sh')
  for line in fileinput.input('run.sh', inplace=True):
    if '#SBATCH --ntasks-per-node=' in line:
      print '#SBATCH --ntasks-per-node=%d' % job.taskspernode
    elif '#SBATCH -n' in line:
      print '#SBATCH -n %d' % job.ranks
    elif '#SBATCH -N' in line:
      print '#SBATCH -N %d' % job.nodes
    elif 'module load' in line and len(job.modules) > 0:
      ln = 'module load'
      for mod in job.modules:
	ln += " %s" % mod
      print ln
    elif 'mpirun solfec-mpi' in line:
      print job.command
    else: print line,

  '''
  if '--post' in sys.argv:
    print job.command.replace('mpirun solfec-mpi', 'solfec')
    os.system('read -s -n 1 -p "Press any key to continue..."')
    print
  else:
    print 'sbatch -J %s run.sh' % job.name
    os.system('read -s -n 1 -p "Press any key to continue..."')
    print
  '''

  if '--post' in sys.argv: process = subprocess.Popen(job.command.replace('mpirun solfec-mpi', 'solfec'), shell=True)
  else: process = subprocess.Popen('sbatch -J %s run.sh' % job.name, shell=True)
  process.wait()

if '--post' in sys.argv: # copy stats to renamed files
  cmd = 'cp out/rotating-drum/scaling/ITERS tr2-dru100-scaling-iters'
  print cmd
  process = subprocess.Popen(cmd, shell=True)
  process.wait()
  cmd = 'cp out/rotating-drum/scaling/RUNTIMES tr2-dru100-scaling-runtimes'
  print cmd
  process = subprocess.Popen(cmd, shell=True)
  process.wait()
  cmd = 'cp out/rotating-drum/scaling/TIMINGS tr2-dru100-scaling-timings'
  print cmd
  process = subprocess.Popen(cmd, shell=True)
  process.wait()
