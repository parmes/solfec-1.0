from shutil import copy
import fileinput
import subprocess
import os, sys

class Job:
  def __init__ (self, name, command, nodes, ranks, taskspernode=24):
    self.modules = ['icc', 'ifort', 'impi/2015.1.133']
    self.name = name
    self.command = command
    self.nodes = nodes
    self.ranks = ranks
    self.taskspernode = taskspernode

jobs = []

jobs.append (Job('dr0ns_1', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS', 1, 24))
jobs.append (Job('dr0ns_2', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS', 2, 48))
jobs.append (Job('dr0ns_4', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS', 4, 96))
jobs.append (Job('dr0ns_8', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS', 8, 192))
jobs.append (Job('dr0ns_16', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS', 16, 384))

jobs.append (Job('dr0gs_1', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS', 1, 24))
jobs.append (Job('dr0gs_2', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS', 2, 48))
jobs.append (Job('dr0gs_4', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS', 4, 96))
jobs.append (Job('dr0gs_8', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS', 8, 192))
jobs.append (Job('dr0gs_16', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo RG -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS', 16, 384))

jobs.append (Job('dr1ns_1', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS', 1, 24))
jobs.append (Job('dr1ns_2', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS', 2, 48))
jobs.append (Job('dr1ns_4', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS', 4, 96))
jobs.append (Job('dr1ns_8', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS', 8, 192))
jobs.append (Job('dr1ns_16', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv NS', 16, 384))

jobs.append (Job('dr1gs_1', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS', 1, 24))
jobs.append (Job('dr1gs_2', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS', 2, 48))
jobs.append (Job('dr1gs_4', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS', 4, 96))
jobs.append (Job('dr1gs_8', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS', 8, 192))
jobs.append (Job('dr1gs_16', 'mpirun solfec-mpi examples/parallel-scaling/rotating-drum.py -npar 100000 -outi 0.3 -kifo PR -lmxi 5 -leps 0.001 -nsdl 0.001 -solv GS', 16, 384))

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
