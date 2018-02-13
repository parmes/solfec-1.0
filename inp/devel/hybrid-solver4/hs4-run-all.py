import os, sys
from subprocess import call
d0 = os.path.dirname(os.path.realpath(__file__))

# run solfec command
def solfec(info, inputfile, args = []):
  print info,
  sys.stdout.flush()
  call (['solfec', '-percent', d0+'/'+inputfile] + args)

# run in 'WRITE' mode (calculate)
print 'Running calculations'
print '===================='
solfec ('Baseline dwell run:', 'hs4-solfec.py')
solfec ('Hybrid dwell run -lwf 1.0:', 'hs4-hybrid.py', ['-lwf', '1.0'])
solfec ('Hybrid dwell run -lwf 0.1:', 'hs4-hybrid.py', ['-lwf', '0.1'])
solfec ('Hybrid dwell run -lwf 0.01:', 'hs4-hybrid.py', ['-lwf', '0.01'])
solfec ('Baseline sweep run:', 'hs4-solfec.py', ['-sweep'])
solfec ('Hybrid sweep run -lwf 1.0:', 'hs4-hybrid.py', ['-lwf', '1.0', '-sweep'])
solfec ('Hybrid sweep run -lwf 0.1:', 'hs4-hybrid.py', ['-lwf', '0.1', '-sweep'])
solfec ('Hybrid sweep run -lwf 0.01:', 'hs4-hybrid.py', ['-lwf', '0.01', '-sweep'])

# run in 'READ' mode (post-process)
print '\nRunning post-processing'
print '======================='
solfec ('Reading baseline dwell;', 'hs4-solfec.py')
solfec ('Reading hybrid dwell at lwf 1.0;', 'hs4-hybrid.py', ['-lwf', '1.0'])
solfec ('Reading hybrid dwell at lwf 0.1;', 'hs4-hybrid.py', ['-lwf', '0.1'])
solfec ('Reading hybrid dwell at lwf 0.01;', 'hs4-hybrid.py', ['-lwf', '0.01'])
solfec ('Reading baseline sweep;', 'hs4-solfec.py', ['-sweep'])
solfec ('Reading hybrid sweep at lwf 1.0;', 'hs4-hybrid.py', ['-lwf', '1.0', '-sweep'])
solfec ('Reading hybrid sweep at lwf 0.1;', 'hs4-hybrid.py', ['-lwf', '0.1', '-sweep'])
solfec ('Reading hybrid sweep at lwf 0.01;', 'hs4-hybrid.py', ['-lwf', '0.01', '-sweep'])

# Plotting
execfile (d0 + '/hs4-postp.py')
