import sys

tests = [
         'tests/pinned-bar.py',
         'tests/math-pendulum.py',
         'tests/double-pendulum.py',
         'tests/projectile.py',
	 'tests/block-sliding.py',
	 'tests/arch.py',
	 'tests/BM01/BM01_pressure.py',
	 'tests/BM01/BM01_force.py',
	 'tests/BM02/BM02_hexa.py',
	 'tests/BM02/BM02_tetra.py',
	 ]

print '------------------------------------------------------------------------------------------'
print 'Solfec serial tests'
print '------------------------------------------------------------------------------------------'

for t in tests:
  print 'Running', t, '...', 
  sys.stdout.flush ()
  f = open (t)
  exec f
  f.close ()
  print '------------------------------------------------------------------------------------------'
