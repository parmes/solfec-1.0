import sys

tests = ['inp/tests/pinned-bar.py',
         'inp/tests/math-pendulum.py',
         'inp/tests/double-pendulum.py',
         'inp/tests/projectile.py',
	 'inp/tests/block-sliding.py',
	 'inp/tests/arch.py']

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
