tests = ['inp/tests/pinned-bar.py',
          'inp/tests/math-pendulum.py',
          'inp/tests/double-pendulum.py']

print '------------------------------------------------------------------------------------------'
print 'Solfec serial tests'
print '------------------------------------------------------------------------------------------'

for t in tests:
  print 'Running', t, '...', 
  f = open (t)
  exec f
  f.close ()
  print '------------------------------------------------------------------------------------------'
