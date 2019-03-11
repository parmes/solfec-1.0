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
	 'tests/BM03/BM03_hexa.py',
	 'tests/BM03/BM03_tetra.py',
	 'tests/BM04/BM04_hexa.py',
	 'tests/BM04/BM04_tetra.py',
	 'tests/BM06/BM06_hexa.py',
	 'tests/BM07/BM07_hexa.py',
	 'tests/BM08/BM08_fine.py',
	 'tests/BM09/BM09_hexa.py',
	 'tests/BM10/BM10_force.py',
	 'tests/BM10/BM10_gravity.py',
	 'tests/BM11/BM11_all.py',
	 'tests/BM12/BM12_hexa.py',
	 'tests/BM13/BM13_coarse.py',
	 'tests/BM13/BM13_fine.py',
	 'tests/BM14/BM14_hexa.py',
	 'tests/BM15/BM15_mixed.py',
	 'tests/BM19/BM19_mixed.py',
	 'tests/BM20/BM20_spring.py',
	 'tests/BM21/BM21_springs.py',
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
