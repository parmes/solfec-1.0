# test suite for abaqusreader.py

# Very skeleton at the moment: doens't display anything, doesn't try to raise any exceptions

""" Usage (from within /tests):
    solfec -w tests.py      # -v doesn't work properly using unittest
"""
import bz2
import sys, os
import unittest
sys.path.append('../')  # module to be tested is in level above
from abaqusreader import AbaqusInput

def bunzip(path):
  """ Uncompresses a .bz2 file (leaving the compressed one there) and returns the path to
      the uncompressed file.  You can delete this when finished using os.remove() """

  # NB: not currently used as test .inp file is no longer compressed...

  # inp/mbfcp/array_sin_sweep.py does the same job by calling the external bunzip executable
  # via commands.getoutput().  The bunzip executable is unlikely to be available on Windows and
  # the commands module is Linux specific.  So this uses python's bz2 module instead for portability.

  inf = bz2.BZ2File(path, 'r')
  outfilename = inf.name.rstrip('bz2')
  with open(outfilename, 'w') as outf:
    outf.writelines(inf.readlines())
  inf.close()  # with context-manager not supported by BZ2File before Py2.7
  return outfilename

class Test(unittest.TestCase):

  def test01(self):
    """ Test we can read the demo .inp which comes with Solfec """
    solfec = SOLFEC('DYNAMIC', 1E-3, 'out')
    m = AbaqusInput(solfec, '../../81array.inp')    # in solfec/inp/mesh
    self.failUnless(isinstance(m, AbaqusInput))
  """
  def test03(self):
    solfec = SOLFEC('DYNAMIC', 1E-3, 'out')
    self.failUnless(isinstance(AbaqusInput(solfec, 'MODEL01.inp'), AbaqusInput))

  def test02(self):
    solfec = SOLFEC('DYNAMIC', 1E-3, 'out')
    self.failUnless(isinstance(AbaqusInput(solfec, 'MODEL02.inp'), AbaqusInput))

  def test03(self):
    solfec = SOLFEC('DYNAMIC', 1E-3, 'out')
    self.failUnless(isinstance(AbaqusInput(solfec, 'MODEL03.inp'), AbaqusInput))
  """
  def test04(self):
    """ C3D10 elements """
    solfec = SOLFEC('DYNAMIC', 1E-3, 'out')
    self.failUnless(isinstance(AbaqusInput(solfec, 'MODEL04.inp'), AbaqusInput))

  def test05(self):
    """ C3D6 elements """
    solfec = SOLFEC('DYNAMIC', 1E-3, 'out')
    m = AbaqusInput(solfec, 'MODEL05.inp')
    self.failUnless(isinstance(m, AbaqusInput))

  def test06(self):
    """ C3D6 and C3D10 elements with surface sets """
    solfec = SOLFEC('DYNAMIC', 1E-3, 'out')
    self.failUnless(isinstance(AbaqusInput(solfec, 'MODEL06.inp'), AbaqusInput))

  def test07(self):
    """ C3D20R elements  with surface sets """
    solfec = SOLFEC('DYNAMIC', 1E-3, 'out')
    self.failUnless(isinstance(AbaqusInput(solfec, 'MODEL07.inp'), AbaqusInput))

  def test08(self):
    """ C3D20R elements  with surface sets & 2x instances """
    solfec = SOLFEC('DYNAMIC', 1E-3, 'out')
    self.failUnless(isinstance(AbaqusInput(solfec, 'MODEL08.inp'), AbaqusInput))

# -- MAIN -----------------------------------------------------------

suite = unittest.TestLoader().loadTestsFromTestCase(Test)
unittest.TextTestRunner(verbosity=4).run(suite)