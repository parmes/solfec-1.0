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
  
  def test09(self):
    """ .nodesets """
    solfec = SOLFEC('DYNAMIC', 1E-3, 'out')
    m = AbaqusInput(solfec, 'MODEL07.inp')
    for inst in m.assembly.instances.values():
      coords = [inst.mesh.node(n) for n in range(inst.mesh.nnod)]
      p = inst.part
      # checked the solfec node IDs below were right using the Viewer
      surf001 = [22, 32, 57, 66, 82, 76, 23, 29, 20, 30, 56, 60, 74, 80, 77, 68, 21, 28, 69, 72]
      surf002 = [10, 22, 42, 48, 57, 76, 104, 11, 14, 18, 8, 20, 43, 49, 53, 55, 68, 96, 88, 78, 100, 9, 12, 16, 70, 84, 92]
      surf003 = [42, 48, 62, 67, 82, 76, 104, 107, 43, 40, 44, 46, 50, 61, 80, 77, 99, 91, 83, 96, 88, 78, 105, 106, 41, 45, 79, 81, 89, 90, 97, 98]
      surf4 = [48, 57, 62, 66, 49, 53, 55, 46, 50, 56, 60, 63, 64, 65, 47, 51, 52, 54, 58, 59]
      self.assertEqual(p.nodesets['SURF001'], surf001)
      self.assertEqual(p.nodesets['SURF002'], surf002)
      self.assertEqual(p.nodesets['SURF003'], surf003)
      self.assertEqual(p.nodesets['SURF4'], surf4)
      

# -- MAIN -----------------------------------------------------------

suite = unittest.TestLoader().loadTestsFromTestCase(Test)
unittest.TextTestRunner(verbosity=4).run(suite)