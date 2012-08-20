# Example of useage of abaqusreader
# run with solfec -v -w example.py

# Import the AbaqusInput class from the module:
import sys
sys.path.append('.')	# add working directory to module search path
from abaqusreader import AbaqusInput

# Create a SOLFEC object for the new model - this is needed by the reader so it can generate materials:
solfec = SOLFEC('DYNAMIC', 1E-3, 'out')

# Create a new AbaqusInput object from the .inp deck:
model = AbaqusInput(solfec, 'tests/MODEL04.inp')

# Create a Finite Element body for each Instance in the Assembly:
for inst in model.assembly.instances.values():	# .instances is a dict
  label = inst.name	              # use Abaqus instance name
  mesh = inst.mesh	              # solfec MESH object at the instance position
  bulkmat = inst.material	        # solfec BULK_MATERIAL object
  bdy = BODY(solfec, 'FINITE_ELEMENT', mesh, bulkmat, label)
