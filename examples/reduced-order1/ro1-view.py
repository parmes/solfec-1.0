import os, sys
dirpath = os.path.dirname(os.path.realpath(__file__))
# juxtapose three simulations in the viewer
percentage = True
print 'Reading/running FEM-TL model ...',
execfile (dirpath + '/ro1-fem-tl.py')
print 'Reading/running FEM-BC model ...',
execfile (dirpath + '/ro1-fem-bc.py')
print 'Reading/running BC-RO model ...',
execfile (dirpath + '/ro1-reduced.py')
