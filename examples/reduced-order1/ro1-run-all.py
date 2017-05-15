# run simulation from within a function
# so that all variables are local to the function
# and are cleared after the function exits
def runsim(input_file):
  import os, sys
  dirpath = os.path.dirname(os.path.realpath(__file__))
  execfile (dirpath + '/' + input_file)

# %-based output progress
percentage = True

# run in 'WRITE' mode (calculate)
print 'Running calculations'
print '===================='
print 'Running FEM-TL ... ',
runsim ('ro1-fem-tl.py')
print 'Running FEM-BC ... ',
runsim ('ro1-fem-bc.py')
print 'Running POD ... '
runsim ('ro1-modred.py')
print 'Running FEM-BC-RO ... ',
runsim ('ro1-reduced.py')

# run in 'READ' mode (post-process)
print '\nRunning post-processing'
print '======================='
print 'Reading FEM-TL ... '
runsim ('ro1-fem-tl.py')
print 'Reading FEM-BC ... '
runsim ('ro1-fem-bc.py')
print 'Reading FEM-BC-RO ... '
runsim ('ro1-reduced.py')
# Plotting
runsim ('ro1-postp.py')
