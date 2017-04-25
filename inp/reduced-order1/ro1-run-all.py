# run simulation from within a function
# so that all variables are local to the function
# and are cleared after the function exits
def runsim(input_file):
  from inspect import currentframe, getframeinfo
  info = getframeinfo(currentframe())
  end = info.filename.rfind('/') + 1
  path = info.filename[0:end]
  percentage = True
  execfile (path + input_file)

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
