# run simulation from within a function
# so that all variables are local to the function
# and are cleared after the function exits
def runsim(input_file):
  import os
  dirpath = os.path.dirname(os.path.realpath(__file__))
  execfile (dirpath + '/' + input_file)

# run simulations
runsim ('ro0-elong.py')
runsim ('ro0-energy.py')
runsim ('ro0-convtest.py')

