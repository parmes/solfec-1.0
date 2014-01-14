""" Example of a viewer script

    Use tools/run python script then at the prompt type    
        scripts/viewscript
    (assuming pwd is in solfec/)
    
    Note that this script is executed in the current analysis script context, i.e. as if that script
    had called execfile()
    
    The script name and arguments are passed in sys.argv
"""

import sys
print 'command:', sys.argv
print 'have %i bodies' % len(solfec.bodies)
