""" Example of a viewer script

    Use tools/run python script then at the prompt type one of the following, assuming
    that the scripts directory is in the pwd:
    
        scripts/viewscript
            (to print the number of bodies)
                
        scripts/viewscript show [prefix]
            (to render the bodies who's labels start with prefix, or all bodies if no prefix given)
        
        scripts/viewscript show1
            (to render the last body)
    
    Note that this script is executed in the current analysis script context, i.e. as if that script
    had called execfile()
    
    The script name and arguments are passed in sys.argv
"""

import sys
print 'command:', sys.argv
print 'have %i bodies' % len(solfec.bodies)

if len(sys.argv) > 1 and sys.argv[1] == "show":
    # test rendering a list of bodies
    if len(sys.argv) > 2:
        prefix = sys.argv[2]
        selection = [b for b in solfec.bodies if b.label.startswith(prefix)]
    else:
        selection = solfec.bodies
    RENDER(solfec, selection)
    print 'rendered %i bodies in a list' % len(selection)

if len(sys.argv) > 1 and sys.argv[1] == "show1":
    # test rendering a single body
    selection = solfec.bodies[-1]
    print selection
    RENDER(solfec, selection)
    print 'rendered a single body object'
    