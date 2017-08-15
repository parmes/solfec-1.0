""" A basic graphical time-history plotter for solfec analyses

    Usage:
        solfec scripts/solplot.py path_to_analysis_script
    
    solfec's "-v" switch should NOT be used - start a separate solfec instance if you want
    to use the viewer at the same time
    
    Requirements:
    - matplotlib must be installed
    - solfec/scripts must be in Python's path (e.g. by putting a .pth file pointing to it in Python26\Lib\site-packages)
    - The name of the analysis script must be a legal python name (so no '-', spaces or leading digits)
    - The analysis script must not create any other Tk windows (so using matplotlib directly is not permitted)
    - matplotlib must have a working tkagg backend (this should be the case by default)
    
    Current Features & GUI commands:
    - Select timehistory to extract from database using the "add a data request" button (2nd-to-right button)
    - Control which timehistories are currently shown using the "show/hide curves" button (right-most button)
    - Highlight curves by clicking them

"""

# stdlib imports:
import os, sys, imp

# site imports:
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as pyplot
# debug above:
#print 'matplotlib version:', matplotlib.__version__
#print pyplot.get_backend() 
#print 'matplotlib:', os.path.abspath(matplotlib.__file__)

# solfec imports:
import solfec
from _solplot.plotter import PlotterWindow

# --- functions ---
def injected_import(module_name, inject_dict, path=None):
    """ Import a module and inject names/values into it.
    
        Either 'path' or sys.path if 'path' is not specified is searched for the named module.
        The names/values specified by inject_dict are then inserted before the module is initialised.
        
        The initialised module is both added to sys.modules so it is globally available, and is
        returned so it can be assigned to a local name.
    """
    
    if path is not None:
        path = [path]
    
    newmodule = imp.new_module('newmodule')
    vars(newmodule).update(inject_dict)
    f, newmodule.__file__, options = imp.find_module(module_name, path)
    exec f.read() in vars(newmodule)
    f.close()
    sys.modules[module_name] = newmodule
    newmodule.__name__ = module_name
    return newmodule

# --- main ---

# get the solfec script to import:
script_path = NON_SOLFEC_ARGV()[0]
script_dir, script_name = os.path.split(script_path)
if script_name.endswith('.py'):
    script_name = script_name[:-3]

# import the script module, injecting all the names from solfec into it:
#  (NB: this is pretty deep hackery and is needed because solfec.exe has only added the names it
#  defines into *THIS* module, not the analysis module  we want to import)
scriptmodule = injected_import(script_name, vars(solfec), script_dir)
#print scriptmodule.__doc__
# actually create a plotter
plotter = PlotterWindow(scriptmodule)
pyplot.show()

