# the main plotter window

from solfec import SOLFEC

# stdlib imports:
import os   
import Tkinter as Tk

# site imports:
import matplotlib.pyplot as pyplot
import matplotlib
from matplotlib.backends.backend_tkagg import ToolTip

# internal imports:
from .requestdialog import RequestDialog
from .curvedialog import CurveDialog

PACKAGEPATH = os.path.dirname(os.path.abspath(__file__)) # relies on this file not being moved!
requestimgpath = os.path.join(PACKAGEPATH, 'requests.gif')
curvesimgpath = os.path.join(PACKAGEPATH, 'curves.gif')
clearrimg = os.path.join(PACKAGEPATH, 'clear.gif')
    
class PlotterWindow(object):
    """ The main plotter window """
    def __init__(self, analysismodule):
        """ analysismodule: module object, the analysis script
        """
        
        self.scriptname = analysismodule.__name__
        
        # get a solfec object:
        self.solfecs = {}   # any solfec objects we find, name:object
        for name, obj in analysismodule.__dict__.items():
            if isinstance(obj, SOLFEC):
                self.solfecs[name] = obj
                if obj.mode=='WRITE':
                    raise TypeError('SOLFEC object "%s" is still in WRITE mode - run me again' % name)
        
        if len(self.solfecs) == 0:
            raise NameError('No solfec object found in %s' % self.scriptname)
            
        self.curves = []    # list of Curve objects        
        
        self.fig = pyplot.figure()
        self.ax = self.fig.add_subplot(111) #, picker=True)
        self.ax.grid(True)
        self.fig.suptitle(self.scriptname + '.py')
        self.fig.canvas.set_window_title('Solplot: %s' % self.scriptname + '.py')
        
        # get a handle to the figure window's toolbar
        toolbar = self.fig.canvas.manager.toolbar
        self.window = self.fig.canvas.manager.window
        
        # add RequestDialog button to toolbar
        self.requestimg = Tk.PhotoImage(master=toolbar, file=requestimgpath)
        b = Tk.Button(master=toolbar, text='requests', image=self.requestimg, command=self.requestclick)
        b.pack(side=Tk.LEFT)
        ToolTip.createToolTip(b, 'Add a data request')
        
        # add CurveDialog button to toolbar
        self.curvesimg =  Tk.PhotoImage(master=toolbar, file=curvesimgpath)
        b = Tk.Button(master=toolbar, text='curves', image=self.curvesimg, command=self.curvesclick)
        b.pack(side=Tk.LEFT)
        ToolTip.createToolTip(b, 'Show/hide curves')
        
        # attach GUI events
        cid = self.fig.canvas.mpl_connect('pick_event', self.onpick)
    
    def onpick(self, event):
        #print event
        if isinstance(event.artist, matplotlib.axes.Axes):
            pass
        if isinstance(event.artist,  matplotlib.lines.Line2D):
            #print event.mouseevent.button
            curve = [c for c in self.curves if c.artist is event.artist]
            if len(curve):
                curve = curve[0]
                curve.highlighted = not curve.highlighted
                self.legend()
                self.fig.canvas.draw()
        
    def requestclick(self):
        newcurves = RequestDialog(self.window, self.solfecs).result
        if newcurves is not None:
            self.curves.extend(newcurves)
            self.plot()
    
    def curvesclick(self):
        if self.curves:
            CurveDialog(self.window, self.curves)
            self.plot()
        else:
            print 'no curves loaded'
    
    def legend(self):
        legend = self.ax.legend(loc='best', prop={'size':'small'})
        if legend is not None: # will be if no curves hence no legend
            legend.draggable(True)
    
    def plot(self):
        
        if self.curves is not None:
            self.ax.cla()
            
            for c in self.curves:
                if c.visible:
                    c.draw(self.ax)
            
            self.legend()
            self.ax.set_xlabel('analysis time')
            self.ax.grid(True)
            
            self.fig.canvas.draw()
        else:
            print 'nothing to plot'         
            