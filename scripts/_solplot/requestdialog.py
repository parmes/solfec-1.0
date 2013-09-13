""" provide a dialog for creating results request to HISTORY() """

# stdlib imports:
import tkSimpleDialog
import Tkinter as Tk

# internal imports:
from .bodydialog import BodyDialog, getbodyname
from .curve import Curve

from solfec import TRANSLATE, DURATION, HISTORY

# define the various items HISTORY() can extract
BASIC_PARAMS = 'CX CY CZ DX DY DZ VX VY VZ SX SY SZ SXY SXZ SYZ MISES'.split()
ENERGY_PARAMS = 'KINETIC INTERNAL EXTERNAL CONTACT FRICTION'.split()
CONTACT_PARAMS = 'GAP R U CR CU'.split()
TIMING_PARAMS = 'TIMINT CONUPD CONDET LOCDYN CONSOL PARBAL'.split()
SOLVER_PARAMS = 'STEP CONS BODS DELBODS NEWBODS GISITERS CSCOLORS GSBOT GSMID GSTOP GSINN GSINIT GSRUN GSCOM GSMCOM MERIT NTITERS'.split()

def str_to_vec(str_vector):
    """ converts the string str_vector to a 3-float (x, y, z) """
    
    vec = str_vector.replace(',',' ').split()
    if len(vec) != 3:
        raise ValueError("can't split %s into 3" % vec)
    return tuple(float(v) for v in vec)

def update_optionmenu(widget, var, items):
    """ updates an optionmeu widget with the new items, setting the var (probably a StringVar) to the 1st item """
        
    menu = widget['menu']
    menu.delete(0, Tk.END)
    for v in items:
        menu.add_command(label=v, command=lambda label=v: var.set(label))
    var.set(items[0])
        
class RequestDialog(tkSimpleDialog.Dialog):
    """ Displays a dialog to request data from Solfec
    
        returns (in .results) a list of new Curve objects
    """
    
    def __init__(self, parent, solfecs):
    
        self.solfecs = solfecs    # dict of solfec objects, by variable name
        self.requests = []
        self.descriptions = []
        self.selectmultiplebodies = False
        
        # call the base classes init:
        tkSimpleDialog.Dialog.__init__(self, parent, title='Request results')
        
    
    def body(self, master):
        """ define the gui - overrides base method """
        
        disabledcolor = 'gray'
        
        # request builder:
        # Note this uses pack geometry manager as it gave a more compact result than grid
        
        fr_solfec = Tk.Frame(master)
        fr_solfec.pack(anchor=Tk.W)
        self.solfecvar = Tk.StringVar()
        Tk.Label(fr_solfec, text='Solfec variable:').pack(side=Tk.LEFT)
        Tk.OptionMenu(fr_solfec, self.solfecvar, *self.solfecs.keys()).pack(side=Tk.LEFT)
        self.solfecvar.set(self.solfecs.keys()[0])
        
        fr_reqbuilder = Tk.LabelFrame(master, text='Build request:')
        fr_reqbuilder.pack()
        
        # type selection:
        fr_type = Tk.LabelFrame(fr_reqbuilder, text='Type:')
        fr_type.pack(side=Tk.LEFT)
        self.typevar = Tk.StringVar()
        Tk.Radiobutton(fr_type, text='Basic', variable=self.typevar, value='Basic', command=self.changereqtype).pack(anchor=Tk.W)
        Tk.Radiobutton(fr_type, text='Energy', variable=self.typevar, value='Energy', command=self.changereqtype).pack(anchor=Tk.W)
        Tk.Radiobutton(fr_type, text='Timing', variable=self.typevar, value='Timing', command=self.changereqtype).pack(anchor=Tk.W)
        Tk.Radiobutton(fr_type, text='Solver', variable=self.typevar, value='Solver', command=self.changereqtype).pack(anchor=Tk.W)
        Tk.Radiobutton(fr_type, text='Contact', variable=self.typevar, value='Contact', command=self.changereqtype).pack(anchor=Tk.W)
        self.typevar.set('Basic')
        
        # item selection:
        fr_item = Tk.LabelFrame(fr_reqbuilder, text='Item:')
        fr_item.pack(side=Tk.LEFT, anchor=Tk.N)
        self.itemtypevar = Tk.StringVar()
        self.bodies = [] # will hold tuples (bodyname, body)
        self.rbsolfec = Tk.Radiobutton(fr_item, text='Solfec', variable=self.itemtypevar, value='Solfec', command=self.changeitemtype)
        self.rbsolfec.pack(anchor=Tk.W)
        self.rbbodies = Tk.Radiobutton(fr_item, text='Bodies', variable=self.itemtypevar, value='Bodies', command=self.changeitemtype)
        self.rbbodies.pack(anchor=Tk.W)
        self.bodyvar = Tk.StringVar()
        self.bodybtn = Tk.Button(fr_item, textvar=self.bodyvar, command=self.getbodies)
        self.bodybtn.pack()
        self.itemtypevar.set('Bodies')
        self.changeitemtype()

        # point selection:
        fr_point = Tk.LabelFrame(fr_reqbuilder, text='Point:')
        fr_point.pack(side=Tk.LEFT, anchor=Tk.N)
        self.pointtypevar = Tk.StringVar()
        self.pointvar = Tk.StringVar()
        self.pointvar.set('0,0,0')
        rboffset = Tk.Radiobutton(fr_point, text='offset', variable=self.pointtypevar, value='offset')
        rboffset.pack(anchor=Tk.W)
        rbglobal = Tk.Radiobutton(fr_point, text='global', variable=self.pointtypevar, value='global')
        rbglobal.pack(anchor=Tk.W)
        self.pointtypevar.set('offset')
        epoint = Tk.Entry(fr_point, textvariable=self.pointvar, disabledbackground=disabledcolor)
        epoint.pack()
        self.pointwidgets = (rboffset, rbglobal, epoint)

        # direction selection:
        fr_dir = Tk.LabelFrame(fr_reqbuilder, text='Direction:')
        fr_dir.pack(side=Tk.LEFT, anchor=Tk.N)
        self.dirvar = Tk.StringVar()
        self.edirection = Tk.Entry(fr_dir, textvariable=self.dirvar, disabledbackground=disabledcolor)
        self.edirection.pack()

        # pair selection:
        fr_pair = Tk.LabelFrame(fr_reqbuilder, text='Pair:')
        fr_pair.pack(side=Tk.LEFT, anchor=Tk.N)
        self.pairvar = Tk.StringVar()
        self.epair = Tk.Entry(fr_pair, textvariable=self.pairvar, width=4, disabledbackground=disabledcolor)
        self.epair.pack()

        # string selection:
        fr_string = Tk.LabelFrame(fr_reqbuilder, text='String:')
        fr_string.pack(side=Tk.LEFT, anchor=Tk.N)
        self.strvar = Tk.StringVar()
        self.omstring = Tk.OptionMenu(fr_string, self.strvar, '')
        self.omstring.pack()
        
        # Add button:
        Tk.Button(fr_reqbuilder, text='Add', command=self.add).pack(fill=Tk.Y, expand=True)
        
        # configure request builder appropriately
        self.changereqtype()
        
        # request list:
        Tk.Label(master, text='Current requests:').pack(anchor=Tk.W)
        self.lst_requests = Tk.Listbox(master)
        self.lst_requests.pack(fill=Tk.BOTH, expand=True)
        
        return fr_type # set initial focus      
    
    def changeitemtype(self):
        """ click on item type radiobutton """

        self.bodies = []
        self.bodyvar.set('<none>')
        
        if self.itemtypevar.get() == 'Solfec':
            self.bodybtn.config(state=Tk.DISABLED)
        else:
            self.bodybtn.config(state=Tk.NORMAL)
            
    def getbodies(self):
        """ select which bodies to get results for """
        
        solfec = self.solfecs[self.solfecvar.get()]
        
        self.bodies = BodyDialog(self, solfec, self.selectmultiplebodies).result
        if not hasattr(self.bodies, '__len__') or len(self.bodies) == 0:
            self.bodyvar.set('<none>')
        elif len(self.bodies) == 1:
            self.bodyvar.set('body %s' % getbodyname(self.bodies[0]))
        else:
            self.bodyvar.set('<%i bodies>' % len(self.bodies))
        #print self.bodies
         
    def add(self):
        """ appends the current status of the request builder to .requests and .descriptions """
        
        if self.typevar.get() == 'Basic':
            if self.bodies == []:
                print 'Invalid request: must pick a body'
                return
            
            body = self.bodies[0]
            bodyname = getbodyname(body)
            
            if self.pointvar.get() == '':
                print 'Invalid request: must pick a point'
                return
            
            point = str_to_vec(self.pointvar.get())
            print 
            if self.pointtypevar.get() == 'global':
                pointname = str(point) + '(global)'
            else:
                vec = point
                #print vec
                point = TRANSLATE(body.center, vec)
                pointname = str(vec) + '(offset)'
            
            entity = self.strvar.get()
            
            self.requests.append((body, point, entity))
            self.descriptions.append((bodyname, pointname, entity))
                
        elif self.typevar.get() == 'Energy':
            
            if self.itemtypevar.get() == 'Bodies':
                if self.bodies == []:
                    print 'Invalid request: must pick a body'
                    return
                elif len(self.bodies) == 1:
                    items = self.bodies[0]
                    itemnames = getbodyname(self.bodies[0])
                else:
                    items = self.bodies
                    itemnames = ','.join(getbodyname(b) for b in self.bodies)
            else:
                items = self.solfecs[self.solfecvar.get()]
                itemnames = 'solfec'
            
            kind = self.strvar.get()
            
            self.requests.append((items, kind))
            self.descriptions.append((itemnames, kind))
            
        elif self.typevar.get() == 'Timing' or self.typevar.get() == 'Solver':
            self.requests.append(self.strvar.get())
            self.descriptions.append(self.strvar.get())
        
        elif self.typevar.get() == 'Contact':
            
            if self.itemtypevar.get() == 'Bodies':
                if self.bodies == []:
                    print 'Invalid request: must pick a body'
                    return
                elif len(self.bodies) == 1:
                    items = self.bodies[0]
                    itemnames = getbodyname(bodies)
                else:
                    items = self.bodies
                    itemnames = ','.join(getbodyname(b) for b in self.bodies)
            else:
                items = self.solfecs[self.solfecvar.get()]
                itemnames = 'solfec'
            
            entity = self.strvar.get()

            if self.dirvar.get() == '': # 2-param form
                self.requests.append((items, entity))
                self.descriptions.append((itemnames, entity))

            else:                       # 4-param form
                directionstr = self.dirvar.get()
                if directionstr.lower() == 'none':
                    direction = None
                else:
                    direction = str_to_vec(directionstr)
                
                if self.pairvar.get() == '':
                    print 'Invalid request: must specify a pair if you specify a direction'
                    return
                elif self.pairvar.get().lower() == 'none':
                    pair = None
                else:
                    pair = tuple(int(v) for v in self.pairvar.get().replace(',',' ').split())
                
                self.requests.append((items, direction, pair, entity))
                self.descriptions.append((itemnames, str(direction), str(pair), entity))

        
        self.lst_requests.insert(Tk.END, str(self.descriptions[-1]))
        
        
    def changereqtype(self):
        """ re-configure the GUI when changing request type """
        
        if self.typevar.get() == 'Basic':
            self.itemtypevar.set('Bodies')
            self.selectmultiplebodies = False
            self.changeitemtype()
            self.rbsolfec.config(state=Tk.DISABLED)
            self.rbbodies.config(state=Tk.NORMAL)
            for w in self.pointwidgets: w.config(state=Tk.NORMAL)
            self.edirection.config(state=Tk.DISABLED)
            self.epair.config(state=Tk.DISABLED)
            self.omstring.config(state=Tk.NORMAL)
            update_optionmenu(self.omstring, self.strvar, BASIC_PARAMS)
            
        elif self.typevar.get() == 'Energy':
            self.itemtypevar.set('Solfec')
            self.selectmultiplebodies = True
            self.changeitemtype()
            self.rbsolfec.config(state=Tk.NORMAL)
            self.rbbodies.config(state=Tk.NORMAL)
            for w in self.pointwidgets: w.config(state=Tk.DISABLED)
            self.edirection.config(state=Tk.DISABLED)
            self.epair.config(state=Tk.DISABLED)
            self.omstring.config(state=Tk.NORMAL)
            update_optionmenu(self.omstring, self.strvar, ENERGY_PARAMS)
            
        elif self.typevar.get() == 'Timing':
            self.itemtypevar.set('Solfec')
            self.changeitemtype()
            self.rbsolfec.config(state=Tk.DISABLED)
            self.rbbodies.config(state=Tk.DISABLED)
            for w in self.pointwidgets: w.config(state=Tk.DISABLED)
            self.edirection.config(state=Tk.DISABLED)
            self.epair.config(state=Tk.DISABLED)
            self.omstring.config(state=Tk.NORMAL)
            update_optionmenu(self.omstring, self.strvar, TIMING_PARAMS)
            
        elif self.typevar.get() == 'Solver':
            self.itemtypevar.set('Solfec')
            self.changeitemtype()
            self.rbsolfec.config(state=Tk.DISABLED)
            self.rbbodies.config(state=Tk.DISABLED)
            for w in self.pointwidgets: w.config(state=Tk.DISABLED)
            self.edirection.config(state=Tk.DISABLED)
            self.epair.config(state=Tk.DISABLED)
            self.omstring.config(state=Tk.NORMAL)
            update_optionmenu(self.omstring, self.strvar, SOLVER_PARAMS)
            
        elif self.typevar.get() == 'Contact':
            self.itemtypevar.set('Solfec')
            self.selectmultiplebodies = True
            self.changeitemtype()
            self.rbsolfec.config(state=Tk.NORMAL)
            self.rbbodies.config(state=Tk.NORMAL)
            for w in self.pointwidgets: w.config(state=Tk.DISABLED)
            self.edirection.config(state=Tk.NORMAL)
            self.epair.config(state=Tk.NORMAL)
            self.omstring.config(state=Tk.NORMAL)
            update_optionmenu(self.omstring, self.strvar, CONTACT_PARAMS)
        
    def apply(self):
        """ define what happens on "OK" - overrides base method """

        if len(self.requests) == 0:
            #print 'no requests'
            self.result = None
            return
            
        newcurves = []

        solfec = self.solfecs[self.solfecvar.get()]
        t0, t1 = DURATION(solfec)
        #print 'getting', self.requests
        #print 'for', t0, t1
        results = HISTORY(solfec, self.requests, t0, t1, progress='ON')
        
        times = results[0]
        for i in range(1, len(results)):
            descr = str(self.descriptions[i-1])
            newcurve = Curve(times, results[i], descr, descr)
            
            newcurves.append(newcurve)
            
        self.result = newcurves