""" provide a dialog for selecting bodies """

import tkSimpleDialog
import Tkinter as Tk

def getbodyname(body):
  """ takes a solfec body object and returns the label if it has one, else the id as a str"""
  return body.label if body.label else str(body.id)

class BodyDialog(tkSimpleDialog.Dialog):
    """ Displays a dialog to select one or more bodies
        .results will be a tuple of Solfec BODY objects
    """
    
    def __init__(self, parent, solfec, allow_multiple):
    
        # store the solfec object
        self.solfec = solfec
        self.selectmode = Tk.EXTENDED if allow_multiple else Tk.BROWSE
        
        # call the base classes init:
        title = 'Select one or more bodies' if allow_multiple else 'Select a body'
        tkSimpleDialog.Dialog.__init__(self, parent, title=title)       
    
    def body(self, master):
        """ define the gui - overrides base method """
        
        frame = Tk.Frame(master)
        frame.pack()
        scrollbar = Tk.Scrollbar(frame, orient=Tk.VERTICAL)
        self.lst_bodies = Tk.Listbox(frame, yscrollcommand=scrollbar.set, selectmode=self.selectmode)
        scrollbar.config(command=self.lst_bodies.yview)
        scrollbar.pack(side=Tk.RIGHT, fill=Tk.Y)
        self.lst_bodies.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=True)
        
        self.bodies = self.solfec.bodies
        
        for b in self.bodies:
            self.lst_bodies.insert(Tk.END, getbodyname(b))

    def apply(self):
        """ define what happens on "OK" - overrides base method """
        
        self.result = [self.bodies[int(i)] for i in self.lst_bodies.curselection()]