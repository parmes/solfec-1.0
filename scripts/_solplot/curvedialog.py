""" provide a dialog for controlling curves """

# stdlib imports:
import tkSimpleDialog
import Tkinter as Tk

class CurveDialog(tkSimpleDialog.Dialog):
    """ Displays a dialog to control display of curves """
    
    def __init__(self, parent, curves):
    
        self.curves = curves
        
        # call the base classes init:
        tkSimpleDialog.Dialog.__init__(self, parent, title='Select Visible Curves')
        
    def body(self, master):
        
        frame = Tk.Frame(master)
        frame.pack()
        scrollbar = Tk.Scrollbar(frame, orient=Tk.VERTICAL)
        self.lst_curves = Tk.Listbox(frame, yscrollcommand=scrollbar.set, selectmode=Tk.EXTENDED)
        scrollbar.config(command=self.lst_curves.yview)
        scrollbar.pack(side=Tk.RIGHT, fill=Tk.Y)
        self.lst_curves.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=True)
        
        for c in self.curves:
            self.lst_curves.insert(Tk.END, c.label)
            if c.visible:
                self.lst_curves.select_set(Tk.END)

    def apply(self):
    
        selected_curves = [int(i) for i in self.lst_curves.curselection()] # curselection() returns a list of strings
        for i, c in enumerate(self.curves):
            if i in selected_curves:
                c.visible = True
            else:
                c.visible = False