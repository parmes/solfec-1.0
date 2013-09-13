class Curve(object):
    """ a time/value series to plot """
    def __init__(self, times, values, label, description):
        self.times = times
        self.values = values
        self.label = label
        self.description = description
        self.fmt = 'default'
        self.visible = True
        self.artist = None
        self.axis = None
        self._highlighted = False

    def draw(self, axis):
        self.axis = axis
        self.artist = axis.plot(self.times, self.values, label=self.label, picker=True)[0] # returns a list of lines
    
    @property
    def highlighted(self):
        return self._highlighted
        
    @highlighted.setter
    def highlighted(self, state):
        if state and self.artist and not self._highlighted:
            self._highlightartist = self.axis.plot(self.times, self.values, label=self.label, alpha=0.5, lw=5, c='y')[0] # returns a list of lines
            self._highlighted = state
        if not state and self.artist and self._highlighted:
            self._highlightartist.remove()
            self._highlighted = state
        