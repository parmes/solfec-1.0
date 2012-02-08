""" Abaqus .inp file reader
    
    Access to the contents of the .inp file is through creating an AbaqusInput object
    which has properties containing the deck contents and some equivalent Solfec objects.
    
    Generation of  BODY() objects needs to be done by the user.
    
    See example.py for example usage.
    
"""

__version__ = '0.2dev'

import sys
from math import sqrt
from solfec import *

modulesolfec = None

def printonce(*args, **kwargs):
  """ Like print statement but only produces output from one CPU when run in parallel.
      Non-keyword arguments are printed separated by space (as for print statement)
      Optionally, pass in mode=value to control in which Solfec modes output occurs:
        'WRITE' = only in write (default)
        'READ'  = only in read
        'BOTH'  = both
  """  
  
  mode = kwargs.get('mode', 'WRITE')  # default to write mode if not specified
  if mode not in ['READ', 'WRITE', 'BOTH']:
    raise ValueError('mode must be one of "READ", "WRITE", or "BOTH", not "%s"' % mode)

  if RANK() != 0:
    return
  else:
    if modulesolfec is None or mode=='BOTH' or mode==modulesolfec.mode:
      # fall back to printing if haven't defined solfec object - useful for debugging
      if isinstance(args, str):
        print args
      else:
        try:
          print ' '.join(str(v) for v in args) # simulate print statement
        except TypeError: # not iterable
          print args

def readinpline(inputfile):
  """ Returns a list-of-strings from the next non-comment line of an .inp file, following the Abaqus format.
      At the end of the file, returns the empty list []
      
      Parsing details (see Abaqus Users Manual Section 1.2.1 "Input syntax rules"):
      1. Entire line converted to uppercase
      2. Line is split on commas
      3. Leading/trailing whitespace is removed from each list element, e.g.
         *Part,     name=test => ['*PART', "NAME=TEST"]
      4. Single- or double-quotes (needed if parameter names include a space) are stripped, e.g.
         *Part,  name="my test" => ['*PART', 'NAME=MY TEST']
      5. If the final field of a dataline is empty it is removed e.g.
         10.8, => ['10.8'] rather than ['10.8', '']
  """
  
  ln = '**'
  while ln.startswith('**'):
    ln = inputfile.readline()
  ln = [e.strip().replace("'", "").replace('"', '') for e in ln.upper().split(',')]
  if ln[-1] == '':  # remove trailing empty fields from datalines (results in emtpy list if no data!)
    ln = ln[:-1]
  return ln

def point (nodes, i):
  return (nodes [3*i], nodes [3*i+1], nodes [3*i+2])

def sub (a, b):
  return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def cros (a, b):
  return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])

def dot (a, b):
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

def l2 (a):
  return sqrt (dot (a, a))

class Part:
  def __init__(self):
    self.name = None          # Abaqus name for part
    self.materialname = None  # string giving name of material
    
    self.nodes = []           # Node data in format for 'node' parameter of Solfec MESH() command
    self.elements = []        # Element data in format for 'element' parameter of Solfec MESH() command
    self.surfids = []         # Surface ID data in format for 'surfids' parameter of the Solfec MESH() command

    # private members:
    self._abqnodes = {}       # key=abaqus node ID (1...), value=(x float, y float, z float)
    self._abqels = {}         # key=abaqus element ID (1...), value= list of abaqus node IDs
    self._abqnsets = {}       # Part-level (non-internal) nodeset dictionary: key=set name, value= list of abaqus node IDs
    
    self.nodemap = {}           # key=abaqus node ID, value=Solfec 0-based continous node ID
    self.elmap = {}           # key=abaqus element ID, value=Solfec 0-based continous element ID
    
class Instance:
  def __init__(self):
    self.name = None          # Abaqus name for instance
    self.mesh = None          # Solfec MESH object for instance - located in assembly space
    self.material = None      # Solfec BULK_MATERIAL object for instance
    self.part = None          # Reference to corresponding Part object
    
    # tuples defining position of instance
    self.translate = (0.0, 0.0, 0.0)
    self.point = (0.0, 0.0, 0.0)
    self.direction = (0.0, 0.0, 0.0)
    self.angle = 0.0
    
class Elset:
  def __init__(self):
    self.name = None
    self.instances = []

class Assembly:
  def __init__(self):
    self.name = None
    self.instances = {}       # Dictionary of Instance objects. Key is Instance.name
    self.elsets = {}

class AbaqusInput:
  
  def __init__(self, solfec, path):
    self.path = path            # String giving path to input file
    self.parts = {}             # Dictionary of Part objects. Key is part name/label
    self.assembly = None        # Single Assembly object (can only have one assembly per input deck)
    self.materials = {}         # Dictionary of Solfec BULK_MATERIAL objects. Key is material name/label
    self.solfec = solfec        # Associated solfec object
    
    # also set module-level solfec variable:
    global modulesolfec
    modulesolfec = solfec
    
    gid = 0   # Global surface ID for unspecified surfaces
    volid = 1 # Volume identifier for all elements
    
    self.parse()  # read the input deck into the instance properties
    
    # Can only do the following after entire deck has been read (else would depend on deck order)
    for p in self.parts.values():
      printonce('processing part "%s":' % p.name)
      
      # re-map node and element numbers
      solfecEID = -1
      lastsolfecNID = -1
      for abaqusEID in p._abqels: # will get these in an arbitrary order - is that a meshing problem??
        solfecEID += 1
        p.elmap[abaqusEID] = solfecEID
        abaqusNIDS = p._abqels[abaqusEID]
        #print 'abaqus:', abaqusEID, abaqusNIDS, '->', solfecEID
        numnodes = len(abaqusNIDS)
        numnodes = {4:4, 6:6, 8:8, 10:4, 15:6, 20:8}[numnodes]  # convert to linear      
        
        # get nodes, creating new ones if needed
        nodeids = []
        for abaqusNID in abaqusNIDS[0:numnodes]:  #  first nodes are corner ones, midside are defined last
          #print 'checking for nid', abaqusNID, 
          try:
            solfecNID = p.nodemap[abaqusNID]
            #print 'got', solfecNID
          except KeyError:
            lastsolfecNID += 1
            solfecNID = lastsolfecNID
            #print 'new', solfecNID
            p.nodemap[abaqusNID] = solfecNID
            nodecoords = p._abqnodes[abaqusNID]
            p.nodes.extend(nodecoords)
          nodeids.append(solfecNID) # push actual node numbers - will be for corner nodes only
        
        # push items to .elements list:
        p.elements.append(numnodes) # number of nodes
        p.elements.extend(nodeids)   # will be corner nodes only
        #print nodeids
        p.elements.append(volid)  # push volume identifier
      
      # generate part.surfids
      p.surfids = [gid,]
      for nsetname in p._abqnsets:
        if nsetname.startswith('SURF'):
          sid = int(nsetname.partition('SURF')[2])
          printonce(' generating surface ID', sid, 'from *NSET', nsetname)
          nsetnodes = set(p._abqnsets[nsetname]) # all nodes in the nset
          
          for eid, elnodes in p._abqels.items():  # potential optimisation: break if we've found all the nodes         
            maxn = {4:4, 6:6, 8:8, 10:4, 15:6, 20:8}[len(elnodes)]
            elnodes = set(elnodes[0:maxn])  # get only corner nodes
            matchnodes = elnodes.intersection(nsetnodes)  # match nodes against each element
            if 2 < len(matchnodes) < 5 :  # i.e. if **one** surface (3 or 4 nodes - always linear now) of this element is in the nodeset
              solfecnodeids = [p.nodemap[v] for v in matchnodes] # convert the matching nodes to solfec ids
              p.surfids.append(len(solfecnodeids))
              p.surfids.extend(solfecnodeids)
              p.surfids.append(sid)

    for i in self.assembly.instances.values():
      print 'processing instance %s:' % i.name
      part = i.part
      try:
        i.material = self.materials[part.materialname]
      except KeyError:
        raise KeyError('Material "%s" for part "%s" not found.' % (part.materialname, part.name)) 
      #print 'part.elements:'
      #for ls in zip(*[iter(part.elements)]*10): print ls

      i.mesh = MESH(part.nodes, part.elements, part.surfids)
      if l2(i.translate) > 0.0: TRANSLATE(i.mesh, i.translate)  # translation is applied before rotation in ABAQUS
      if l2(i.direction) > 0.0: ROTATE(i.mesh, i.point, i.direction, i.angle)

  def parse(self):
    
    # parse-wide variables
    nsetwarnings = ['ELSET', 'INSTANCE']  # will throw an error on finding these parameters used for *nset
    # --------
    
    with open(self.path, 'r') as inp:
      
      # require *HEADING on first line:
      lin = readinpline(inp)
      if not lin[0] == '*HEADING':
        raise NameError('*HEADING not found in first line:' + lin)
      lin = readinpline(inp)
      
      while lin != []:
        # ------------------------------ Part ---------------------------------
        
        if lin[0] == '*PART':
          p = Part()
          p.name = lin[1].partition('NAME=')[2]
          printonce('reading *PART "%s":' % p.name)
          lin = readinpline(inp)
          while not lin[0] == '*END PART':
          
            if lin[0] == '*NODE':
              ## ADD ERROR TRAPPING FOR OTHER PARAMETERS ON *NODE HERE
              lin = readinpline(inp)
              while lin[0].isdigit():
                nodenum = int(lin[0])
                coords = (float(v) for v in lin[1:])
                p._abqnodes[nodenum] = coords
                lin = readinpline(inp)
              
            elif lin[0] == '*ELEMENT':  # might be multiple of these per part
              type = lin[1].partition('TYPE=')[2]
              if 'C3D' not in type:   # All 3D solid elements are ok (i.e. C3D, DC3D.., DCC3D.., AC3D..)
                raise NameError('Invalid element type %s' % type)
              numnodes = int(type.partition('3D')[2].rstrip('RHT')) # remove optional bits to just get number of nodes
              if numnodes not in (4, 6, 8, 10, 15, 20):
                raise NameError('Invalid number of elements %i' % numnodes)
              if numnodes > 8:
                printonce(' WARNING: mid-side nodes will be deleted from %s elements' % type)
              lin = readinpline(inp)            
              while lin[0].isdigit():
                elnum = int(lin[0])
                nodes = [int(v) for v in lin[1:]]
                while len(nodes) < numnodes:
                  lin = readinpline(inp)
                  nodes.extend(int(v) for v in lin)
                p._abqels[elnum] = nodes
                lin = readinpline(inp)
                
            elif lin[0] == '*ELSET':
              # TODO => implement element sets
              lin = readinpline(inp)
              
            elif lin[0] == '*NSET':
              setname = lin[1].partition('NSET=')[2]  #*Nset, nset=SURF1
              generate = True if 'GENERATE' in lin[2:] else False
  
              if 'INTERNAL' in lin[2:]: # check this first  
                printonce(' internal *NSET "%s" ignored' % setname)
                lin = readinpline(inp)
              else:
                printonce(' reading *NSET "%s"' % setname)
                for param in lin[2:]:
                  if param in nsetwarnings:
                    raise NameError(" Can't handle parameter %s on non-internal *NSET %s" % (param, setname))
                nodes = []
                lin = readinpline(inp)
                while lin[0].isdigit():
                  if generate:
                    first, last, increment = [int(v) for v in lin]  # NB "last" is different from range()'s stop!
                    nodes.extend(range(first, last + 1 , increment)) # abaqus node numbers
                  else:
                    nodes.extend(int(n) for n in lin)   # abaqus node numbers
                  lin = readinpline(inp)
                p._abqnsets[setname] = nodes
  
            elif lin[0]  == '*SOLID SECTION':
              # material to element set mapping can be done by the user once elset capability added
              for l in lin:
                if l.startswith('MATERIAL='):
                  p.materialname = l.partition('MATERIAL=')[2]  # NB: only use one material per part (the last one found)
              lin = readinpline(inp)
            
            else:
              lin = readinpline(inp)
              
          self.parts[p.name] = p
          
          # -------------------------- Assembly ---------------------------------
        
        elif lin[0] =='*ASSEMBLY':
        
          if self.assembly != None:
            raise ValueError, 'can only have a single *ASSEMBLY definition'
          self.assembly = Assembly()
          self.assembly.name = lin[1].partition('NAME=')[2]
          
          lin = readinpline(inp)
          while not lin[0] == '*END ASSEMBLY':
            
            if lin[0] == '*INSTANCE':
              i = Instance ()
              i.name = lin[1].partition('NAME=')[2].replace(',', '')
              partname = lin[2].partition('PART=')[2]
              i.part = self.parts[partname]
              lin = readinpline(inp)
              while not lin[0] == '*END INSTANCE':
                if len(lin) == 3:
                  i.translate = (float(lin[0]), float(lin[1]), float(lin[2]))
                elif len (lin) == 7:
                  p0 = (float (lin[0]), float (lin[1]), float (lin[2]))
                  p1 = (float (lin[3]), float (lin[4]), float (lin[5]))
                  i.point = p0
                  i.direction = sub(p1, p0)
                  i.angle = float(lin[6])
                else:
                  raise NameError('Invalid line in instance definition:\n' + ','.join(lin))
                  
                lin = readinpline(inp)
                
              self.assembly.instances[i.name] = i
              
              # ---- Instance end ----  ## REALLY? INDENTATION doesn't look right
              
            elif lin[0] == '*ELSET':
              # TODO => implement element sets
              lin = readinpline(inp)
                  
              # ---- Elset end ----
                  
            lin = readinpline(inp)
          printonce('Read', len(self.assembly.instances), 'assembly instances')
            
          # ------------------------------ Material ---------------------------------
        
        elif lin[0] == '*MATERIAL':
          materialname = lin[1].partition('NAME=')[2]
          printonce('reading *MATERIAL "%s"' % materialname)
          lin = readinpline(inp)
          if lin[0] == '*DENSITY':
            lin = readinpline(inp)
            matdensity = float(lin[0])
          else:
            raise NameError('Density definition is missing: %s' % ','.join(lin))
          lin = readinpline(inp)
          if lin[0] == '*ELASTIC':
            lin = readinpline(inp)
            matyoung = float(lin[0])
            matpoisson = float(lin[1])
          else:
            raise NameError('Elastic parameters definition is missing: %s' % ','.join(lin))
          self.materials[materialname] = BULK_MATERIAL(self.solfec, 'KIRCHHOFF', materialname, matyoung, matpoisson, matdensity)
        
        lin = readinpline(inp)    
      # ----------------------------------------------------------------
    

