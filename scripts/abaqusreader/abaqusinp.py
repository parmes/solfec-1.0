""" Abaqus .inp file reader
    
    Access to the contents of the .inp file is through creating an AbaqusInput object
    which has properties containing the deck contents and some equivalent Solfec objects.
    
    Generation of BODY() objects needs to be done by the user.
    
    See example.py for example usage.
    
"""

# Developer's notes:
#    The design of the api for this is hampered by maintaining backward compatibility.
#    Ideally, it would separate out the various things which occur, i.e. have parse() which
#    just reads the input deck, and convert() which takes the various solfec parameters needed
#    and actually generates the solfec-objects. It could then have a convenience function which
#    does "everything".
#    
#    Code is being converted to PEP8 compliance as and when it is touched.    

import sys
from math import sqrt
from solfec import *

# Constants:
NSET_ILLEGAL_PARAMS = ['ELSET', 'INSTANCE']  # parameters not allowed for *nset - will raise an err

# Globals:
_solfecmode = None # global to control solfec mode, if available
_verbosity = 0

def printonce(*args, **kwargs):
  """ Like print statement but only produces output from one CPU when run in parallel.
      Non-keyword arguments are printed separated by space (as for print statement)
      
      Optionally, use the following keyword args to control output:
      
      mode (str): Control in which Solfec modes output occurs:
        'WRITE' = only in write (default)
        'READ'  = only in read
        'BOTH'  = both
      
      pri (int, default=0): Output is suppressed if pri <= global _verbosity
                            So by default, no output occurs.
  """  
  
  pri = kwargs.get('pri', 0)
  if pri <= _verbosity:
    return
  
  mode = kwargs.get('mode', 'WRITE')  # default to write mode if not specified
  if mode not in ['READ', 'WRITE', 'BOTH']:
    raise ValueError('"mode" must be one of "READ", "WRITE", or "BOTH", not "%s"' % mode)

  if RANK() != 0:
    return
  else:
    if _solfecmode is None or mode=='BOTH' or mode==_solfecmode:
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

def _sub (a, b):
  return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def _dot (a, b):
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

def _l2 (a):
  return sqrt (_dot (a, a))

class ElasticMaterial(object):
    def __init__(self):
        
        # always available:
        self.name = None    # float
        self.density = None # float
        self.Ym = None      # float
        self.poisson = None # float
        
        # available after .convert() has been called:
        self.bulkmat = None

    def convert(self, solfec):
        self.bulkmat = BULK_MATERIAL(solfec, 'KIRCHHOFF', self.name, self.Ym, self.poisson, self.density)
        
   
    

    
#class Elset(object): # UNUSED at present
#  def __init__(self):
#    self.name = None
#    self.instances = []

class Assembly(object):
    def __init__(self):
        self.name = None    # str, Abaqus name for assembly
        self.instances = {} # Instance objects, key=Instance.name
        self.elsets = {}    # empty at present

class Part(object):
    def __init__(self):
    
        # always available:
        self.name = None            # str, Abaqus name for part
        self.materialname = None    # str, Abaqus name for material
        self.abqnodes = {}          # key=abaqus node ID (1...), value=(x float, y float, z float)
        self.abqelements = {}       # key=abaqus element ID (1...), value=(typestr, abaqus node IDs, ...)
        self.abqnodesets = {}       # Part-level (non-internal) nodeset dictionary: key=set name, value= list of abaqus node IDs
        
        # available after .convert() has been called:    
        self.nodes = []     # Node data in format for 'node' parameter of Solfec MESH()
        self.elements = []  # Element data in format for 'element' parameter of Solfec MESH()
        self.surfids = []   # Surface ID data in format for 'surfids' parameter of the Solfec MESH()
        self.nodesets = {}  # Node IDs (solfec node numbers) for *NODESET
        self.nodemap = {}   # key=abaqus node ID, value=Solfec 0-based continous node ID
        self.elmap = {}     # key=abaqus element ID, value=Solfec 0-based continous element ID

    def convert(self, volid, gid):
      printonce('processing part "%s":' % self.name)
      
      # create self.elements:
      # re-map node and element numbers
      solfecEID = -1
      lastsolfecNID = -1
      for abaqusEID in self.abqelements: # will get these in an arbitrary order - is that a meshing problem??
        solfecEID += 1
        self.elmap[abaqusEID] = solfecEID
        typestr = self.abqelements[abaqusEID][0]
        if 'C3D' not in typestr:   # All 3D solid elements are ok (i.e. C3D, DC3D.., DCC3D.., AC3D..)
            raise ValueError('Cannot convert "%s" elements for part "%s", must be 3D solids' % (typestr, self.name))
        abaqusNIDS = self.abqelements[abaqusEID][1:]
        numnodes = len(abaqusNIDS)
        if numnodes > 8:
            printonce(' WARNING: mid-side nodes will be deleted from "%s" elements' % typestr)  
        try:
            numnodes = {4:4, 6:6, 8:8, 10:4, 15:6, 20:8}[numnodes]  # convert to linear      
        except KeyError:
            raise ValueError('Cannot convert "%s" elements for part "%s", must be 4, 6, 8, 10, 15 or 20-noded' % (typestr, self.name))
        
        # get nodes, creating new ones if needed
        nodeids = []
        for abaqusNID in abaqusNIDS[0:numnodes]:  #  first nodes are corner ones, midside are defined last
          #print 'checking for nid', abaqusNID, 
          try:
            solfecNID = self.nodemap[abaqusNID]
            #print 'got', solfecNID
          except KeyError:
            lastsolfecNID += 1
            solfecNID = lastsolfecNID
            #print 'new', solfecNID
            self.nodemap[abaqusNID] = solfecNID
            nodecoords = self.abqnodes[abaqusNID]
            self.nodes.extend(nodecoords)
          nodeids.append(solfecNID) # push actual node numbers - will be for corner nodes only
        
        # push items onto .elements
        self.elements.append(numnodes) # number of nodes
        self.elements.extend(nodeids)   # will be corner nodes only
        self.elements.append(volid)  # push volume identifier
      
      # create self.nodesets:
      for nsetname in self.abqnodesets:
        solfecnodes = [self.nodemap[nid] for nid in self.abqnodesets[nsetname] if self.nodemap.get(nid, None) is not None]
        self.nodesets[nsetname] = solfecnodes
      
      # create self.surfids:
      self.surfids = [gid,]
      for nsetname in self.abqnodesets:
        if nsetname.startswith('SURF'):
          sid = int(nsetname.partition('SURF')[2])
          printonce(' generating surface ID', sid, 'from *NSET', nsetname)
          nsetnodes = set(self.abqnodesets[nsetname]) # all nodes in the nset
          
          for eid in self.abqelements:  # potential optimisation: break if we've found all the nodes
            elnodes = self.abqelements[eid][1:]
            maxn = {4:4, 6:6, 8:8, 10:4, 15:6, 20:8}[len(elnodes)]
            elnodes = set(elnodes[0:maxn])  # get only corner nodes
            matchnodes = elnodes.intersection(nsetnodes)  # match nodes against each element
            if 2 < len(matchnodes) < 5 :  # i.e. if **one** surface (3 or 4 nodes - always linear now) of this element is in the nodeset
              solfecnodeids = [self.nodemap[v] for v in matchnodes] # convert the matching nodes to solfec ids
              self.surfids.append(len(solfecnodeids))
              self.surfids.extend(solfecnodeids)
              self.surfids.append(sid)   

class Instance(object):
    def __init__(self):

        # always available:
        self.name = None          # Abaqus name for instance
        self.part = None          # Reference to corresponding Part object
        # tuples defining position of instance ...
        self.translate = (0.0, 0.0, 0.0)
        self.point = (0.0, 0.0, 0.0)
        self.direction = (0.0, 0.0, 0.0)
        self.angle = 0.0
    
        # available after .convert() has been called:
        self.mesh = None          # Solfec MESH object for instance - located in assembly space
        self.material = None      # Solfec BULK_MATERIAL object for instance - name kept for backward compatability
        
    def convert(self, bulkmat=None):
        
        part = self.part # must already be converted
        
        mesh = MESH(part.nodes, part.elements, part.surfids)
        if _l2(self.translate) > 0.0:
            TRANSLATE(mesh, self.translate)  # translation is applied before rotation in ABAQUS
        if _l2(self.direction) > 0.0:
            ROTATE(mesh, self.point, self.direction, self.angle)
        self.mesh = mesh
        
        self.material = bulkmat
        
class AbaqusInput(object):
  
  def __init__(self, path, solfec=None, gid=0, volid=1, verbose=0, convert=True):
    """ Access information from an Abaqus .inp file
    
        parameters:
            path (str): path to .inp file

        optional parameters:
            solfec: Solfec object, only needed if material properties are required.
            gid (int): Global surface ID for unspecified surfaces (default 0)
            volid (int): Volume identifier for all elements (default 1)
            verbose (int): Higher numbers give increasing detail of processing on stdout
            convert (bool): Create solfec objects (default=True)
    """
    
    print 'using local ver'
    
    self.ver = 'NEW'
    global _verbosity
    _verbosity = verbose
    self.path = path            # str, path to input file
    self.parts = {}             # key (str):part name/label, value:Part object
    self.assembly = None        # Single Assembly object (can only have one assembly per input deck)
    self.materials = {}         # key (str):material name/label, value:ElasticMaterial object
    self.solfec = solfec        # Associated solfec object
    
    # set the solfec mode, if we have a solfec object
    if solfec:
        global _solfecmode
        _solfecmode = solfec.mode
    
    self.gid = gid
    self.volid = volid
    
    self.parse()  # read the input deck into the instance properties
    if convert:
        self.convert()
    
  def convert(self):
    """ Convert instances, and their referenced parts and possibly materials to equivalent
        solfec objects.
    
        This can only be done after the entire deck has been read, else it would require a specific
        deck order,
    """
    for i in self.assembly.instances.values():
        printonce('processing instance %s:' % i.name)
        
        # convert its part:
        i.part.convert(self.volid, self.gid)
        
        # convert its material
        bulkmat = None
        if self.solfec:
            
            matname = i.part.materialname
            try:
                mat = self.materials[matname]
            except KeyError:
                raise KeyError('Material "%s" for instance of part "%s" not found.' % (matname, i.part.name)) 
            bulkmat = mat.convert(self.solfec).bulkmat
        
        # convert the instance itself
        i.convert(bulkmat)
                
  def __str__(self):
    return 'AbaqusObject from "%s" with %s parts, %s instances, %s materials' % (self.path, len(self.parts),
            len(self.assembly.instances), len(self.materials))

  def parse(self):
    
    """ fills out self.parts, self.assembly, self.materials """
        
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
                coords = tuple(float(v) for v in lin[1:])
                p.abqnodes[nodenum] = coords
                lin = readinpline(inp)
              
            elif lin[0] == '*ELEMENT':  # might be multiple of these per part
              
              typestr = lin[1].partition('TYPE=')[2]
              numnodes = int(typestr.partition('3D')[2].rstrip('RHT')) # remove optional bits to just get number of nodes
              lin = readinpline(inp)
              while lin[0].isdigit():
                elnum = int(lin[0])
                nodes = [int(v) for v in lin[1:]]
                while len(nodes) < numnodes:
                  lin = readinpline(inp)
                  nodes.extend(int(v) for v in lin)
                p.abqelements[elnum] = tuple([typestr,] + nodes)
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
                  if param in NSET_ILLEGAL_PARAMS:
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
                p.abqnodesets[setname] = nodes
  
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
                  i.direction = _sub(p1, p0)
                  i.angle = float(lin[6])
                else:
                  raise NameError('Invalid line in instance definition:\n' + ','.join(lin))
                  
                lin = readinpline(inp)
                
              self.assembly.instances[i.name] = i
              
            elif lin[0] == '*ELSET':
              # TODO => implement element sets
              lin = readinpline(inp)
                  
              # ---- Elset end ----
                  
            lin = readinpline(inp)
          printonce('Read', len(self.assembly.instances), 'assembly instances', pri=1)
            
          # ------------------------------ Material ---------------------------------
        
        elif lin[0] == '*MATERIAL':
          mat = ElasticMaterial()
          mat.name = lin[1].partition('NAME=')[2]
          printonce('reading *MATERIAL "%s"' % mat.name)
          lin = readinpline(inp)
          if lin[0] == '*DENSITY':
            lin = readinpline(inp)
            mat.density = float(lin[0])
          else:
            raise NameError('Density definition is missing: %s' % ','.join(lin))
          lin = readinpline(inp)
          if lin[0] == '*ELASTIC':
            lin = readinpline(inp)
            mat.Ym = float(lin[0])
            mat.poisson = float(lin[1])
          else:
            raise NameError('Elastic parameters definition is missing: %s' % ','.join(lin))
          self.materials[mat.name] = mat
        
        lin = readinpline(inp)    
      # ----------------------------------------------------------------
    

