# Abaqus model reading (partial)
from math import sqrt
from solfec import *

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
  name = 'Part0'
  nodes = []
  elements = []
  material = 'Material0' # XXX => single material
  def __init__(self):
    self.name = 'Part0'
    self.nodes = []
    self.elements = []
    self.material = 'Material0' # XXX => single material

class Instance:
  name = 'Item0'
  part = 'Part0'
  translate = (0.0, 0.0, 0.0)
  point = (0.0, 0.0, 0.0)
  direction = (0.0, 0.0, 0.0)
  angle = 0.0
  def __init__(self):
    self.name = 'Item0'
    self.part = 'Part0'
    self.translate = (0.0, 0.0, 0.0)
    self.point = (0.0, 0.0, 0.0)
    self.direction = (0.0, 0.0, 0.0)
    self.angle = 0.0

class Elset:
  name = 'Elset0'
  instances = [] # XXX => complete instances only
  def __init__(self):
    self.name = 'Elset0'
    self.instances = []

class Assembly:
  name = 'Assmebly0'
  instances = dict ()
  elsets = dict ()
  rigbods = []
  def __init__(self):
    self.name = 'Assmebly0'
    self.instances = dict ()
    self.elsets = dict ()
    self.rigbods = []

class Material:
  name = 'Material0'
  density = 1.0
  young = 1.0
  poisson = 0.0
  def __init__(self):
    self.name = 'Material0'
    self.density = 1.0
    self.young = 1.0
    self.poisson = 0.0

def ABAQUS_PARSE (path):

  inp = open (path, 'r')
  lin = inp.readline ()
  volid = 1 # spare volume id
  parts = dict () # all parts
  assems = dict () # all assemblies
  elsets = dict () # all element sets
  materials = dict () # all materials

  if lin.find ('*Heading') < 0:
    print 'Invalid Abaqus formar: *Heading not in first line'
    raise NameError ('*Heading')

  while lin != "":

    lin = inp.readline ()

# ------------------------------ Part ---------------------------------

    if lin.find ('*Part') >= 0:

      p = Part ()
      lst = lin.split ()
      p.name = lst [1].lstrip ('name=')

      lin = inp.readline ()
      while lin.find ('*End Part') < 0:

	if  lin.find ('*Node') >= 0:
	  lin = inp.readline ()
	  lst = lin.split ()
	  while lst [0].replace (',', '').isdigit():
	    p.nodes.append (float (lst [1].replace (',', '')))
	    p.nodes.append (float (lst [2].replace (',', '')))
	    p.nodes.append (float (lst [3].replace (',', '')))
	    lin = inp.readline ()
	    lst = lin.split ()
	elif lin.find ('*Element') >= 0:
	  lst = lin.split ()
	  type = lst [1].lstrip ('type=')
	  if type == 'C3D8R': type = 8
	  else:
	    print 'Invalid element type:', type
	    raise NameError (str(type))
          lin = inp.readline ()
	  lst = lin.split ()
	  while lst [0].replace (',', '').isdigit():
	    p.elements.append (type)
	    for i in range (type):
	      p.elements.append (int (lst [i+1].replace (',', '')) - 1) # zero based
	    p.elements.append (volid)
	    lin = inp.readline ()
	    lst = lin.split ()
	elif lin.find ('*Elset') >= 0:
	  # TODO => implement element sets
          lin = inp.readline ()
	elif lin.find ('*Solid Section') >= 0:
	  # TODO => implement material to element set mapping
	  lst = lin.split ()
	  for l in lst:
	    if l.find ('material=') >= 0:
	      p.material = l.lstrip ('material=') # XXX => one material per part
          lin = inp.readline ()
	else: lin = inp.readline ()

      parts [p.name] = p
      volid = volid + 1

# -------------------------- Assembly ---------------------------------

    elif lin.find ('*Assembly') >= 0:

      lst = lin.split ()
      a = Assembly ()
      a.name = lst [1].lstrip ('name=')

      lin = inp.readline ()
      while lin.find ('*End Assembly') < 0:

	if lin.find ('*Instance') >= 0:
	  i = Instance ()
	  lst = lin.split ()
	  i.name = lst [1].lstrip ('name=').replace (',', '')
	  i.part = lst [2].lstrip ('part=')
	  lin = inp.readline ()
	  while lin.find ('*End Instance') < 0:
	    lst = lin.split ()
	    if len (lst) == 3:
	      i.translate = (float (lst [0].replace (',', '')),
			     float (lst [1].replace (',', '')),
			     float (lst [2].replace (',', '')))
            elif len (lst) == 7:
	      p0 = (float (lst [0].replace (',', '')),
		    float (lst [1].replace (',', '')),
		    float (lst [2].replace (',', '')))
	      p1 = (float (lst [3].replace (',', '')),
		    float (lst [4].replace (',', '')),
		    float (lst [5].replace (',', '')))
	      i.point = p0
	      i.direction = sub (p1, p0)
	      i.angle = float (lst [6].replace (',', ''))
	    else:
	      print 'Invalid line in instance definition:', lin
	      raise NameError (lin)

	    lin = inp.readline ()

	  a.instances [i.name] = i

	  # ---- Instance end ----

        elif lin.find ('*Elset') >= 0:
	  lst = lin.split ()
	  name = 'null'
	  inst = 'null'
	  for l in lst:
	    if l.find ('elset=') >= 0:
              name = l.lstrip ('elset=').replace (',', '')
	    elif l.find ('instance=') >= 0:
              inst = l.lstrip ('instance=').replace (',', '')
	  if name != 'null' and inst != 'null':
	    if name in a.elsets:
	      a.elsets [name].instances.append (inst)
	    else:
	      e = Elset ()
	      e.name = name
	      a.elsets [name] = e
	      e.instances.append (inst)

	  # TODO => implement element group selections

	  # ---- Elset end ----

        elif lin.find ('*Rigid Body') >= 0:
	  lst = lin.split ()
	  for l in lst:
	    if l.find ('elset=') >= 0:
	      name = l.lstrip ('elset=')
	      if name in a.elsets:
		a.rigbods.append (name)
	      else:
		print 'Invalid element set for a rigid body:', name
		raise NameError (name)

	  # ---- Rigid Body end ----

	lin = inp.readline ()

      assems [a.name] = a

# ------------------------------ Material ---------------------------------

    elif lin.find ('*Material') >= 0:

      lst = lin.split ()
      m = Material ()
      m.name = lst [1].lstrip ('name=')
      lin = inp.readline ()
      if lin.find ('*Density') >= 0:
        lin = inp.readline ()
	lst = lin.split ()
	m.density = float (lst [0].replace (',', ''))
      else:
	print 'Density definition is missing'
	raise NameError (lin)
      lin = inp.readline ()
      if lin.find ('*Elastic') >= 0:
        lin = inp.readline ()
	lst = lin.split ()
	m.young = float (lst [0].replace (',', ''))
	m.poisson = float (lst [1].replace (',', ''))
      else:
	print 'Elastic parameters definition is missing'
	raise NameError (lin)

      materials [m.name] = m
 
# ----------------------------------------------------------------

  inp.close ()

  ret = dict ()
  ret ['PARTS'] = parts
  ret ['ASSEMS'] = assems
  ret ['MATERIALS'] = materials
  return ret

def ABAQUS_READ (path, solfec):

  o = ABAQUS_PARSE (path)

  surfid = 1 # XXX => not used in ABAQUS

  parts = o ['PARTS']
  assems = o ['ASSEMS']
  materials = o ['MATERIALS']

  for mkey in materials:
    m = materials [mkey]
    BULK_MATERIAL (solfec, model = 'KIRCHHOFF', label = m.name, young = m.young, poisson = m.poisson, density = m.density)

  for akey in assems:
    a = assems [akey]
    for ekey in a.rigbods:
      e = a.elsets [ekey]
      shp = []
      material = 'Mateiral0'
      for ikey in e.instances:
	i = a.instances [ikey]
	p = parts [i.part]
	material = p.material # XXX => pick the last one
	s = MESH (p.nodes, p.elements, surfid)
	if l2 (i.translate) > 0.0: TRANSLATE (s, i.translate)
	if l2 (i.direction) > 0.0: ROTATE (s, i.point, i.direction, i.angle)
	shp.append (s)

      bulkmat = BYLABEL (solfec, 'BULK_MATERIAL',  material)
      if bulkmat != None: BODY (solfec, 'RIGID', shp, bulkmat)
      else:
	print 'Invalid bulk material'
	raise NameError (material)
