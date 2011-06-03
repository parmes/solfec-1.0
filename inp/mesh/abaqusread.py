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

def SPLIT (str):
  lst = str.split ()
  out = []
  cnt = 0
  for l in lst:
    if l.find ('"') >= 0: # "some names are like that"
      if cnt == 0:
	out.append (l)
	cnt = cnt + 1
      elif cnt == 1:
	l0 = out [len (out)-1]
	out [len (out)-1] = l0 + l
	cnt = 0
    elif cnt == 1:
      l0 = out [len (out)-1]
      out [len (out)-1] = l0 + l
    else: out.append (l)

  return out

def READLINE (inp):
  return inp.readline().lower()

def ABAQUS_PARSE (path):

  inp = open (path, 'r')
  lin = READLINE (inp)
  volid = 1 # spare volume id
  parts = dict () # all parts
  assems = dict () # all assemblies
  elsets = dict () # all element sets
  materials = dict () # all materials

  if lin.find ('*heading') < 0:
    print 'Invalid Abaqus formar: *Heading not in first line'
    raise NameError ('*Heading')

  while lin != "":

    lin = READLINE (inp)

# ------------------------------ Part ---------------------------------

    if lin.find ('*part') >= 0:

      skip = 0
      p = Part ()
      lst = SPLIT (lin)
      p.name = lst [1].replace ('name=', '')
      nmap = dict ()

      lin = READLINE (inp)
      while lin.find ('*end part') < 0:

	if  lin.find ('*node') >= 0:
	  lin = READLINE (inp)
	  lst = SPLIT (lin)
	  while lst [0].replace (',', '').isdigit():
	    p.nodes.append (float (lst [1].replace (',', '')))
	    p.nodes.append (float (lst [2].replace (',', '')))
	    p.nodes.append (float (lst [3].replace (',', '')))
	    lin = READLINE (inp)
	    lst = SPLIT (lin)
	elif lin.find ('*element') >= 0:
	  lst = SPLIT (lin)
	  type = lst [1].replace ('type=', '').replace (',', '')
	  deflen = 0
	  if type == 'c3d8r':
	    type = 8
	    deflen = 8
	  elif type == 'c3d20':
	    print 'ABAQUS_READ WARNING: converting 20 node bricks into 8 node ones!'
	    type = 8
	    deflen = 20
	  else:
	    print 'ABAQUS_READ WARNING: Invalid element type:', type, ' => The part "', p.name, '" will be skipped!'
	    skip = 1
	    break

          lin = READLINE (inp)
	  lst = SPLIT (lin)
	  while lst [0].replace (',', '').isdigit():
	    p.elements.append (type)
	    for i in range (len (lst)-1):
	      if i < type:
		num = int (lst [i+1].replace (',', '')) - 1 # zero based
		p.elements.append (num)
		if num not in nmap:
		  nmap [num] = len (nmap)

	    while i+1 < deflen:
	      lin = READLINE (inp)
	      lst = SPLIT (lin)
	      for j in range (len (lst)):
		i = i+1
		if i < type:
		  num = int (lst [j].replace (',', '')) - 1 # zero based
		  p.elements.append (num)
		  if num not in nmap:
		    nmap [num] = len (nmap)

	    p.elements.append (volid)
	    lin = READLINE (inp)
	    lst = SPLIT (lin)
	elif lin.find ('*elset') >= 0:
	  # TODO => implement element sets
          lin = READLINE (inp)
	elif lin.find ('*solid section') >= 0:
	  # TODO => implement material to element set mapping
	  lst = SPLIT (lin)
	  for l in lst:
	    if l.find ('material=') >= 0:
	      p.material = l.replace ('material=', '').replace (',', '') # XXX => one material per part
          lin = READLINE (inp)
	else: lin = READLINE (inp)

      if not skip:

	rmap = dict ()
	for n in nmap: # reversed node map (new to old)
	  rmap [nmap [n]] = n

	nodes = []
	for n in rmap: # remap nodes (e.g. for 20 node bricks converetd into 8 node ones)
	  i = rmap [n]
	  nodes.append (p.nodes [i*3])
	  nodes.append (p.nodes [i*3+1])
	  nodes.append (p.nodes [i*3+2])
        p.nodes = nodes

        i = 0
	while i < len (p.elements): # remap element node numbers
	  for j in range (p.elements [i]):
	    p.elements [i+1+j] = nmap [p.elements [i+1+j]]
	  i += p.elements [i]+2

	parts [p.name] = p
	volid = volid + 1

# -------------------------- Assembly ---------------------------------

    elif lin.find ('*assembly') >= 0:

      lst = SPLIT (lin)
      a = Assembly ()
      a.name = lst [1].replace ('name=', '')

      lin = READLINE (inp)
      while lin.find ('*end assembly') < 0:

	if lin.find ('*instance') >= 0:
	  i = Instance ()
	  lst = SPLIT (lin)
	  if lst [1].find ('name=') >= 0 and lst [2].find ('part=') >= 0:
	    i.name = lst [1].replace ('name=', '').replace (',', '')
	    i.part = lst [2].replace ('part=', '').replace (',', '')
	  elif lst [2].find ('name=') >= 0 and lst [1].find ('part=') >= 0:
	    i.name = lst [2].replace ('name=', '').replace (',', '')
	    i.part = lst [1].replace ('part=', '').replace (',', '')
	  else:
	    print 'Invalid line in instance definition:', lin
	    raise NameError (lin)

	  if i.part not in parts:
	    print 'ABAQUS_READ WARNING: Skipping assembly part:', i.part, '!'
	    break

	  lin = READLINE (inp)
	  while lin.find ('*end instance') < 0:
	    lst = SPLIT (lin)
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

	    lin = READLINE (inp)

	  a.instances [i.name] = i

	  # ---- Instance end ----

        elif lin.find ('*elset') >= 0:
	  lst = SPLIT (lin)
	  name = 'null'
	  inst = 'null'
	  for l in lst:
	    if l.find ('elset=') >= 0:
              name = l.replace ('elset=', '').replace (',', '')
	    elif l.find ('instance=') >= 0:
              inst = l.replace ('instance=', '').replace (',', '')
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

        elif lin.find ('*rigid body') >= 0:
	  lst = SPLIT (lin)
	  for l in lst:
	    if l.find ('elset=') >= 0:
	      name = l.replace ('elset=', '')
	      if name in a.elsets:
		a.rigbods.append (name)
	      else:
		print 'Invalid element set for a rigid body:', name
		raise NameError (name)

	  # ---- Rigid Body end ----

	lin = READLINE (inp)

      assems [a.name] = a

# ------------------------------ Material ---------------------------------

    elif lin.find ('*material') >= 0:

      lst = SPLIT (lin)
      m = Material ()
      m.name = lst [1].replace ('name=', '')
      lin = READLINE (inp)
      if lin.find ('*density') >= 0:
        lin = READLINE (inp)
	lst = SPLIT (lin)
	m.density = float (lst [0].replace (',', ''))
      else:
	print 'ABAQUS_READ WARNING: Density definition is missing => Assuming 1.0!'
	m.density = 1.0
      lin = READLINE (inp)
      if lin.find ('*elastic') >= 0:
        lin = READLINE (inp)
	lst = SPLIT (lin)
	m.young = float (lst [0].replace (',', ''))
	m.poisson = float (lst [1].replace (',', ''))
      else:
	print 'ABAQUS_READ WARNING: Elastic parameters definition is missing => Assuming: 1.0 and 0.0!'
	m.young = 1.0
	m.poisson = 0.0

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
    used_insts  = dict ()

    # rigid bodies
    for ekey in a.rigbods:
      e = a.elsets [ekey]
      shp = []
      material = 'Mateiral0'
      blabel = 'Body0'
      for ikey in e.instances:
	used_insts [ikey] = 1 # mark as used
	i = a.instances [ikey]
	p = parts [i.part]
	blabel = p.name
	material = p.material # XXX => pick the last one
	s = MESH (p.nodes, p.elements, surfid)
	if l2 (i.translate) > 0.0: TRANSLATE (s, i.translate)
	if l2 (i.direction) > 0.0: ROTATE (s, i.point, i.direction, i.angle)
	shp.append (s)

      bulkmat = BYLABEL (solfec, 'BULK_MATERIAL',  material)
      if bulkmat != None: BODY (solfec, 'RIGID', shp, bulkmat, label = blabel)
      else:
	print 'Invalid bulk material'
	raise NameError (material)

    for ikey in a.instances:
      if ikey not in used_insts:
	material = 'Mateiral0'
	blabel = 'Body0'
	i = a.instances [ikey]
	p = parts [i.part]
	blabel = p.name
	material = p.material # XXX => pick the last one
	s = MESH (p.nodes, p.elements, surfid)
	if l2 (i.translate) > 0.0: TRANSLATE (s, i.translate)
	if l2 (i.direction) > 0.0: ROTATE (s, i.point, i.direction, i.angle)

	bulkmat = BYLABEL (solfec, 'BULK_MATERIAL',  material)
	if bulkmat != None: BODY (solfec, 'FINITE_ELEMENT', s, bulkmat, label = blabel)
	else:
	  print 'Invalid bulk material'
	  raise NameError (material)
