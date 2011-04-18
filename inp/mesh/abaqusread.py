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
  elset = dict ()
  materials = dict ()
  numele = 0
  def __init__(self):
    self.name = 'Part0'
    self.nodes = []
    self.elements = []
    self.elset = dict ()
    self.materials = dict ()
    self.numele = 0

class Assembly_Item:
  name = 'Item0'
  part = 'Part0'
  translate = (0.0, 0.0, 0.0)
  point = (0.0, 0.0, 0.0)
  direction = (0.0, 0.0, 0.0)
  angle = 0.0
  def __init__(self):
    self.name = 'Item0'
    self.npart = 'Part0'
    self.ntranslate = (0.0, 0.0, 0.0)
    self.npoint = (0.0, 0.0, 0.0)
    self.ndirection = (0.0, 0.0, 0.0)
    self.nangle = 0.0

class Element_Set:
  name = 'Elset0'
  assems = []
  def __init__(self):
    self.name = 'Elset0'
    self.assems = []

class Assembly:
  name = 'Assmebly0'
  items = []
  elsets = []
  rigbods = []
  def __init__(self):
    self.name = 'Assmebly0'
    self.items = []
    self.elsets = []
    self.rigbods = []

def ABAQUS_PARSE (path):

  inp = open (path, 'r')
  lin = inp.readline ()
  volid = 1 # spare volume id
  parts = dict () # all parts
  assems = dict () # all assemblies
  elsets = dict () # all element sets

  if lin.find ('*Heading') < 0:
    print 'Invalid Abaqus formar: *Heading not in first line'
    return None

  while lin != "":

    lin = inp.readline ()

# ------------------------------ Part ---------------------------------

    if lin.find ('*Part') >= 0:

      p = Part ()
      lst = lin.split ()
      p.name = lst [1].strip ('name=')

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
	  type = lst [1].strip ('type=')
	  if type == 'C3D8R': type = 8
	  else:
	    print 'Invalid element type:', type
	    return None
          lin = inp.readline ()
	  lst = lin.split ()
	  while lst [0].replace (',', '').isdigit():
	    p.elements.append (type)
	    for i in range (type):
	      p.elements.append (int (lst [i+1].replace (',', '')) - 1) # zero based
	    p.elements.append (volid)
	    p.numele = p.numele + 1
	    lin = inp.readline ()
	    lst = lin.split ()
	elif lin.find ('*Elset') >= 0:
	  # TODO
          lin = inp.readline ()
	elif lin.find ('*Solid Section') >= 0:
	  # TODO
          lin = inp.readline ()
	else: lin = inp.readline ()

      parts [p.name] = p
      volid = volid + 1

# -------------------------- Assembly ---------------------------------

    elif lin.find ('*Assembly') >= 0:

      lst = lin.split ()
      a = Assembly ()
      a.name = lst [1].strip ('name=')

      lin = inp.readline ()
      while lin.find ('*End Assembly') < 0:

	if lin.find ('*Instance') >= 0:
	  i = Assembly_Item ()
	  lst = lin.split ()
	  i.name = lst [1].strip ('name=').replace (',', '')
	  i.part = lst [2].strip ('part=')
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
	      return None

	    lin = inp.readline ()

	  a.items.append (i)

	  # ---- Instance end ----

        #elif lin.find ('*Elset') >= 0:
	  # TODO
	  #lin = inp.readline ()
	  # ---- Elset end ----

        #elif lin.find ('*Rigid Body') >= 0:
	  # TODO
	  #lin = inp.readline ()
	  # ---- Rigid Body end ----

	lin = inp.readline ()

      assems [a.name] = a

# ----------------------------------------------------------------

  inp.close ()

  ret = dict ()
  ret ['PARTS'] = parts
  ret ['ASSEMS'] = assems
  return ret

def ABAQUS_READ (path, step, outpath):

  solfec = SOLFEC ('DYNAMIC', step, outpath)

  bulkmat = BULK_MATERIAL (solfec, model = 'KIRCHHOFF', young = 15E9, poisson = 0.3, density = 1.8E3)

  surfmat = SURFACE_MATERIAL (solfec, model = 'SIGNORINI_COULOMB', friction = 0.5, restitution = 0.25)

  o = ABAQUS_PARSE (path)

  parts = o ['PARTS']
  assems = o ['ASSEMS']

  for key in assems:
    a = assems [key]
    for item in a.items:
      p = parts [item.part]
      shp = MESH (p.nodes, p.elements, [1])
      if l2 (item.translate) > 0.0: TRANSLATE (shp, item.translate)
      if l2 (item.direction) > 0.0: ROTATE (shp, item.point, item.direction, item.angle)
      BODY (solfec, 'RIGID', shp, bulkmat)

  return solfec
