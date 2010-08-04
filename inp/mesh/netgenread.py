# netgen mesh reading
from solfec import MESH

def point (nodes, i):
  return (nodes [3*i], nodes [3*i+1], nodes [3*i+2])

def sub (a, b):
  return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def cros (a, b):
  return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])

def dot (a, b):
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

def NETGEN_READ (path, volid, surfid):

  inp = open (path, 'r')
  nodes = []
  elements = []

  lin = inp.readline ()
  lst = lin.split ()

  if lst [0] != 'mesh3d':
    print 'Invalid netgen formar'
    return None

  while lin != "":

    lin = inp.readline ()
    lst = lin.split ()

    if len(lst) > 0:

      if lst [0] == 'volumeelements':

	lin = inp.readline ()
	lst = lin.split ()
	m = int (lst [0])

	for i in range (0, m):

	  lin = inp.readline ()
	  lst = lin.split ()

	  if lst [1] != '4':
	    print 'Invalid netgen formar'
	    return None

	  elements.append (4)
	  elements.append (int (lst[2]) - 1)
	  elements.append (int (lst[3]) - 1)
	  elements.append (int (lst[4]) - 1)
	  elements.append (int (lst[5]) - 1)
	  elements.append (volid)

      elif lst [0] == 'points':

	lin = inp.readline ()
	lst = lin.split ()
	m = int (lst [0])

	for i in range (0, m):

	  lin = inp.readline ()
	  lst = lin.split ()

	  nodes.append (float (lst [0]))
	  nodes.append (float (lst [1]))
	  nodes.append (float (lst [2]))

  inp.close ()

  # ensure element orientation
  m = len (elements) / 6
  for i in range (0, m):

    n1 = elements [6*i+1]
    n2 = elements [6*i+2]
    n3 = elements [6*i+3]
    n4 = elements [6*i+4]

    a = point (nodes, n1)
    b = point (nodes, n2)
    c = point (nodes, n3)
    d = point (nodes, n4)

    ab = sub (b, a)
    bc = sub (c, b)
    ad = sub (d, a)

    abc = cros (ab, bc)

    if dot (abc, ad) < 0.0:

      elements [6*i+1] = n3
      elements [6*i+2] = n2
      elements [6*i+3] = n1

  return MESH (nodes, elements, [surfid])
