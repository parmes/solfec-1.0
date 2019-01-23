""" Test MESH.get_data() including restoring MESH objects from pickled data

    Raises a MeshError if the tests fail, so completion indicates tests passed """

import pickle

# constants:
PKLFILE = 'msh.pkl'

# --- helper stuff ---

class MeshError(Exception): pass

def compare_meshes(nodes, elements, surfids, mesh):
    """ Compare original data with what's returned by mesh.get_data()
        
        Raises a MeshError if there is a mismatch
    """

    new_nodes, new_els, new_surfids = mesh.get_data()
    
    check_nodes(nodes, new_nodes)
    
    check_elements(elements, new_els)
    
    check_surfids(surfids, new_surfids)

def check_nodes(nodes1, nodes2):
    """ this is easy as the nodes have to be in the same order """
    if nodes1 != nodes2:
        raise MeshError("""nodes don't match:
nodes1:%s
nodes2:%s
            """ % (nodes1, nodes2))
    
def check_elements(els1, els2):
    """ elements can be given in any order so have to sort them """
    if sorted(nest(els1)) != sorted(nest(els2)):
        raise MeshError("""elements don't match:
els1:%s
els2:%s
                """ % (els1, els2))

def check_surfids(surfids1, surfids2):
    """ just checks WHICH surfids (including GID) exist, not the nodes in each face,
        as the choice of GID is arbitrary
    """
    sids1 = sorted(get_sids(surfids1))
    sids2 = sorted(get_sids(surfids2))
    if sids1 != sids2:
        raise MeshError("""surfids don't match:
surfids1 has:%s
surfids2 has:%s
                """ % (sids1, sids2))
    
def get_sids(surfids):
    """ returns a list with all sids i.e. gid, s1, s2 etc """
    if isinstance(surfids, int): # surfids = gid
        return [surfids]
    elif len(surfids) == 1:  # surfids = [gid]
        return surfids
    else:
        sids = []
        sids.append(surfids[0]) # gid
        i = 1
        while i < len(surfids):
            fn = surfids[i]
            nodes = surfids[i+1:i+fn+1]
            sid = surfids[i+fn+1]
            sids.append(sid)
            i += fn + 1 + 1
            
    return sids

def nest(elements):
    """ Convert flat sequence of elements into a a list-of-tuples.
    
        Each tuple represents an element as (n0, n1, ..., nX, vid) """

    i = 0
    zippedels = []
    while i < len(elements):
        num_nodes = elements[i]
        nodes = elements[i+1:i+num_nodes+1]
        vid = elements[i+num_nodes+1]
        zipel = tuple(nodes + [vid])
        zippedels.append(zipel)
        i += num_nodes+1+1
    return zippedels

# --- define the actual tests ---        
def test1():
    """ using mesh from inp/prbpinbar.py  - has no internal elements """

    nodes = [-0.05, -0.05, 0.0,
            0.05, -0.05, 0.0,
            0.05,  0.05, 0.0,
       -0.05,  0.05, 0.0,
       -0.05, -0.05, 1.0,
        0.05, -0.05, 1.0,
        0.05,  0.05, 1.0,
       -0.05,  0.05, 1.0]

    elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, 0]
    surfids = 0
    
    m = MESH(nodes, elements, surfids)
    m2 = MESH(*m.get_data())
    
    compare_meshes(nodes, elements, surfids, m)
    
    b1 = BODY(sol, 'RIGID', m, mat, 'original')
    b2 = BODY(sol, 'FINITE_ELEMENT', TRANSLATE(m2, (2, 2, 0)), mat, 'new')
    
def test2():
    """ using a HPB mesh with surfids, pickling mesh """

    from abaqusreader import AbaqusInput
    model = AbaqusInput("inp/mesh/81fbi.inp")
    part = model.parts["FB"]
    # for a better check with multiple sids:
    #model = AbaqusInput("hpb_draft8.inp")
    #part = model.parts['P204_noAxKey_V0'.upper()]
    m1 = MESH(part.nodes, part.elements, part.surfids)
    
    with open(PKLFILE, 'wb') as pkl:
        pickle.dump(m1.get_data(), pkl, pickle.HIGHEST_PROTOCOL)
    
    with open(PKLFILE, 'r') as pkl:
        n, e, s = pickle.load(pkl)
    
    m2 = MESH(n, e, s)
    
    # check:
    compare_meshes(part.nodes, part.elements, part.surfids, m2)
    
    # make some bodies for a visual check:
    b1 = BODY(sol, 'RIGID', m1, mat, 'm1')
    b2 = BODY(sol, 'RIGID', TRANSLATE(m2, (1, 0, 0)), mat, 'm2')

def test3():
    """ unpickle the test2 mesh """
    with open(PKLFILE, 'r') as pkl:
        n, e, s = pickle.load(pkl)
    m = MESH(n, e, s)
    b = BODY(sol, 'RIGID', m, mat, 'restored')
    
# --- globals ---
sol = SOLFEC('DYNAMIC', 1E-4, 'out')
mat = BULK_MATERIAL(sol)
ns = NEWTON_SOLVER()


# -- run the tests ---
passes = []
tests = (test1, test2, test3)
for test in tests:
    print '\nrunning test:', test.__name__
    try:
        test()
        passes.append('%s : PASS' % test.__name__)
    except MeshError as e:
        errmsg = 'Test %s: %s' % (test.__name__, e.args[0])
        raise MeshError(errmsg)
print "\nTest summary:"
for v in passes:
    print ' ', v
print len(passes), 'out of', len(tests), 'tests passed'
print