""" test overwriting "default" surface material pairs """
nodes = [-1.0, -1.0, 0.0,
          1.0, -1.0, 0.0,
          1.0, 1.0, 0.0,
         -1.0, 1.0, 0.0,
         -1.0, -1.0, 1.0,
          1.0, -1.0, 1.0,
          1.0, 1.0, 1.0,
         -1.0, 1.0, 1.0]

def block(vid, sid):
    mesh = HEX(nodes, 3, 3, 2, vid, [sid, sid, sid, sid, sid, sid]) #(nodes, i, j, k, volid, surfids | dx, dy, dz)    
    return mesh

solfec = SOLFEC('DYNAMIC', 1E-4, 'out')
bmat = BULK_MATERIAL(solfec)

# create 3x FEA blocks, stacked on top of each other (b1=bottom, b3=top):
bottom = block(1, 0) # surfid = 0
middle = block(1, 2) # surfid = 2
TRANSLATE(middle, (0, 0, 1))
top = block(1, 60) # surfid = 60
TRANSLATE(top, (0, 0, 2))
bottom = BODY(solfec, 'FINITE_ELEMENT', bottom, bmat, 'bottom')
middle = BODY(solfec, 'FINITE_ELEMENT', middle, bmat, 'middle')
top = BODY(solfec, 'FINITE_ELEMENT', top, bmat, 'top')

# constrain the bottom block::
FIX_POINT(bottom, (0.5, 0.5, 0.5))
FIX_POINT(bottom, (-0.5, 0.5, 0.5))
FIX_POINT(bottom, (-0.5, -0.5, 0.5))

# Now apply surface materials: *desired* state is that bottom/middle contact has tangential
# friction but middle/top contact does not, regardless of the order of definitions;
# Experiment with setting True or False before and looking at the RT comonentn in the viewer
if True:
    # first pairing with surface 60, then default pairing
    SURFACE_MATERIAL(solfec, 60, model='SIGNORINI_COULOMB', friction=0, restitution=0.0)
    SURFACE_MATERIAL(solfec, model='SIGNORINI_COULOMB', friction=0.5, restitution=0.0)
else:
    # first default pairing, then pairing with surface 60
    SURFACE_MATERIAL(solfec, model='SIGNORINI_COULOMB', friction=0.5, restitution=0.0)
    SURFACE_MATERIAL(solfec, 60, model='SIGNORINI_COULOMB', friction=0, restitution=0.0)
    

GRAVITY(solfec, (0.0,0.0,-9.81))
solver = NEWTON_SOLVER(maxiter=250)
RUN(solfec, solver, 0.01)
