#Procedure defining base plate obstacle
def base (bl, bw, bh, base_vid, base_sid):
    base = HULL([-bw, -bl, -bh,
           bw, -bl, -bh,
           bw,  bl, -bh,
          -bw,  bl, -bh,
          -bw, -bl,  0,
           bw, -bl,  0,
           bw,  bl,  0,
          -bw,  bl,  0], base_vid, base_sid)
    TRANSLATE (base, (0.5*l, 1.5*w, 0))
    BODY(solfec,'OBSTACLE', base,bulk)

#Procedure defining individual brick shape.
def single_brick (length, width, height, brick_vid, brick_sid):
    shape = CONVEX([0, 0, 0,
                l, 0, 0,
                l, w, 0,
                0, w, 0,
                0, 0, h,
                l, 0, h,
                l, w, h,
                0, w, h],
               [4, 0, 3, 2, 1, brick_sid,
                4, 1, 2, 6, 5, brick_sid,
                4, 2, 3, 7, 6, brick_sid,
                4, 3, 0, 4, 7, brick_sid,
                4, 0, 1, 5, 4, brick_sid,
                4, 5, 6, 7, 4, brick_sid], brick_vid)
    return shape

#Procedure defining brick as body. Translates shape as necessary before converting to object
def instanceSingleBrick(brick, solfec, material, position, type='PSEUDO_RIGID'):
    brick = TRANSLATE(COPY(brick), position)
    if x%2==0:
        ROTATE (brick, (0.5*l, 1.5*w, h), (0, 0, 1), 90)
    brickBody = BODY(solfec, type, brick, material)
    return brickBody

#Procedure defining a single row of bricks. 'brick' calls 'single_brick' procedure
def single_row(solfec, brick, material, w, h, n):
    row = []
    for i in range(n):
        position = (0,w*i,h)
        brickBody = instanceSingleBrick(brick, solfec, material, position)
        row.append(brickBody)
    return row

#Beginning of main body of code

step=1E-3
scheme = 'DEF_EXP'
solfec=SOLFEC ('DYNAMIC', step, 'out/jenga/' + scheme)
#Brick material defined as wood
bulk = BULK_MATERIAL (solfec,
                      model = 'KIRCHHOFF',
                      young = 11E9,
                      poisson = 0.0,
                      density = 0.8e3)

#Ball material defined as steel
bulk_ball = BULK_MATERIAL (solfec,
                      model = 'KIRCHHOFF',
                      young = 209E9,
                      poisson = 0.3,
                      density = 8e3)

#Brick dimensions
h=1.5
w=2.
l=6.
#No. of bricks in row
n=3
#No. of rows in tower
m=16

#Calling procedure to create base
base (20, 20, 1, 5, 5)

#Calling procedures to create tower
brickShape = single_brick(l, w, h, 1, 1)
for x in range(m):
    tower = single_row (solfec, brickShape, bulk, w, h*x, n)
    for b in tower: b.scheme = scheme

#Procedure creating a ball to impact the tower
sphere = SPHERE((4,10,16), 3, 4, 4)
body = BODY (solfec, 'RIGID', sphere, bulk_ball)
INITIAL_VELOCITY (body, (0, -20, 0), (0, 0, 0))

#Miscellaneous coding required to complete model
SURFACE_MATERIAL (solfec, model='SIGNORINI_COULOMB', friction =0.5, restitution=0.0)
GRAVITY (solfec, (0, 0, -10))
gs= GAUSS_SEIDEL_SOLVER (1E-3 , 1000)
RUN (solfec, gs, 1.0)
