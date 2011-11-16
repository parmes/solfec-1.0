# pinned box example
# just one body and two fixed point constraints
# ==================================================
# run with mpirun -np 2 and you will see divergence
# same for PRB or FEM and either solver; not for rigid
# ------------
# Guess: inconsitency between child and parent body data

a = 0.125
b = 0.125
c = 0.25
step = 0.002
stop = 0.18

nodes = [-a, -b, 0,
          a, -b, 0,
          a,  b, 0,
         -a,  b, 0,
         -a, -b, c,
          a, -b, c,
          a,  b, c,
         -a,  b, c]

msh = HEX (nodes, 2, 2, 2, 0, [0, 1, 2, 3, 4, 5])

sol = SOLFEC ('DYNAMIC', step, 'out/bug1')
bulk = BULK_MATERIAL (sol, model = 'KIRCHHOFF', young = 15E9, poisson = 0.25, density = 1.8E3)
bod = BODY (sol, 'PSEUDO_RIGID', msh, bulk)
bod.scheme = 'DEF_LIM'
FIX_POINT (bod, (-a, -b, 0))
FIX_POINT (bod, (-a, b, 0))

#gs = GAUSS_SEIDEL_SOLVER (1E-5, 1000)
gs = NEWTON_SOLVER ()
GRAVITY (sol, (0, 0, -10))
RUN (sol, gs, stop)
