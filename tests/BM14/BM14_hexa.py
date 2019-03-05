# CONTACT OF 2 SIMPLE PLATES UNDER PRESSURE

def create_simulation (step, stop, outstep):
    Dir_out='out/tests/BM14_hexa_'+str(step)+'s'
    GEOMETRIC_EPSILON (1E-5)
  
    solfec=SOLFEC ('QUASI_STATIC', step, Dir_out)
    
    nodes = [0.00000000000000E+00,  -2.50000000000000E+00,   5.00000000000000E-01,
             1.00000000000000E+01,  -2.50000000000000E+00,   5.00000000000000E-01,
             1.00000000000000E+01,   2.50000000000000E+00,   5.00000000000000E-01,
             0.00000000000000E+01,   2.50000000000000E+00,   5.00000000000000E-01,
             0.00000000000000E+00,  -2.50000000000000E+00,   0.00000000000000E-01,
             1.00000000000000E+01,  -2.50000000000000E+00,   0.00000000000000E-01,
             1.00000000000000E+01,   2.50000000000000E+00,   0.00000000000000E-01,
             0.00000000000000E+01,   2.50000000000000E+00,   0.00000000000000E-01,]

    Msh_M1 = HEX (nodes, 16, 8, 1, 0, [1, 2, 3, 4, 5,6])
                                       
    bulk = BULK_MATERIAL (solfec,
                          model = 'KIRCHHOFF',
                          young = 2.0E+11,
                          poisson = 0.3,
                          density = 7.850E+03)

    Msh_M2=COPY(Msh_M1)
    TRANSLATE(Msh_M2,(0,0,-0.501))

    Mesh_final = Msh_M1 
    Mesh_final_2 = Msh_M2

    FIX_YZ=[]
    for i in Mesh_final.nodes_on_surface(3):
      FIX_YZ.append(Mesh_final.node(i))

    for i in Mesh_final.nodes_on_surface(5):
      FIX_YZ.append(Mesh_final.node(i))


    FIX_XZ=[]
    for i in Mesh_final.nodes_on_surface(2):
      FIX_XZ.append(Mesh_final.node(i)) 

    for i in Mesh_final.nodes_on_surface(4):
      FIX_XZ.append(Mesh_final.node(i)) 


    FIX_YZ_2=[]
    for i in Mesh_final_2.nodes_on_surface(3):
      FIX_YZ_2.append(Mesh_final_2.node(i))

    for i in Mesh_final_2.nodes_on_surface(5):
      FIX_YZ_2.append(Mesh_final_2.node(i))

    FIX_XZ_2=[]
    for i in Mesh_final_2.nodes_on_surface(2):
      FIX_XZ_2.append(Mesh_final_2.node(i)) 

    for i in Mesh_final_2.nodes_on_surface(4):
      FIX_XZ_2.append(Mesh_final_2.node(i))   
      
        
    N1=(5.0, 0.0, 0.5)
    N2=(5.0, 0.0, 0.0)
    N3=(5.0, 0.0, -0.001)
    N4=(5.0, 0.0, -0.501)
    
    Point_A=N1    
    No_Trait=[]
    No_Trait.append(N1)
    No_Trait.append(N2)

    No_Trait_2=[]
    No_Trait_2.append(N3)
    No_Trait_2.append(N4)
        
    SURFACE_MATERIAL (solfec, 
                      model = 'SIGNORINI_COULOMB', 
                      friction = 0., 
                      restitution = 0.0)
        
    Bod_M1 = BODY (solfec, 'FINITE_ELEMENT', Mesh_final, bulk, label = 'bod0')
    Bod_M2 = BODY (solfec, 'FINITE_ELEMENT', Mesh_final_2, bulk, label = 'bod1')

    Fix_direc_X=(1.,0.,0.)
    Fix_direc_Y=(0.,1.,0.)
    Fix_direc_Z=(0.,0.,1.)    

    for fix in FIX_YZ:
       FIX_DIRECTION (Bod_M1, fix, Fix_direc_Y)
       FIX_DIRECTION (Bod_M1, fix, Fix_direc_Z)
       
    for fix in FIX_XZ:
       FIX_DIRECTION (Bod_M1, fix, Fix_direc_X)
       FIX_DIRECTION (Bod_M1, fix, Fix_direc_Z)

    for fix in No_Trait:
       FIX_DIRECTION (Bod_M1, fix, Fix_direc_X)
       FIX_DIRECTION (Bod_M1, fix, Fix_direc_Y)

    for fix in FIX_YZ_2:
       FIX_DIRECTION (Bod_M2, fix, Fix_direc_Y)
       FIX_DIRECTION (Bod_M2, fix, Fix_direc_Z)
       
    for fix in FIX_XZ_2:
       FIX_DIRECTION (Bod_M2, fix, Fix_direc_X)
       FIX_DIRECTION (Bod_M2, fix, Fix_direc_Z)

    for fix in No_Trait_2:
       FIX_DIRECTION (Bod_M2, fix, Fix_direc_X)
       FIX_DIRECTION (Bod_M2, fix, Fix_direc_Y)
       
    Pressure = -2.5E+08
    PRESSURE (Bod_M1, 6, Pressure)  

    nt = NEWTON_SOLVER (1E-10, locdyn = 'OFF')
    OUTPUT(solfec, outstep)
    if not VIEWER() and solfec.mode == 'WRITE':
      solfec.verbose = '%'
    RUN (solfec, nt, stop)
    return (solfec, N1, N2)

# data
step = 2E-4 # any of: 1E-4, 2E-4, 5E-5
stop = 0.1
outstep = 0.05

# calculate
(solfec, N1, N2) = create_simulation (step, stop, outstep)

# read results
if not VIEWER():
  (solfec, N1, N2) = create_simulation (step, stop, outstep)
  if solfec.mode != 'READ':
    print 'READ ERROR'
    import sys
    sys.exit(1)

  body = BYLABEL(solfec, 'BODY', 'bod0')

  data_dz = {
#       time step
        5E-05 : [
#       point        DX, DY, DZ                  DX, DY, DZ
#       coords     Solfec 841200b                Code_Aster 
        (N1,  (0.0,  0.0,  -5.9060E-02),    (0.0, 0.0, -5.87E-02)),
	(N2,  (0.0,  0.0,  -5.8663E-02),    (0.0, 0.0, -5.78E-02))],
	2E-04 : [
        (N1,  (0.0,  0.0, -5.9840E-02),    (0.0, 0.0, -5.87E-02)),
	(N2,  (0.0,  0.0, -5.9401E-02),    (0.0, 0.0, -5.78E-02))],
	1E-04 : [
        (N1,  (0.0,  0.0, -5.8391E-02),    (0.0, 0.0, -5.87E-02)),
	(N2,  (0.0,  0.0, -5.8017E-02),    (0.0, 0.0, -5.78E-02))] }

# tolerance Solfec c060f54
  tole_dz = {
        5E-05 : [
        (N1,  (1E-6,  1E-6,  0.013)),
	(N2,  (1E-6,  1E-6,  0.0052))],
	2E-04 : [
        (N1,  (1E-6,  1E-6,  0.021)),
	(N2,  (1E-6,  1E-6,  0.028))],
	1E-04 : [
        (N1,  (1E-6,  1E-6,  0.0079)),
	(N2,  (1E-6,  1E-6,  0.016))] }

  def tolerance_epsilon (tolerance):
    if tolerance > 1.0: return 0.01*tolerance
    tstr = "%.15f" % tolerance
    teps = 0.000000000000001
    j = len(tstr)-1
    while j > 0 and tstr[j] == '0':
      teps = teps * 10.
      j = j - 1
    return 0.0 if j == 1 else teps

  failed = False
  SEEK (solfec, stop)
  label = ('DX', 'DY', 'DZ')
  for (item, tole) in zip(data_dz[step], tole_dz[step]):
    pnt = item[0]
    disp = DISPLACEMENT (body, pnt)
    for i in range(0,3):
      Value = disp[i]
      Code_Aster_Ref = item[2][i]
      div =  1.0 if abs(Code_Aster_Ref) == 0.0 else abs(Code_Aster_Ref) 
      error = abs(Value-Code_Aster_Ref)/div 
      tolerance = tole[1][i] + tolerance_epsilon(tole[1][i])

      if error > tolerance:
	failed = True
	print '\b\b\b\bFAILED', 
	print '(Computed displacement %s at point (%g, %g, %g) was %g' % (label[i], pnt[0], pnt[1], pnt[2], Value), 'while reference is %g)' % Code_Aster_Ref,
	print '-->(error %g, while tolerance %g)' % (error, tolerance)

  data_p = {
#       time step
        5E-05 : [
#       point       SX, SY, SZ                                SX, SY, SZ
#       coords     Solfec 841200b                             Code_Aster 
        (N1,  ( 1.04E+09,  2.15E+09,  8.09E+08),    ( 1.01E+09,  2.06E+09,  7.82E+08)),
	(N2,  (-1.14E+09, -2.15E+09, -1.14E+09),    (-1.16E+09, -2.20E+09, -1.15E+09))],
	2E-04 : [
        (N1,  ( 1.04E+09,  2.17E+09,  7.93E+08),    (1.01E+09,   2.06E+09,  7.82E+08)),
	(N2,  (-1.16E+09, -2.15E+09, -1.16E+09),    (-1.16E+09, -2.20E+09, -1.15E+09))],
	1E-04 : [
        (N1,  ( 1.03E+09,  2.13E+09,  8.10E+08),    ( 1.01E+09,  2.06E+09,  7.82E+08)),
	(N2,  (-1.14E+09, -2.14E+09, -1.13E+09),    (-1.16E+09, -2.20E+09, -1.15E+09))] }

# tolerance Solfec e58e6a4
  tole_p = {
        5E-05 : [
        (N1,  (0.11,  0.015,  0.42)),
	(N2,  (0.12,  0.041,  0.31))],
	2E-04 : [
        (N1,  (0.14, 0.046,  0.49)), # XXX: way larger discrepancies than those reported in 2018-NG-GRA-D3.2_V1.2 (Table 53)
	(N2,  (0.096,  0.011,  0.31))],
	1E-04 : [
        (N1,  (0.14,  0.041,  0.48)),
	(N2,  (0.11,  0.017,  0.31))] }

  def tolerance_epsilon (tolerance):
    if tolerance > 1.0: return 0.01*tolerance
    tstr = "%.15f" % tolerance
    teps = 0.000000000000001
    j = len(tstr)-1
    while j > 0 and tstr[j] == '0':
      teps = teps * 10.
      j = j - 1
    return 0.0 if j == 1 else teps

  failed = False
  SEEK (solfec, stop)
  label = ('SX', 'SY', 'SZ')
  for (item, tole) in zip(data_p[step], tole_p[step]):
    pnt = item[0]
    stre = STRESS (body, pnt)
    for i in range(0,3):
      Value = -stre[i] # XXX: note the need to reverse sign
      Code_Aster_Ref = item[2][i]
      div =  1.0 if abs(Code_Aster_Ref) == 0.0 else abs(Code_Aster_Ref) 
      error = abs(Value-Code_Aster_Ref)/div 
      tolerance = tole[1][i] + tolerance_epsilon(tole[1][i])

      if error > tolerance:
	failed = True
	print '\b\b\b\bFAILED', 
	print '(Computed stress %s at point (%g, %g, %g) was %g' % (label[i], pnt[0], pnt[1], pnt[2], Value), 'while reference is %g)' % Code_Aster_Ref,
	print '-->(error %g, while tolerance %g)' % (error, tolerance)

  if not failed: print '\b\b\b\bPASSED'
