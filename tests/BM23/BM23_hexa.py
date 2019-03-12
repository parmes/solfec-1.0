# IMPACT DYNAMIC OF A BEAM ON A RIGID WALL

def create_simulation (step, stop, outstep):  
    Dir_out='out/tests/BM23_hexa_'+str(step)+'s'
    GEOMETRIC_EPSILON (1E-8)
    
    solfec=SOLFEC ('DYNAMIC', step, Dir_out)
    
    Nodes_FONDA=[
         -0.2000000E-01,-0.2000000E-01, 0.0000000E+00,
          0.0000000E+00,-0.2000000E-01, 0.0000000E+00,
         -0.2000000E-01, 0.0000000E+00, 0.0000000E+00,
          0.0000000E+00, 0.0000000E+00, 0.0000000E+00,
          0.2000000E-01,-0.2000000E-01, 0.0000000E+00,
          0.2000000E-01, 0.0000000E+00, 0.0000000E+00,
          0.0000000E+00, 0.2000000E-01, 0.0000000E+00,
          0.2000000E-01, 0.2000000E-01, 0.0000000E+00,
         -0.2000000E-01, 0.2000000E-01, 0.0000000E+00,
         -0.2000000E-01,-0.2000000E-01,-0.1000000E-01,
          0.0000000E+00,-0.2000000E-01,-0.1000000E-01,
         -0.2000000E-01, 0.0000000E+00,-0.1000000E-01,
          0.0000000E+00, 0.0000000E+00,-0.1000000E-01,
          0.2000000E-01,-0.2000000E-01,-0.1000000E-01,
          0.2000000E-01, 0.0000000E+00,-0.1000000E-01,
          0.0000000E+00, 0.2000000E-01,-0.1000000E-01,
          0.2000000E-01, 0.2000000E-01,-0.1000000E-01,
         -0.2000000E-01, 0.2000000E-01,-0.1000000E-01,]
    
    M35      = [8,    9,     10,     12,     11,     0,      1,      3,      2,        1]
    M36      = [8,   10,     13,     14,     12,     1,      4,      5,      3,        2]
    M37      = [8,   11,     12,     15,     17,     2,      3,      6,      8,        3]
    M38      = [8,   12,     14,     16,     15,     3,      5,      7,      6,        4]
    
    Mesh_FONDA = M35      + \
            M36      + \
            M37      + \
            M38 
            
    
    Nodes_POUTRE=[
         -0.1000000E-01,-0.1000000E-01, 0.2000000E-02,
         -0.1000000E-01,-0.1000000E-01, 0.1200000E-01,
         -0.1000000E-01,-0.1000000E-01, 0.2200000E-01,
         -0.1000000E-01,-0.1000000E-01, 0.3200000E-01,
         -0.1000000E-01,-0.1000000E-01, 0.4200000E-01,
         -0.1000000E-01,-0.1000000E-01, 0.5200000E-01,
         -0.1000000E-01,-0.1000000E-01, 0.6200000E-01,
         -0.1000000E-01,-0.1000000E-01, 0.7200000E-01,
         -0.1000000E-01,-0.1000000E-01, 0.8200000E-01,
         -0.1000000E-01,-0.1000000E-01, 0.9200000E-01,
         -0.1000000E-01,-0.1000000E-01, 0.1020000E+00,
         -0.1000000E-01,-0.1000000E-01, 0.1120000E+00,
         -0.1000000E-01,-0.1000000E-01, 0.1220000E+00,
         -0.1000000E-01,-0.1000000E-01, 0.1320000E+00,
         -0.1000000E-01,-0.1000000E-01, 0.1420000E+00,
         -0.1000000E-01,-0.1000000E-01, 0.1520000E+00,
         -0.1000000E-01,-0.1000000E-01, 0.1620000E+00,
         -0.1000000E-01,-0.1000000E-01, 0.1720000E+00,
         -0.1000000E-01,-0.1000000E-01, 0.1820000E+00,
         -0.1000000E-01,-0.1000000E-01, 0.1920000E+00,
          0.1000000E-01,-0.1000000E-01, 0.2000000E-02,
          0.1000000E-01,-0.1000000E-01, 0.1200000E-01,
          0.1000000E-01,-0.1000000E-01, 0.2200000E-01,
          0.1000000E-01,-0.1000000E-01, 0.3200000E-01,
          0.1000000E-01,-0.1000000E-01, 0.4200000E-01,
          0.1000000E-01,-0.1000000E-01, 0.5200000E-01,
          0.1000000E-01,-0.1000000E-01, 0.6200000E-01,
          0.1000000E-01,-0.1000000E-01, 0.7200000E-01,
          0.1000000E-01,-0.1000000E-01, 0.8200000E-01,
          0.1000000E-01,-0.1000000E-01, 0.9200000E-01,
          0.1000000E-01,-0.1000000E-01, 0.1020000E+00,
          0.1000000E-01,-0.1000000E-01, 0.1120000E+00,
          0.1000000E-01,-0.1000000E-01, 0.1220000E+00,
          0.1000000E-01,-0.1000000E-01, 0.1320000E+00,
          0.1000000E-01,-0.1000000E-01, 0.1420000E+00,
          0.1000000E-01,-0.1000000E-01, 0.1520000E+00,
          0.1000000E-01,-0.1000000E-01, 0.1620000E+00,
          0.1000000E-01,-0.1000000E-01, 0.1720000E+00,
          0.1000000E-01,-0.1000000E-01, 0.1820000E+00,
          0.1000000E-01,-0.1000000E-01, 0.1920000E+00,
          0.1000000E-01, 0.1000000E-01, 0.2000000E-02,
          0.1000000E-01, 0.1000000E-01, 0.1200000E-01,
          0.1000000E-01, 0.1000000E-01, 0.2200000E-01,
          0.1000000E-01, 0.1000000E-01, 0.3200000E-01,
          0.1000000E-01, 0.1000000E-01, 0.4200000E-01,
          0.1000000E-01, 0.1000000E-01, 0.5200000E-01,
          0.1000000E-01, 0.1000000E-01, 0.6200000E-01,
          0.1000000E-01, 0.1000000E-01, 0.7200000E-01,
          0.1000000E-01, 0.1000000E-01, 0.8200000E-01,
          0.1000000E-01, 0.1000000E-01, 0.9200000E-01,
          0.1000000E-01, 0.1000000E-01, 0.1020000E+00,
          0.1000000E-01, 0.1000000E-01, 0.1120000E+00,
          0.1000000E-01, 0.1000000E-01, 0.1220000E+00,
          0.1000000E-01, 0.1000000E-01, 0.1320000E+00,
          0.1000000E-01, 0.1000000E-01, 0.1420000E+00,
          0.1000000E-01, 0.1000000E-01, 0.1520000E+00,
          0.1000000E-01, 0.1000000E-01, 0.1620000E+00,
          0.1000000E-01, 0.1000000E-01, 0.1720000E+00,
          0.1000000E-01, 0.1000000E-01, 0.1820000E+00,
          0.1000000E-01, 0.1000000E-01, 0.1920000E+00,
         -0.1000000E-01, 0.1000000E-01, 0.2000000E-02,
         -0.1000000E-01, 0.1000000E-01, 0.1200000E-01,
         -0.1000000E-01, 0.1000000E-01, 0.2200000E-01,
         -0.1000000E-01, 0.1000000E-01, 0.3200000E-01,
         -0.1000000E-01, 0.1000000E-01, 0.4200000E-01,
         -0.1000000E-01, 0.1000000E-01, 0.5200000E-01,
         -0.1000000E-01, 0.1000000E-01, 0.6200000E-01,
         -0.1000000E-01, 0.1000000E-01, 0.7200000E-01,
         -0.1000000E-01, 0.1000000E-01, 0.8200000E-01,
         -0.1000000E-01, 0.1000000E-01, 0.9200000E-01,
         -0.1000000E-01, 0.1000000E-01, 0.1020000E+00,
         -0.1000000E-01, 0.1000000E-01, 0.1120000E+00,
         -0.1000000E-01, 0.1000000E-01, 0.1220000E+00,
         -0.1000000E-01, 0.1000000E-01, 0.1320000E+00,
         -0.1000000E-01, 0.1000000E-01, 0.1420000E+00,
         -0.1000000E-01, 0.1000000E-01, 0.1520000E+00,
         -0.1000000E-01, 0.1000000E-01, 0.1620000E+00,
         -0.1000000E-01, 0.1000000E-01, 0.1720000E+00,
         -0.1000000E-01, 0.1000000E-01, 0.1820000E+00,
         -0.1000000E-01, 0.1000000E-01, 0.1920000E+00,
         -0.1000000E-01,-0.1000000E-01, 0.2020000E+00,
          0.1000000E-01,-0.1000000E-01, 0.2020000E+00,
          0.1000000E-01, 0.1000000E-01, 0.2020000E+00,
         -0.1000000E-01, 0.1000000E-01, 0.2020000E+00,]
    #ELEM_TETRA4
    ############
    #ELEM_TETRA5
    ############
    #ELEM_TETRA6
    ############
    #ELEM_TETRA8
    ############
    M92      = [8,      0,     20,     40,     60,      1,     21,     41,     61,      1]
    M93      = [8,      1,     21,     41,     61,      2,     22,     42,     62,      2]
    M94      = [8,      2,     22,     42,     62,      3,     23,     43,     63,      3]
    M95      = [8,      3,     23,     43,     63,      4,     24,     44,     64,      4]
    M96      = [8,      4,     24,     44,     64,      5,     25,     45,     65,      5]
    M97      = [8,      5,     25,     45,     65,      6,     26,     46,     66,      6]
    M98      = [8,      6,     26,     46,     66,      7,     27,     47,     67,      7]
    M99      = [8,      7,     27,     47,     67,      8,     28,     48,     68,      8]
    M100     = [8,      8,     28,     48,     68,      9,     29,     49,     69,      9]
    M101     = [8,      9,     29,     49,     69,     10,     30,     50,     70,     10]
    M102     = [8,     10,     30,     50,     70,     11,     31,     51,     71,     11]
    M103     = [8,     11,     31,     51,     71,     12,     32,     52,     72,     12]
    M104     = [8,     12,     32,     52,     72,     13,     33,     53,     73,     13]
    M105     = [8,     13,     33,     53,     73,     14,     34,     54,     74,     14]
    M106     = [8,     14,     34,     54,     74,     15,     35,     55,     75,     15]
    M107     = [8,     15,     35,     55,     75,     16,     36,     56,     76,     16]
    M108     = [8,     16,     36,     56,     76,     17,     37,     57,     77,     17]
    M109     = [8,     17,     37,     57,     77,     18,     38,     58,     78,     18]
    M110     = [8,     18,     38,     58,     78,     19,     39,     59,     79,     19]
    M111     = [8,     19,     39,     59,     79,     80,     81,     82,     83,     20]
      #ELEM_TETRA4
      ############
      #ELEM_TETRA5
      ############
      #ELEM_TETRA6
      ############
      #ELEM_TETRA8
      ############
    Mesh_POUTRE = M92      + \
            M93      + \
            M94      + \
            M95      + \
            M96      + \
            M97      + \
            M98      + \
            M99      + \
            M100     + \
            M101     + \
            M102     + \
            M103     + \
            M104     + \
            M105     + \
            M106     + \
            M107     + \
            M108     + \
            M109     + \
            M110     + \
            M111           
    
                                          
    Msh_FONDA = MESH(Nodes_FONDA,Mesh_FONDA,1)
    Msh_POUTRE = MESH(Nodes_POUTRE,Mesh_POUTRE,2)
    
    Node=[]
    for i in Msh_POUTRE.nodes_on_surface(2):
      Node.append(Msh_POUTRE.node(i))
      
    
    FIX_XYZ=[]
    FIX_XYZ.append(Msh_FONDA.node(9))
    FIX_XYZ.append(Msh_FONDA.node(10))
    FIX_XYZ.append(Msh_FONDA.node(11))
    FIX_XYZ.append(Msh_FONDA.node(12))
    FIX_XYZ.append(Msh_FONDA.node(13))      
    FIX_XYZ.append(Msh_FONDA.node(14))
    FIX_XYZ.append(Msh_FONDA.node(15))
    FIX_XYZ.append(Msh_FONDA.node(16))
    FIX_XYZ.append(Msh_FONDA.node(17)) 
    
    FIX_XY_FONDA=[]
    FIX_XY_FONDA.append(Msh_FONDA.node(1))
    FIX_XY_FONDA.append(Msh_FONDA.node(2))
    FIX_XY_FONDA.append(Msh_FONDA.node(3))
    FIX_XY_FONDA.append(Msh_FONDA.node(5))
    FIX_XY_FONDA.append(Msh_FONDA.node(6))  
    
    
    FIX_XY_POUTRE=[]
    FIX_XY_POUTRE.append(Msh_POUTRE.node(0))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(1))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(2))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(3))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(4))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(5))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(6))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(7))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(8))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(9))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(10))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(11))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(12))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(13))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(14))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(15))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(16))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(17))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(18))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(19))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(20))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(21))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(22))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(23))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(24))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(25))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(26))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(27))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(28))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(29))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(30))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(31))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(32))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(33))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(34))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(35))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(36))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(37))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(38))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(39))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(40))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(41))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(42))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(43))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(44))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(45))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(46))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(47))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(48))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(49))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(50))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(51))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(52))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(53))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(54))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(55))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(56))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(57))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(58))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(59))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(60))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(61))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(62))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(63))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(64))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(65))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(66))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(67))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(68))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(69))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(70))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(71))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(72))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(73))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(74))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(75))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(76))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(77))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(78))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(79))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(80))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(81))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(82))
    FIX_XY_POUTRE.append(Msh_POUTRE.node(83))
    
    Point_C=Msh_POUTRE.node(0)
    
    N5  =(-1.00000000000000E-02, -1.00000000000000E-02,  2.00000000000000E-03)
    N100=(1.00000000000000E-02,  1.00000000000000E-02,  2.02000000000000E-01)
    
    No_Trait=[]  
    No_Trait.append(Msh_POUTRE.node(0))
                          
    ACIER = BULK_MATERIAL (solfec,
                          model = 'KIRCHHOFF',
                          young = 2.0E+11,
                          poisson = 0.,
                          density = 8.E+03)
    
    ACONTA = BULK_MATERIAL (solfec,
                          model = 'KIRCHHOFF',
                          young = 1.0E+16,
                          poisson = 0.0,
                          density = 1.E+00)                      
                                                
                          
    Bod_FONDA = BODY (solfec, 'OBSTACLE',Msh_FONDA,ACONTA)
    Bod_POUTRE = BODY (solfec, 'FINITE_ELEMENT',Msh_POUTRE,ACIER)
    
    
    for i in range (0,len(Node)):
      DISPLAY_POINT (Bod_POUTRE, Node[i], 'N%d'%i) # press 'D' in the viewer to see display points
    
    for fix in FIX_XYZ:
      FIX_POINT (Bod_FONDA, fix)    
        
    Fix_direc_X=(1.,0.,0.)
    Fix_direc_Y=(0.,1.,0.)
    Fix_direc_Z=(0.,0.,1.)  
        
    for fix in FIX_XY_FONDA:
      FIX_DIRECTION (Bod_FONDA, fix, Fix_direc_X) 
      FIX_DIRECTION (Bod_FONDA, fix, Fix_direc_Y) 
      
    for fix in FIX_XY_POUTRE:
      FIX_DIRECTION (Bod_POUTRE, fix, Fix_direc_X) 
      FIX_DIRECTION (Bod_POUTRE, fix, Fix_direc_Y)  
                              
    INITIAL_VELOCITY (Bod_POUTRE, (0., 0., -100), (0, 0, 0))
    
    gs= GAUSS_SEIDEL_SOLVER (1E-3 , 1000)
    gs_imp= GAUSS_SEIDEL_SOLVER (1E-6 , 1000, 1E-6)
    nt = NEWTON_SOLVER ()
    nt_imp= NEWTON_SOLVER (1E-9, 1000, delta = 1E-6)
    OUTPUT(solfec, outstep)
    if not VIEWER() and solfec.mode == 'WRITE':
      solfec.verbose = '%'
    RUN (solfec, gs_imp, stop)
    return (solfec, N5)
    
# data
step=2E-6
stop=0.00014
outstep=stop/10.

# calculate
(solfec, point) = create_simulation (2.E-6, stop, outstep)

# read results
if not VIEWER():
  (solfec, point) = create_simulation (step, stop, outstep)
  if solfec.mode != 'READ':
    print 'READ ERROR'
    import sys
    sys.exit(1)

  bod = solfec.bodies[0]
  his = HISTORY (solfec, [(bod, point, 'DZ'),(bod, point, 'VZ')],0, stop)

  #print his[1]
  #print his[2]

  # Solfec 24e88bd reference values
  refdz = [0.0, -0.001414213562373094, -0.0020501026727926534, -0.002050092655094636, -0.0020500955903457302,
    -0.002050096920308624, -0.002050095977140452, -0.002050188874306146, -0.0010794837493955337, 0.0002299152583814213]
  refvz = [-100.0, -99.99999999999798, 0.01749781148292584, 0.000845545473026732, -0.0017951314248847439,
    -8.96604709055282e-05, 0.000734002630039754, -0.04262958243722892, 103.56510775613184, 88.22422074488902]

  vals = [his[1], his[2]]
  refs = [refdz, refvz]
  failed = False
  for (val, ref) in zip (vals, refs):
    for (v, r) in zip (val, ref):
      if abs(v-r) > 1E-14:
	failed = True
	break

  if not failed:
    print '\b\b\b\bPASSED'
  else:
    print '\b\b\b\bFAILED'
