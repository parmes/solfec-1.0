# INFINITE PLATE WITH A CIRCULAR HOLE

def create_simulation (step, stop, outstep):
    Dir_out='out/tests/BM09_hexa_step_'+str(step)+'s'
    GEOMETRIC_EPSILON (1E-10)
  
    solfec=SOLFEC ('QUASI_STATIC', step, Dir_out)
    Nodes_M1=[
         -0.7406400E-01,-0.5343600E-01, 0.0000000E+00,
         -0.7204400E-01,-0.5343600E-01, 0.0000000E+00,
         -0.6950000E-01,-0.5343600E-01, 0.0000000E+00,
         -0.6600400E-01,-0.5343600E-01, 0.0000000E+00,
         -0.6116700E-01,-0.5343600E-01, 0.0000000E+00,
         -0.5461800E-01,-0.5343600E-01, 0.0000000E+00,
         -0.4600000E-01,-0.5343600E-01, 0.0000000E+00,
         -0.3496800E-01,-0.5343600E-01, 0.0000000E+00,
         -0.2118700E-01,-0.5343600E-01, 0.0000000E+00,
         -0.4326000E-02,-0.5343600E-01, 0.0000000E+00,
          0.1593600E-01,-0.5343600E-01, 0.0000000E+00,
         -0.7418700E-01,-0.5187300E-01, 0.0000000E+00,
         -0.7216400E-01,-0.5145900E-01, 0.0000000E+00,
         -0.6961700E-01,-0.5093800E-01, 0.0000000E+00,
         -0.6611700E-01,-0.5022200E-01, 0.0000000E+00,
         -0.6127300E-01,-0.4923100E-01, 0.0000000E+00,
         -0.5471500E-01,-0.4788900E-01, 0.0000000E+00,
         -0.4608500E-01,-0.4612400E-01, 0.0000000E+00,
         -0.3503800E-01,-0.4386400E-01, 0.0000000E+00,
         -0.2123800E-01,-0.4104100E-01, 0.0000000E+00,
         -0.4354000E-02,-0.3758700E-01, 0.0000000E+00,
          0.1593600E-01,-0.3343600E-01, 0.0000000E+00,
         -0.7455300E-01,-0.5034600E-01, 0.0000000E+00,
         -0.7252200E-01,-0.4951800E-01, 0.0000000E+00,
         -0.6996400E-01,-0.4847400E-01, 0.0000000E+00,
         -0.6645000E-01,-0.4704100E-01, 0.0000000E+00,
         -0.6158600E-01,-0.4505700E-01, 0.0000000E+00,
         -0.5500100E-01,-0.4237100E-01, 0.0000000E+00,
         -0.4633700E-01,-0.3883700E-01, 0.0000000E+00,
         -0.3524500E-01,-0.3431200E-01, 0.0000000E+00,
         -0.2138900E-01,-0.2866000E-01, 0.0000000E+00,
         -0.4436000E-02,-0.2174600E-01, 0.0000000E+00,
          0.1593600E-01,-0.1343600E-01, 0.0000000E+00,
         -0.7515400E-01,-0.4889700E-01, 0.0000000E+00,
         -0.7310900E-01,-0.4765200E-01, 0.0000000E+00,
         -0.7053400E-01,-0.4608400E-01, 0.0000000E+00,
         -0.6699600E-01,-0.4393000E-01, 0.0000000E+00,
         -0.6210100E-01,-0.4094900E-01, 0.0000000E+00,
         -0.5547200E-01,-0.3691300E-01, 0.0000000E+00,
         -0.4675000E-01,-0.3160300E-01, 0.0000000E+00,
         -0.3558500E-01,-0.2480500E-01, 0.0000000E+00,
         -0.2163600E-01,-0.1631200E-01, 0.0000000E+00,
         -0.4572000E-02,-0.5922000E-02, 0.0000000E+00,
          0.1593600E-01, 0.6564000E-02, 0.0000000E+00,
         -0.7597400E-01,-0.4755900E-01, 0.0000000E+00,
         -0.7391100E-01,-0.4589500E-01, 0.0000000E+00,
         -0.7131200E-01,-0.4379900E-01, 0.0000000E+00,
         -0.6774300E-01,-0.4092100E-01, 0.0000000E+00,
         -0.6280300E-01,-0.3693700E-01, 0.0000000E+00,
         -0.5611500E-01,-0.3154300E-01, 0.0000000E+00,
         -0.4731400E-01,-0.2444500E-01, 0.0000000E+00,
         -0.3604900E-01,-0.1536000E-01, 0.0000000E+00,
         -0.2197500E-01,-0.4010000E-02, 0.0000000E+00,
         -0.4756000E-02, 0.9877000E-02, 0.0000000E+00,
          0.1593600E-01, 0.2656400E-01, 0.0000000E+00,
         -0.7699300E-01,-0.4636500E-01, 0.0000000E+00,
         -0.7490700E-01,-0.4427900E-01, 0.0000000E+00,
         -0.7228000E-01,-0.4165200E-01, 0.0000000E+00,
         -0.6867100E-01,-0.3804300E-01, 0.0000000E+00,
         -0.6367700E-01,-0.3304800E-01, 0.0000000E+00,
         -0.5691400E-01,-0.2628600E-01, 0.0000000E+00,
         -0.4801600E-01,-0.1738700E-01, 0.0000000E+00,
         -0.3662500E-01,-0.5997000E-02, 0.0000000E+00,
         -0.2239500E-01, 0.8233000E-02, 0.0000000E+00,
         -0.4986000E-02, 0.2564200E-01, 0.0000000E+00,
          0.1593600E-01, 0.4656400E-01, 0.0000000E+00,
         -0.7818700E-01,-0.4534600E-01, 0.0000000E+00,
         -0.7652300E-01,-0.4328200E-01, 0.0000000E+00,
         -0.7442800E-01,-0.4068400E-01, 0.0000000E+00,
         -0.7154900E-01,-0.3711500E-01, 0.0000000E+00,
         -0.6756500E-01,-0.3217500E-01, 0.0000000E+00,
         -0.6217100E-01,-0.2548700E-01, 0.0000000E+00,
         -0.5507400E-01,-0.1668600E-01, 0.0000000E+00,
         -0.4598800E-01,-0.5420000E-02, 0.0000000E+00,
         -0.3463800E-01, 0.8654000E-02, 0.0000000E+00,
         -0.2075200E-01, 0.2587200E-01, 0.0000000E+00,
         -0.4064000E-02, 0.4656400E-01, 0.0000000E+00,
         -0.7952500E-01,-0.4452500E-01, 0.0000000E+00,
         -0.7828000E-01,-0.4248100E-01, 0.0000000E+00,
         -0.7671200E-01,-0.3990600E-01, 0.0000000E+00,
         -0.7455800E-01,-0.3636800E-01, 0.0000000E+00,
         -0.7157700E-01,-0.3147200E-01, 0.0000000E+00,
         -0.6754200E-01,-0.2484400E-01, 0.0000000E+00,
         -0.6223100E-01,-0.1612100E-01, 0.0000000E+00,
         -0.5543300E-01,-0.4956000E-02, 0.0000000E+00,
         -0.4694000E-01, 0.8992000E-02, 0.0000000E+00,
         -0.3655000E-01, 0.2605700E-01, 0.0000000E+00,
         -0.2406400E-01, 0.4656400E-01, 0.0000000E+00,
         -0.8097400E-01,-0.4392500E-01, 0.0000000E+00,
         -0.8014600E-01,-0.4189400E-01, 0.0000000E+00,
         -0.7910200E-01,-0.3933500E-01, 0.0000000E+00,
         -0.7766900E-01,-0.3582100E-01, 0.0000000E+00,
         -0.7568500E-01,-0.3095800E-01, 0.0000000E+00,
         -0.7299900E-01,-0.2437300E-01, 0.0000000E+00,
         -0.6946500E-01,-0.1570800E-01, 0.0000000E+00,
         -0.6494100E-01,-0.4617000E-02, 0.0000000E+00,
         -0.5928900E-01, 0.9240000E-02, 0.0000000E+00,
         -0.5237400E-01, 0.2619200E-01, 0.0000000E+00,
         -0.4406400E-01, 0.4656400E-01, 0.0000000E+00,
         -0.8250100E-01,-0.4355900E-01, 0.0000000E+00,
         -0.8208800E-01,-0.4153600E-01, 0.0000000E+00,
         -0.8156600E-01,-0.3898800E-01, 0.0000000E+00,
         -0.8085000E-01,-0.3548800E-01, 0.0000000E+00,
         -0.7985900E-01,-0.3064500E-01, 0.0000000E+00,
         -0.7851800E-01,-0.2408600E-01, 0.0000000E+00,
         -0.7675200E-01,-0.1545600E-01, 0.0000000E+00,
         -0.7449200E-01,-0.4410000E-02, 0.0000000E+00,
         -0.7166900E-01, 0.9391000E-02, 0.0000000E+00,
         -0.6821500E-01, 0.2627400E-01, 0.0000000E+00,
         -0.6406400E-01, 0.4656400E-01, 0.0000000E+00,
         -0.8406400E-01,-0.4343600E-01, 0.0000000E+00,
         -0.8406400E-01,-0.4141600E-01, 0.0000000E+00,
         -0.8406400E-01,-0.3887100E-01, 0.0000000E+00,
         -0.8406400E-01,-0.3537600E-01, 0.0000000E+00,
         -0.8406400E-01,-0.3053900E-01, 0.0000000E+00,
         -0.8406400E-01,-0.2399000E-01, 0.0000000E+00,
         -0.8406400E-01,-0.1537200E-01, 0.0000000E+00,
         -0.8406400E-01,-0.4340000E-02, 0.0000000E+00,
         -0.8406400E-01, 0.9442000E-02, 0.0000000E+00,
         -0.8406400E-01, 0.2630200E-01, 0.0000000E+00,
         -0.8406400E-01, 0.4656400E-01, 0.0000000E+00,
         -0.7406400E-01,-0.5343600E-01, 0.1000000E-02,
         -0.7204400E-01,-0.5343600E-01, 0.1000000E-02,
         -0.6950000E-01,-0.5343600E-01, 0.1000000E-02,
         -0.6600400E-01,-0.5343600E-01, 0.1000000E-02,
         -0.6116700E-01,-0.5343600E-01, 0.1000000E-02,
         -0.5461800E-01,-0.5343600E-01, 0.1000000E-02,
         -0.4600000E-01,-0.5343600E-01, 0.1000000E-02,
         -0.3496800E-01,-0.5343600E-01, 0.1000000E-02,
         -0.2118700E-01,-0.5343600E-01, 0.1000000E-02,
         -0.4326000E-02,-0.5343600E-01, 0.1000000E-02,
          0.1593600E-01,-0.5343600E-01, 0.1000000E-02,
         -0.7418700E-01,-0.5187300E-01, 0.1000000E-02,
         -0.7216400E-01,-0.5145900E-01, 0.1000000E-02,
         -0.6961700E-01,-0.5093800E-01, 0.1000000E-02,
         -0.6611700E-01,-0.5022200E-01, 0.1000000E-02,
         -0.6127300E-01,-0.4923100E-01, 0.1000000E-02,
         -0.5471500E-01,-0.4788900E-01, 0.1000000E-02,
         -0.4608500E-01,-0.4612400E-01, 0.1000000E-02,
         -0.3503800E-01,-0.4386400E-01, 0.1000000E-02,
         -0.2123800E-01,-0.4104100E-01, 0.1000000E-02,
         -0.4354000E-02,-0.3758700E-01, 0.1000000E-02,
          0.1593600E-01,-0.3343600E-01, 0.1000000E-02,
         -0.7455300E-01,-0.5034600E-01, 0.1000000E-02,
         -0.7252200E-01,-0.4951800E-01, 0.1000000E-02,
         -0.6996400E-01,-0.4847400E-01, 0.1000000E-02,
         -0.6645000E-01,-0.4704100E-01, 0.1000000E-02,
         -0.6158600E-01,-0.4505700E-01, 0.1000000E-02,
         -0.5500100E-01,-0.4237100E-01, 0.1000000E-02,
         -0.4633700E-01,-0.3883700E-01, 0.1000000E-02,
         -0.3524500E-01,-0.3431200E-01, 0.1000000E-02,
         -0.2138900E-01,-0.2866000E-01, 0.1000000E-02,
         -0.4436000E-02,-0.2174600E-01, 0.1000000E-02,
          0.1593600E-01,-0.1343600E-01, 0.1000000E-02,
         -0.7515400E-01,-0.4889700E-01, 0.1000000E-02,
         -0.7310900E-01,-0.4765200E-01, 0.1000000E-02,
         -0.7053400E-01,-0.4608400E-01, 0.1000000E-02,
         -0.6699600E-01,-0.4393000E-01, 0.1000000E-02,
         -0.6210100E-01,-0.4094900E-01, 0.1000000E-02,
         -0.5547200E-01,-0.3691300E-01, 0.1000000E-02,
         -0.4675000E-01,-0.3160300E-01, 0.1000000E-02,
         -0.3558500E-01,-0.2480500E-01, 0.1000000E-02,
         -0.2163600E-01,-0.1631200E-01, 0.1000000E-02,
         -0.4572000E-02,-0.5922000E-02, 0.1000000E-02,
          0.1593600E-01, 0.6564000E-02, 0.1000000E-02,
         -0.7597400E-01,-0.4755900E-01, 0.1000000E-02,
         -0.7391100E-01,-0.4589500E-01, 0.1000000E-02,
         -0.7131200E-01,-0.4379900E-01, 0.1000000E-02,
         -0.6774300E-01,-0.4092100E-01, 0.1000000E-02,
         -0.6280300E-01,-0.3693700E-01, 0.1000000E-02,
         -0.5611500E-01,-0.3154300E-01, 0.1000000E-02,
         -0.4731400E-01,-0.2444500E-01, 0.1000000E-02,
         -0.3604900E-01,-0.1536000E-01, 0.1000000E-02,
         -0.2197500E-01,-0.4010000E-02, 0.1000000E-02,
         -0.4756000E-02, 0.9877000E-02, 0.1000000E-02,
          0.1593600E-01, 0.2656400E-01, 0.1000000E-02,
         -0.7699300E-01,-0.4636500E-01, 0.1000000E-02,
         -0.7490700E-01,-0.4427900E-01, 0.1000000E-02,
         -0.7228000E-01,-0.4165200E-01, 0.1000000E-02,
         -0.6867100E-01,-0.3804300E-01, 0.1000000E-02,
         -0.6367700E-01,-0.3304800E-01, 0.1000000E-02,
         -0.5691400E-01,-0.2628600E-01, 0.1000000E-02,
         -0.4801600E-01,-0.1738700E-01, 0.1000000E-02,
         -0.3662500E-01,-0.5997000E-02, 0.1000000E-02,
         -0.2239500E-01, 0.8233000E-02, 0.1000000E-02,
         -0.4986000E-02, 0.2564200E-01, 0.1000000E-02,
          0.1593600E-01, 0.4656400E-01, 0.1000000E-02,
         -0.7818700E-01,-0.4534600E-01, 0.1000000E-02,
         -0.7652300E-01,-0.4328200E-01, 0.1000000E-02,
         -0.7442800E-01,-0.4068400E-01, 0.1000000E-02,
         -0.7154900E-01,-0.3711500E-01, 0.1000000E-02,
         -0.6756500E-01,-0.3217500E-01, 0.1000000E-02,
         -0.6217100E-01,-0.2548700E-01, 0.1000000E-02,
         -0.5507400E-01,-0.1668600E-01, 0.1000000E-02,
         -0.4598800E-01,-0.5420000E-02, 0.1000000E-02,
         -0.3463800E-01, 0.8654000E-02, 0.1000000E-02,
         -0.2075200E-01, 0.2587200E-01, 0.1000000E-02,
         -0.4064000E-02, 0.4656400E-01, 0.1000000E-02,
         -0.7952500E-01,-0.4452500E-01, 0.1000000E-02,
         -0.7828000E-01,-0.4248100E-01, 0.1000000E-02,
         -0.7671200E-01,-0.3990600E-01, 0.1000000E-02,
         -0.7455800E-01,-0.3636800E-01, 0.1000000E-02,
         -0.7157700E-01,-0.3147200E-01, 0.1000000E-02,
         -0.6754200E-01,-0.2484400E-01, 0.1000000E-02,
         -0.6223100E-01,-0.1612100E-01, 0.1000000E-02,
         -0.5543300E-01,-0.4956000E-02, 0.1000000E-02,
         -0.4694000E-01, 0.8992000E-02, 0.1000000E-02,
         -0.3655000E-01, 0.2605700E-01, 0.1000000E-02,
         -0.2406400E-01, 0.4656400E-01, 0.1000000E-02,
         -0.8097400E-01,-0.4392500E-01, 0.1000000E-02,
         -0.8014600E-01,-0.4189400E-01, 0.1000000E-02,
         -0.7910200E-01,-0.3933500E-01, 0.1000000E-02,
         -0.7766900E-01,-0.3582100E-01, 0.1000000E-02,
         -0.7568500E-01,-0.3095800E-01, 0.1000000E-02,
         -0.7299900E-01,-0.2437300E-01, 0.1000000E-02,
         -0.6946500E-01,-0.1570800E-01, 0.1000000E-02,
         -0.6494100E-01,-0.4617000E-02, 0.1000000E-02,
         -0.5928900E-01, 0.9240000E-02, 0.1000000E-02,
         -0.5237400E-01, 0.2619200E-01, 0.1000000E-02,
         -0.4406400E-01, 0.4656400E-01, 0.1000000E-02,
         -0.8250100E-01,-0.4355900E-01, 0.1000000E-02,
         -0.8208800E-01,-0.4153600E-01, 0.1000000E-02,
         -0.8156600E-01,-0.3898800E-01, 0.1000000E-02,
         -0.8085000E-01,-0.3548800E-01, 0.1000000E-02,
         -0.7985900E-01,-0.3064500E-01, 0.1000000E-02,
         -0.7851800E-01,-0.2408600E-01, 0.1000000E-02,
         -0.7675200E-01,-0.1545600E-01, 0.1000000E-02,
         -0.7449200E-01,-0.4410000E-02, 0.1000000E-02,
         -0.7166900E-01, 0.9391000E-02, 0.1000000E-02,
         -0.6821500E-01, 0.2627400E-01, 0.1000000E-02,
         -0.6406400E-01, 0.4656400E-01, 0.1000000E-02,
         -0.8406400E-01,-0.4343600E-01, 0.1000000E-02,
         -0.8406400E-01,-0.4141600E-01, 0.1000000E-02,
         -0.8406400E-01,-0.3887100E-01, 0.1000000E-02,
         -0.8406400E-01,-0.3537600E-01, 0.1000000E-02,
         -0.8406400E-01,-0.3053900E-01, 0.1000000E-02,
         -0.8406400E-01,-0.2399000E-01, 0.1000000E-02,
         -0.8406400E-01,-0.1537200E-01, 0.1000000E-02,
         -0.8406400E-01,-0.4340000E-02, 0.1000000E-02,
         -0.8406400E-01, 0.9442000E-02, 0.1000000E-02,
         -0.8406400E-01, 0.2630200E-01, 0.1000000E-02,
         -0.8406400E-01, 0.4656400E-01, 0.1000000E-02,]
    #ELEM_TETRA4
    ############
    #ELEM_TETRA5
    ############
    #ELEM_TETRA6
    ############
    #ELEM_TETRA8
    ############
    M101     = [8,      0,      1,     12,     11,    121,    122,    133,    132,      1]
    M102     = [8,      1,      2,     13,     12,    122,    123,    134,    133,      2]
    M103     = [8,      2,      3,     14,     13,    123,    124,    135,    134,      3]
    M104     = [8,      3,      4,     15,     14,    124,    125,    136,    135,      4]
    M105     = [8,      4,      5,     16,     15,    125,    126,    137,    136,      5]
    M106     = [8,      5,      6,     17,     16,    126,    127,    138,    137,      6]
    M107     = [8,      6,      7,     18,     17,    127,    128,    139,    138,      7]
    M108     = [8,      7,      8,     19,     18,    128,    129,    140,    139,      8]
    M109     = [8,      8,      9,     20,     19,    129,    130,    141,    140,      9]
    M110     = [8,      9,     10,     21,     20,    130,    131,    142,    141,     10]
    M111     = [8,     11,     12,     23,     22,    132,    133,    144,    143,     11]
    M112     = [8,     12,     13,     24,     23,    133,    134,    145,    144,     12]
    M113     = [8,     13,     14,     25,     24,    134,    135,    146,    145,     13]
    M114     = [8,     14,     15,     26,     25,    135,    136,    147,    146,     14]
    M115     = [8,     15,     16,     27,     26,    136,    137,    148,    147,     15]
    M116     = [8,     16,     17,     28,     27,    137,    138,    149,    148,     16]
    M117     = [8,     17,     18,     29,     28,    138,    139,    150,    149,     17]
    M118     = [8,     18,     19,     30,     29,    139,    140,    151,    150,     18]
    M119     = [8,     19,     20,     31,     30,    140,    141,    152,    151,     19]
    M120     = [8,     20,     21,     32,     31,    141,    142,    153,    152,     20]
    M121     = [8,     22,     23,     34,     33,    143,    144,    155,    154,     21]
    M122     = [8,     23,     24,     35,     34,    144,    145,    156,    155,     22]
    M123     = [8,     24,     25,     36,     35,    145,    146,    157,    156,     23]
    M124     = [8,     25,     26,     37,     36,    146,    147,    158,    157,     24]
    M125     = [8,     26,     27,     38,     37,    147,    148,    159,    158,     25]
    M126     = [8,     27,     28,     39,     38,    148,    149,    160,    159,     26]
    M127     = [8,     28,     29,     40,     39,    149,    150,    161,    160,     27]
    M128     = [8,     29,     30,     41,     40,    150,    151,    162,    161,     28]
    M129     = [8,     30,     31,     42,     41,    151,    152,    163,    162,     29]
    M130     = [8,     31,     32,     43,     42,    152,    153,    164,    163,     30]
    M131     = [8,     33,     34,     45,     44,    154,    155,    166,    165,     31]
    M132     = [8,     34,     35,     46,     45,    155,    156,    167,    166,     32]
    M133     = [8,     35,     36,     47,     46,    156,    157,    168,    167,     33]
    M134     = [8,     36,     37,     48,     47,    157,    158,    169,    168,     34]
    M135     = [8,     37,     38,     49,     48,    158,    159,    170,    169,     35]
    M136     = [8,     38,     39,     50,     49,    159,    160,    171,    170,     36]
    M137     = [8,     39,     40,     51,     50,    160,    161,    172,    171,     37]
    M138     = [8,     40,     41,     52,     51,    161,    162,    173,    172,     38]
    M139     = [8,     41,     42,     53,     52,    162,    163,    174,    173,     39]
    M140     = [8,     42,     43,     54,     53,    163,    164,    175,    174,     40]
    M141     = [8,     44,     45,     56,     55,    165,    166,    177,    176,     41]
    M142     = [8,     45,     46,     57,     56,    166,    167,    178,    177,     42]
    M143     = [8,     46,     47,     58,     57,    167,    168,    179,    178,     43]
    M144     = [8,     47,     48,     59,     58,    168,    169,    180,    179,     44]
    M145     = [8,     48,     49,     60,     59,    169,    170,    181,    180,     45]
    M146     = [8,     49,     50,     61,     60,    170,    171,    182,    181,     46]
    M147     = [8,     50,     51,     62,     61,    171,    172,    183,    182,     47]
    M148     = [8,     51,     52,     63,     62,    172,    173,    184,    183,     48]
    M149     = [8,     52,     53,     64,     63,    173,    174,    185,    184,     49]
    M150     = [8,     53,     54,     65,     64,    174,    175,    186,    185,     50]
    M151     = [8,     55,     56,     67,     66,    176,    177,    188,    187,     51]
    M152     = [8,     56,     57,     68,     67,    177,    178,    189,    188,     52]
    M153     = [8,     57,     58,     69,     68,    178,    179,    190,    189,     53]
    M154     = [8,     58,     59,     70,     69,    179,    180,    191,    190,     54]
    M155     = [8,     59,     60,     71,     70,    180,    181,    192,    191,     55]
    M156     = [8,     60,     61,     72,     71,    181,    182,    193,    192,     56]
    M157     = [8,     61,     62,     73,     72,    182,    183,    194,    193,     57]
    M158     = [8,     62,     63,     74,     73,    183,    184,    195,    194,     58]
    M159     = [8,     63,     64,     75,     74,    184,    185,    196,    195,     59]
    M160     = [8,     64,     65,     76,     75,    185,    186,    197,    196,     60]
    M161     = [8,     66,     67,     78,     77,    187,    188,    199,    198,     61]
    M162     = [8,     67,     68,     79,     78,    188,    189,    200,    199,     62]
    M163     = [8,     68,     69,     80,     79,    189,    190,    201,    200,     63]
    M164     = [8,     69,     70,     81,     80,    190,    191,    202,    201,     64]
    M165     = [8,     70,     71,     82,     81,    191,    192,    203,    202,     65]
    M166     = [8,     71,     72,     83,     82,    192,    193,    204,    203,     66]
    M167     = [8,     72,     73,     84,     83,    193,    194,    205,    204,     67]
    M168     = [8,     73,     74,     85,     84,    194,    195,    206,    205,     68]
    M169     = [8,     74,     75,     86,     85,    195,    196,    207,    206,     69]
    M170     = [8,     75,     76,     87,     86,    196,    197,    208,    207,     70]
    M171     = [8,     77,     78,     89,     88,    198,    199,    210,    209,     71]
    M172     = [8,     78,     79,     90,     89,    199,    200,    211,    210,     72]
    M173     = [8,     79,     80,     91,     90,    200,    201,    212,    211,     73]
    M174     = [8,     80,     81,     92,     91,    201,    202,    213,    212,     74]
    M175     = [8,     81,     82,     93,     92,    202,    203,    214,    213,     75]
    M176     = [8,     82,     83,     94,     93,    203,    204,    215,    214,     76]
    M177     = [8,     83,     84,     95,     94,    204,    205,    216,    215,     77]
    M178     = [8,     84,     85,     96,     95,    205,    206,    217,    216,     78]
    M179     = [8,     85,     86,     97,     96,    206,    207,    218,    217,     79]
    M180     = [8,     86,     87,     98,     97,    207,    208,    219,    218,     80]
    M181     = [8,     88,     89,    100,     99,    209,    210,    221,    220,     81]
    M182     = [8,     89,     90,    101,    100,    210,    211,    222,    221,     82]
    M183     = [8,     90,     91,    102,    101,    211,    212,    223,    222,     83]
    M184     = [8,     91,     92,    103,    102,    212,    213,    224,    223,     84]
    M185     = [8,     92,     93,    104,    103,    213,    214,    225,    224,     85]
    M186     = [8,     93,     94,    105,    104,    214,    215,    226,    225,     86]
    M187     = [8,     94,     95,    106,    105,    215,    216,    227,    226,     87]
    M188     = [8,     95,     96,    107,    106,    216,    217,    228,    227,     88]
    M189     = [8,     96,     97,    108,    107,    217,    218,    229,    228,     89]
    M190     = [8,     97,     98,    109,    108,    218,    219,    230,    229,     90]
    M191     = [8,     99,    100,    111,    110,    220,    221,    232,    231,     91]
    M192     = [8,    100,    101,    112,    111,    221,    222,    233,    232,     92]
    M193     = [8,    101,    102,    113,    112,    222,    223,    234,    233,     93]
    M194     = [8,    102,    103,    114,    113,    223,    224,    235,    234,     94]
    M195     = [8,    103,    104,    115,    114,    224,    225,    236,    235,     95]
    M196     = [8,    104,    105,    116,    115,    225,    226,    237,    236,     96]
    M197     = [8,    105,    106,    117,    116,    226,    227,    238,    237,     97]
    M198     = [8,    106,    107,    118,    117,    227,    228,    239,    238,     98]
    M199     = [8,    107,    108,    119,    118,    228,    229,    240,    239,     99]
    M200     = [8,    108,    109,    120,    119,    229,    230,    241,    240,    100]
      #ELEM_TETRA4
      ############
      #ELEM_TETRA5
      ############
      #ELEM_TETRA6
      ############
      #ELEM_TETRA8
      ############
    Mesh_M1 = M101     + \
            M102     + \
            M103     + \
            M104     + \
            M105     + \
            M106     + \
            M107     + \
            M108     + \
            M109     + \
            M110     + \
            M111     + \
            M112     + \
            M113     + \
            M114     + \
            M115     + \
            M116     + \
            M117     + \
            M118     + \
            M119     + \
            M120     + \
            M121     + \
            M122     + \
            M123     + \
            M124     + \
            M125     + \
            M126     + \
            M127     + \
            M128     + \
            M129     + \
            M130     + \
            M131     + \
            M132     + \
            M133     + \
            M134     + \
            M135     + \
            M136     + \
            M137     + \
            M138     + \
            M139     + \
            M140     + \
            M141     + \
            M142     + \
            M143     + \
            M144     + \
            M145     + \
            M146     + \
            M147     + \
            M148     + \
            M149     + \
            M150     + \
            M151     + \
            M152     + \
            M153     + \
            M154     + \
            M155     + \
            M156     + \
            M157     + \
            M158     + \
            M159     + \
            M160     + \
            M161     + \
            M162     + \
            M163     + \
            M164     + \
            M165     + \
            M166     + \
            M167     + \
            M168     + \
            M169     + \
            M170     + \
            M171     + \
            M172     + \
            M173     + \
            M174     + \
            M175     + \
            M176     + \
            M177     + \
            M178     + \
            M179     + \
            M180     + \
            M181     + \
            M182     + \
            M183     + \
            M184     + \
            M185     + \
            M186     + \
            M187     + \
            M188     + \
            M189     + \
            M190     + \
            M191     + \
            M192     + \
            M193     + \
            M194     + \
            M195     + \
            M196     + \
            M197     + \
            M198     + \
            M199     + \
            M200 
            
                   
                                           
    Msh_M1 = MESH(Nodes_M1,Mesh_M1,1)

    Mesh_Final = Msh_M1

    Node=[]
    for i in Mesh_Final.nodes_on_surface(1):
      Node.append(Mesh_Final.node(i))
     

    FIX_X=[]
    FIX_X.append(Mesh_Final.node(110))
    FIX_X.append(Mesh_Final.node(111))
    FIX_X.append(Mesh_Final.node(112))
    FIX_X.append(Mesh_Final.node(113))
    FIX_X.append(Mesh_Final.node(114))
    FIX_X.append(Mesh_Final.node(115))
    FIX_X.append(Mesh_Final.node(116))
    FIX_X.append(Mesh_Final.node(117))
    FIX_X.append(Mesh_Final.node(118))
    FIX_X.append(Mesh_Final.node(119))
    FIX_X.append(Mesh_Final.node(120))
    FIX_X.append(Mesh_Final.node(231))
    FIX_X.append(Mesh_Final.node(232))
    FIX_X.append(Mesh_Final.node(233))
    FIX_X.append(Mesh_Final.node(234))
    FIX_X.append(Mesh_Final.node(235))
    FIX_X.append(Mesh_Final.node(236))
    FIX_X.append(Mesh_Final.node(237))
    FIX_X.append(Mesh_Final.node(238))
    FIX_X.append(Mesh_Final.node(239))
    FIX_X.append(Mesh_Final.node(240))
    FIX_X.append(Mesh_Final.node(241))


    FIX_Y=[]
    FIX_Y.append(Mesh_Final.node(0))
    FIX_Y.append(Mesh_Final.node(1))
    FIX_Y.append(Mesh_Final.node(2))
    FIX_Y.append(Mesh_Final.node(3))
    FIX_Y.append(Mesh_Final.node(4))
    FIX_Y.append(Mesh_Final.node(5))
    FIX_Y.append(Mesh_Final.node(6))
    FIX_Y.append(Mesh_Final.node(7))
    FIX_Y.append(Mesh_Final.node(8))
    FIX_Y.append(Mesh_Final.node(9))
    FIX_Y.append(Mesh_Final.node(10))
    FIX_Y.append(Mesh_Final.node(121))
    FIX_Y.append(Mesh_Final.node(122))
    FIX_Y.append(Mesh_Final.node(123))
    FIX_Y.append(Mesh_Final.node(124))
    FIX_Y.append(Mesh_Final.node(125))
    FIX_Y.append(Mesh_Final.node(126))
    FIX_Y.append(Mesh_Final.node(127))
    FIX_Y.append(Mesh_Final.node(128))
    FIX_Y.append(Mesh_Final.node(129))
    FIX_Y.append(Mesh_Final.node(130))
    FIX_Y.append(Mesh_Final.node(131))

    FIX_Z=[]
    FIX_Z.append(Mesh_Final.node(11))
    FIX_Z.append(Mesh_Final.node(12))
    FIX_Z.append(Mesh_Final.node(13))
    FIX_Z.append(Mesh_Final.node(14))
    FIX_Z.append(Mesh_Final.node(15))
    FIX_Z.append(Mesh_Final.node(16))
    FIX_Z.append(Mesh_Final.node(17))
    FIX_Z.append(Mesh_Final.node(18))
    FIX_Z.append(Mesh_Final.node(19))
    FIX_Z.append(Mesh_Final.node(20))
    FIX_Z.append(Mesh_Final.node(21))
    FIX_Z.append(Mesh_Final.node(22))
    FIX_Z.append(Mesh_Final.node(23))
    FIX_Z.append(Mesh_Final.node(24))
    FIX_Z.append(Mesh_Final.node(25))
    FIX_Z.append(Mesh_Final.node(26))
    FIX_Z.append(Mesh_Final.node(27))
    FIX_Z.append(Mesh_Final.node(28))
    FIX_Z.append(Mesh_Final.node(29))
    FIX_Z.append(Mesh_Final.node(30))
    FIX_Z.append(Mesh_Final.node(31))
    FIX_Z.append(Mesh_Final.node(32))
    FIX_Z.append(Mesh_Final.node(33))
    FIX_Z.append(Mesh_Final.node(34))
    FIX_Z.append(Mesh_Final.node(35))
    FIX_Z.append(Mesh_Final.node(36))
    FIX_Z.append(Mesh_Final.node(37))
    FIX_Z.append(Mesh_Final.node(38))
    FIX_Z.append(Mesh_Final.node(39))
    FIX_Z.append(Mesh_Final.node(40))
    FIX_Z.append(Mesh_Final.node(41))
    FIX_Z.append(Mesh_Final.node(42))
    FIX_Z.append(Mesh_Final.node(43))
    FIX_Z.append(Mesh_Final.node(44))
    FIX_Z.append(Mesh_Final.node(45))
    FIX_Z.append(Mesh_Final.node(46))
    FIX_Z.append(Mesh_Final.node(47))
    FIX_Z.append(Mesh_Final.node(48))
    FIX_Z.append(Mesh_Final.node(49))
    FIX_Z.append(Mesh_Final.node(50))
    FIX_Z.append(Mesh_Final.node(51))
    FIX_Z.append(Mesh_Final.node(52))
    FIX_Z.append(Mesh_Final.node(53))
    FIX_Z.append(Mesh_Final.node(54))
    FIX_Z.append(Mesh_Final.node(55))
    FIX_Z.append(Mesh_Final.node(56))
    FIX_Z.append(Mesh_Final.node(57))
    FIX_Z.append(Mesh_Final.node(58))
    FIX_Z.append(Mesh_Final.node(59))
    FIX_Z.append(Mesh_Final.node(60))
    FIX_Z.append(Mesh_Final.node(61))
    FIX_Z.append(Mesh_Final.node(62))
    FIX_Z.append(Mesh_Final.node(63))
    FIX_Z.append(Mesh_Final.node(64))
    FIX_Z.append(Mesh_Final.node(65))
    FIX_Z.append(Mesh_Final.node(66))
    FIX_Z.append(Mesh_Final.node(67))
    FIX_Z.append(Mesh_Final.node(68))
    FIX_Z.append(Mesh_Final.node(69))
    FIX_Z.append(Mesh_Final.node(70))
    FIX_Z.append(Mesh_Final.node(71))
    FIX_Z.append(Mesh_Final.node(72))
    FIX_Z.append(Mesh_Final.node(73))
    FIX_Z.append(Mesh_Final.node(74))
    FIX_Z.append(Mesh_Final.node(75))
    FIX_Z.append(Mesh_Final.node(76))
    FIX_Z.append(Mesh_Final.node(77))
    FIX_Z.append(Mesh_Final.node(78))
    FIX_Z.append(Mesh_Final.node(79))
    FIX_Z.append(Mesh_Final.node(80))
    FIX_Z.append(Mesh_Final.node(81))
    FIX_Z.append(Mesh_Final.node(82))
    FIX_Z.append(Mesh_Final.node(83))
    FIX_Z.append(Mesh_Final.node(84))
    FIX_Z.append(Mesh_Final.node(85))
    FIX_Z.append(Mesh_Final.node(86))
    FIX_Z.append(Mesh_Final.node(87))
    FIX_Z.append(Mesh_Final.node(88))
    FIX_Z.append(Mesh_Final.node(89))
    FIX_Z.append(Mesh_Final.node(90))
    FIX_Z.append(Mesh_Final.node(91))
    FIX_Z.append(Mesh_Final.node(92))
    FIX_Z.append(Mesh_Final.node(93))
    FIX_Z.append(Mesh_Final.node(94))
    FIX_Z.append(Mesh_Final.node(95))
    FIX_Z.append(Mesh_Final.node(96))
    FIX_Z.append(Mesh_Final.node(97))
    FIX_Z.append(Mesh_Final.node(98))
    FIX_Z.append(Mesh_Final.node(99))
    FIX_Z.append(Mesh_Final.node(100))
    FIX_Z.append(Mesh_Final.node(101))
    FIX_Z.append(Mesh_Final.node(102))
    FIX_Z.append(Mesh_Final.node(103))
    FIX_Z.append(Mesh_Final.node(104))
    FIX_Z.append(Mesh_Final.node(105))
    FIX_Z.append(Mesh_Final.node(106))
    FIX_Z.append(Mesh_Final.node(107))
    FIX_Z.append(Mesh_Final.node(108))
    FIX_Z.append(Mesh_Final.node(109))
    FIX_Z.append(Mesh_Final.node(132))
    FIX_Z.append(Mesh_Final.node(133))
    FIX_Z.append(Mesh_Final.node(134))
    FIX_Z.append(Mesh_Final.node(135))
    FIX_Z.append(Mesh_Final.node(136))
    FIX_Z.append(Mesh_Final.node(137))
    FIX_Z.append(Mesh_Final.node(138))
    FIX_Z.append(Mesh_Final.node(139))
    FIX_Z.append(Mesh_Final.node(140))
    FIX_Z.append(Mesh_Final.node(141))
    FIX_Z.append(Mesh_Final.node(142))
    FIX_Z.append(Mesh_Final.node(143))
    FIX_Z.append(Mesh_Final.node(144))
    FIX_Z.append(Mesh_Final.node(145))
    FIX_Z.append(Mesh_Final.node(146))
    FIX_Z.append(Mesh_Final.node(147))
    FIX_Z.append(Mesh_Final.node(148))
    FIX_Z.append(Mesh_Final.node(149))
    FIX_Z.append(Mesh_Final.node(150))
    FIX_Z.append(Mesh_Final.node(151))
    FIX_Z.append(Mesh_Final.node(152))
    FIX_Z.append(Mesh_Final.node(153))
    FIX_Z.append(Mesh_Final.node(154))
    FIX_Z.append(Mesh_Final.node(155))
    FIX_Z.append(Mesh_Final.node(156))
    FIX_Z.append(Mesh_Final.node(157))
    FIX_Z.append(Mesh_Final.node(158))
    FIX_Z.append(Mesh_Final.node(159))
    FIX_Z.append(Mesh_Final.node(160))
    FIX_Z.append(Mesh_Final.node(161))
    FIX_Z.append(Mesh_Final.node(162))
    FIX_Z.append(Mesh_Final.node(163))
    FIX_Z.append(Mesh_Final.node(164))
    FIX_Z.append(Mesh_Final.node(165))
    FIX_Z.append(Mesh_Final.node(166))
    FIX_Z.append(Mesh_Final.node(167))
    FIX_Z.append(Mesh_Final.node(168))
    FIX_Z.append(Mesh_Final.node(169))
    FIX_Z.append(Mesh_Final.node(170))
    FIX_Z.append(Mesh_Final.node(171))
    FIX_Z.append(Mesh_Final.node(172))
    FIX_Z.append(Mesh_Final.node(173))
    FIX_Z.append(Mesh_Final.node(174))
    FIX_Z.append(Mesh_Final.node(175))
    FIX_Z.append(Mesh_Final.node(176))
    FIX_Z.append(Mesh_Final.node(177))
    FIX_Z.append(Mesh_Final.node(178))
    FIX_Z.append(Mesh_Final.node(179))
    FIX_Z.append(Mesh_Final.node(180))
    FIX_Z.append(Mesh_Final.node(181))
    FIX_Z.append(Mesh_Final.node(182))
    FIX_Z.append(Mesh_Final.node(183))
    FIX_Z.append(Mesh_Final.node(184))
    FIX_Z.append(Mesh_Final.node(185))
    FIX_Z.append(Mesh_Final.node(186))
    FIX_Z.append(Mesh_Final.node(187))
    FIX_Z.append(Mesh_Final.node(188))
    FIX_Z.append(Mesh_Final.node(189))
    FIX_Z.append(Mesh_Final.node(190))
    FIX_Z.append(Mesh_Final.node(191))
    FIX_Z.append(Mesh_Final.node(192))
    FIX_Z.append(Mesh_Final.node(193))
    FIX_Z.append(Mesh_Final.node(194))
    FIX_Z.append(Mesh_Final.node(195))
    FIX_Z.append(Mesh_Final.node(196))
    FIX_Z.append(Mesh_Final.node(197))
    FIX_Z.append(Mesh_Final.node(198))
    FIX_Z.append(Mesh_Final.node(199))
    FIX_Z.append(Mesh_Final.node(200))
    FIX_Z.append(Mesh_Final.node(201))
    FIX_Z.append(Mesh_Final.node(202))
    FIX_Z.append(Mesh_Final.node(203))
    FIX_Z.append(Mesh_Final.node(204))
    FIX_Z.append(Mesh_Final.node(205))
    FIX_Z.append(Mesh_Final.node(206))
    FIX_Z.append(Mesh_Final.node(207))
    FIX_Z.append(Mesh_Final.node(208))
    FIX_Z.append(Mesh_Final.node(209))
    FIX_Z.append(Mesh_Final.node(210))
    FIX_Z.append(Mesh_Final.node(211))
    FIX_Z.append(Mesh_Final.node(212))
    FIX_Z.append(Mesh_Final.node(213))
    FIX_Z.append(Mesh_Final.node(214))
    FIX_Z.append(Mesh_Final.node(215))
    FIX_Z.append(Mesh_Final.node(216))
    FIX_Z.append(Mesh_Final.node(217))
    FIX_Z.append(Mesh_Final.node(218))
    FIX_Z.append(Mesh_Final.node(219))
    FIX_Z.append(Mesh_Final.node(220))
    FIX_Z.append(Mesh_Final.node(221))
    FIX_Z.append(Mesh_Final.node(222))
    FIX_Z.append(Mesh_Final.node(223))
    FIX_Z.append(Mesh_Final.node(224))
    FIX_Z.append(Mesh_Final.node(225))
    FIX_Z.append(Mesh_Final.node(226))
    FIX_Z.append(Mesh_Final.node(227))
    FIX_Z.append(Mesh_Final.node(228))
    FIX_Z.append(Mesh_Final.node(229))
    FIX_Z.append(Mesh_Final.node(230))

    N66 = Mesh_Final.node(65)
    N77 = Mesh_Final.node(76)
    N88 = Mesh_Final.node(87)
    N99 = Mesh_Final.node(98)
    N110 = Mesh_Final.node(109)
    N121 = Mesh_Final.node(120)
    N187 = Mesh_Final.node(186)
    N242 = Mesh_Final.node(241)
    N198 = Mesh_Final.node(197)
    N209 = Mesh_Final.node(208)
    N220 = Mesh_Final.node(219)
    N231 = Mesh_Final.node(230)

    Point_A=Mesh_Final.node(110)
    Point_B=Mesh_Final.node(55)
    Point_C=Mesh_Final.node(0)

    No_Trait=[]  
    No_Trait.append(Mesh_Final.node(0  ))
    No_Trait.append(Mesh_Final.node(55 ))
    No_Trait.append(Mesh_Final.node(110))
    
     
    bulk = BULK_MATERIAL (solfec,
                          model = 'KIRCHHOFF',
                          young = 3.0E+10,
                          poisson = 0.25,
                          density = 7.82E+03)
                          
    Bod_M1 = BODY (solfec, 'FINITE_ELEMENT',Mesh_Final,bulk)

    for i in range (0,len(Node)):
      DISPLAY_POINT (Bod_M1, Node[i], 'N%d'%i) # press 'D' in the viewer to see display points

    Fix_direc_X=(1.,0.,0.)
    Fix_direc_Y=(0.,1.,0.)
    Fix_direc_Z=(0.,0.,1.)

    for fix in FIX_Y:
       FIX_DIRECTION (Bod_M1, fix, Fix_direc_Y)

    for fix in FIX_X:
      FIX_DIRECTION (Bod_M1, fix, Fix_direc_X)

    for fix in FIX_Z:
      FIX_DIRECTION (Bod_M1, fix, Fix_direc_Z)


    F1=12.5
    F2=25

    FORCE(Bod_M1,'SPATIAL', N77,(0.,1.,0.),F2)
    FORCE(Bod_M1,'SPATIAL', N88,(0.,1.,0.),F2)
    FORCE(Bod_M1,'SPATIAL', N99,(0.,1.,0.),F2)
    FORCE(Bod_M1,'SPATIAL', N110,(0.,1.,0.),F2)
    FORCE(Bod_M1,'SPATIAL', N198,(0.,1.,0.),F2)
    FORCE(Bod_M1,'SPATIAL', N209,(0.,1.,0.),F2)
    FORCE(Bod_M1,'SPATIAL', N220,(0.,1.,0.),F2)
    FORCE(Bod_M1,'SPATIAL', N231,(0.,1.,0.),F2)

    FORCE(Bod_M1,'SPATIAL', N66,(0.,1.,0.),F1)
    FORCE(Bod_M1,'SPATIAL', N121,(0.,1.,0.),F1)
    FORCE(Bod_M1,'SPATIAL', N187,(0.,1.,0.),F1)
    FORCE(Bod_M1,'SPATIAL', N242,(0.,1.,0.),F1)
  
    nt = NEWTON_SOLVER (1E-9)
    OUTPUT(solfec, outstep)
    if not VIEWER() and solfec.mode == 'WRITE':
      solfec.verbose = '%'
    RUN (solfec, nt, stop)
    return (solfec, Point_A, Point_B, Point_C)
      
# data
step = 2E-4
stop = 0.004
outstep = 0.001

# calculate
(solfec, A, B, C) = create_simulation (step, stop, outstep)

# read results
if not VIEWER():
  (solfec, A, B, C) = create_simulation (step, stop, outstep)
  if solfec.mode != 'READ':
    print 'READ ERROR'
    import sys
    sys.exit(1)

  body = solfec.bodies[0]
  mesh0 = COPY(body.mesh)

  failed = False
  SEEK (solfec, stop)
  label = ('SX', 'SY', 'SZ')
  points = (C, B, A)
  #      ((5.26E+05, 7.46E+06, 9.22E+05), (5.27E+05, 2.40E+06, 7.31E+05), (2.31E+06, -1.33E+05, 2.95E+05)) # Code_Aster
  vref = ((5.26E+05, 7.46E+06, -9.22E+05), (5.27E+05, 2.40E+06, 7.31E+05), (-2.31E+06, -1.33E+05, 2.95E+05)) # XXX: modified signs
  vtol = ((0.13, 0.035, 0.76), (0.36, 0.096, 0.84), (0.0064, 0.43, 0.83)) # Solfec c060f54 tolerance
  for (point, refe, tole) in zip(points, vref, vtol):
    stre = STRESS (body, point)
    for i in range (0,3):
      Value = stre[i]
      Code_Aster_Ref = refe[i]
      error = abs(Value-Code_Aster_Ref)/abs(Code_Aster_Ref) 
      tolerance = tole[i]

      if error > tolerance:
	failed = True
	print '\b\b\b\bFAILED', 
	print '(Computed stress %s at point (%g,%g,%g) was %g' % (label[i], point[0], point[1], point[2], Value),
	print 'while reference is %g)' % Code_Aster_Ref, '-->(error %g, while tolerance %g)' % (error, tolerance)

  if not failed: print '\b\b\b\bPASSED'
