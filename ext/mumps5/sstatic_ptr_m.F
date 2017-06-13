C
C  This file is part of MUMPS 5.1.1, released
C  on Mon Mar 20 14:34:33 UTC 2017
C
C
C  Copyright 1991-2017 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license:
C  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
C
      MODULE SMUMPS_STATIC_PTR_M
      PUBLIC :: SMUMPS_TMP_PTR, SMUMPS_GET_TMP_PTR
      REAL, DIMENSION(:), POINTER, SAVE :: SMUMPS_TMP_PTR
      CONTAINS
      SUBROUTINE SMUMPS_SET_STATIC_PTR(ARRAY)
      REAL, DIMENSION(:), TARGET :: ARRAY
      SMUMPS_TMP_PTR => ARRAY
      RETURN
      END SUBROUTINE SMUMPS_SET_STATIC_PTR
      SUBROUTINE SMUMPS_GET_TMP_PTR(PTR)
#if defined(MUMPS_F2003)
      REAL, DIMENSION(:), POINTER, INTENT(OUT) :: PTR
#else
      REAL, DIMENSION(:), POINTER :: PTR
#endif
      PTR => SMUMPS_TMP_PTR
      RETURN
      END SUBROUTINE SMUMPS_GET_TMP_PTR
      END MODULE SMUMPS_STATIC_PTR_M
