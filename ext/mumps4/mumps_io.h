/*
 *
 *  This file is part of MUMPS 4.9.2, built on Thu Nov  5 07:05:08 UTC 2009
 *
 *
 *  This version of MUMPS is provided to you free of charge. It is public
 *  domain, based on public domain software developed during the Esprit IV
 *  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL.
 *  Since this first public domain version in 1999, the developments are
 *  supported by the following institutions: CERFACS, CNRS, INPT(ENSEEIHT)-
 *  IRIT, and INRIA.
 *
 *  Current development team includes Patrick Amestoy, Alfredo Buttari,
 *  Abdou Guermouche, Jean-Yves L'Excellent, Bora Ucar.
 *
 *  Up-to-date copies of the MUMPS package can be obtained
 *  from the Web pages:
 *  http://mumps.enseeiht.fr/  or  http://graal.ens-lyon.fr/MUMPS
 *
 *
 *   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
 *   EXPRESSED OR IMPLIED. ANY USE IS AT YOUR OWN RISK.
 *
 *
 *  User documentation of any code that uses this software can
 *  include this complete notice. You can acknowledge (using
 *  references [1] and [2]) the contribution of this package
 *  in any scientific publication dependent upon the use of the
 *  package. You shall use reasonable endeavours to notify
 *  the authors of the package of this publication.
 *
 *   [1] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
 *   A fully asynchronous multifrontal solver using distributed dynamic
 *   scheduling, SIAM Journal of Matrix Analysis and Applications,
 *   Vol 23, No 1, pp 15-41 (2001).
 *
 *   [2] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
 *   S. Pralet, Hybrid scheduling for the parallel solution of linear
 *   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
 *
 */
#ifndef MUMPS_IO_H
#define MUMPS_IO_H
#include "mumps_common.h"
#include "mumps_c_types.h"
/*
 *  Two character arrays that are used by low_level_init_prefix
 *  and low_level_init_tmpdir to store intermediate file names.
 *  They are passed to mumps_io_basic.c inside the routine
 *  mumps_low_level_init_ooc_c.
 *  Note that both low_level_init_prefix and low_level_init_tmpdir
 *  MUST be called before low_level_init_ooc_c.
 * 
 */
#define MUMPS_OOC_PREFIX_MAX_LENGTH 63
#define MUMPS_OOC_TMPDIR_MAX_LENGTH 255
static char MUMPS_OOC_STORE_PREFIX[MUMPS_OOC_PREFIX_MAX_LENGTH];
static int  MUMPS_OOC_STORE_PREFIXLEN=-1;
static char MUMPS_OOC_STORE_TMPDIR[MUMPS_OOC_TMPDIR_MAX_LENGTH];
static int  MUMPS_OOC_STORE_TMPDIRLEN=-1;
#define MUMPS_LOW_LEVEL_INIT_PREFIX \
    F_SYMBOL(low_level_init_prefix,LOW_LEVEL_INIT_PREFIX)
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_PREFIX(int * dim, char * str, mumps_ftnlen l1);
#define MUMPS_LOW_LEVEL_INIT_TMPDIR \
    F_SYMBOL(low_level_init_tmpdir,LOW_LEVEL_INIT_TMPDIR)
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_TMPDIR(int * dim, char * str, mumps_ftnlen l1);
#define MUMPS_LOW_LEVEL_INIT_OOC_C \
    F_SYMBOL(low_level_init_ooc_c,LOW_LEVEL_INIT_OOC_C)
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_OOC_C(int* _myid, int* total_size_io,int* size_element,
                           int* async, int* k211, int * nb_file_type,
                           int * flag_tab , int* ierr);
#define MUMPS_TEST_REQUEST_C \
    F_SYMBOL(test_request_c,TEST_REQUEST_C)
void MUMPS_CALL
MUMPS_TEST_REQUEST_C(int* request_id,int* flag,int* ierr);
#define MUMPS_WAIT_REQUEST \
    F_SYMBOL(wait_request,WAIT_REQUEST)
void MUMPS_CALL
MUMPS_WAIT_REQUEST(int* request_id,int* ierr);
#define MUMPS_LOW_LEVEL_WRITE_OOC_C \
    F_SYMBOL(low_level_write_ooc_c,LOW_LEVEL_WRITE_OOC_C)
void MUMPS_CALL
MUMPS_LOW_LEVEL_WRITE_OOC_C( const int * strat_IO, 
                             void * address_block,
                             int * block_size_int1,
                             int * block_size_int2,
                             int * inode,
                             int * request_arg,
                             int * type,
                             int * vaddr_int1,
                             int * vaddr_int2,
                             int * ierr);
#define MUMPS_LOW_LEVEL_READ_OOC_C \
    F_SYMBOL(low_level_read_ooc_c,LOW_LEVEL_READ_OOC_C)
void MUMPS_CALL
MUMPS_LOW_LEVEL_READ_OOC_C( const int * strat_IO, 
                            void * address_block,
                            int * block_size_int1,
                            int * block_size_int2,
                            int * inode,
                            int * request_arg,
                            int * type,
                            int * vaddr_int1,
                            int * vaddr_int2,
                            int * ierr);
#define MUMPS_LOW_LEVEL_DIRECT_READ \
    F_SYMBOL(low_level_direct_read,LOW_LEVEL_DIRECT_READ)
void MUMPS_CALL
MUMPS_LOW_LEVEL_DIRECT_READ(void * address_block,
                            int * block_size_int1,
                            int * block_size_int2,
                            int * type,
                            int * vaddr_int1,
                            int * vaddr_int2,
                            int * ierr);
#define MUMPS_CLEAN_IO_DATA_C \
    F_SYMBOL(clean_io_data_c,CLEAN_IO_DATA_C)
void MUMPS_CALL
MUMPS_CLEAN_IO_DATA_C(int* myid,int* step,int* ierr);
#define MUMPS_GET_MAX_NB_REQ_C \
    F_SYMBOL(get_max_nb_req_c,GET_MAX_NB_REQ_C)
void MUMPS_CALL
MUMPS_GET_MAX_NB_REQ_C(int *max,int* ierr);
#define MUMPS_GET_MAX_FILE_SIZE_C \
    F_SYMBOL(get_max_file_size_c,GET_MAX_FILE_SIZE_C)
void MUMPS_CALL
MUMPS_GET_MAX_FILE_SIZE_C(double * max_ooc_file_size);
#define MUMPS_OOC_GET_NB_FILES_C \
    F_SYMBOL(ooc_get_nb_files_c,OOC_GET_NB_FILES_C)
void MUMPS_CALL
MUMPS_OOC_GET_NB_FILES_C(const int* type, int* nb_files);
#define MUMPS_OOC_GET_FILE_NAME_C \
    F_SYMBOL(ooc_get_file_name_c,OOC_GET_FILE_NAME_C)
void MUMPS_CALL
MUMPS_OOC_GET_FILE_NAME_C(int* type, int* indice, int * length,
                          char* name, mumps_ftnlen l1);
#define MUMPS_OOC_SET_FILE_NAME_C \
    F_SYMBOL(ooc_set_file_name_c,OOC_SET_FILE_NAME_C)
void MUMPS_CALL
MUMPS_OOC_SET_FILE_NAME_C(int* type, int* indice, int* length, int* ierr,
                          char* name, mumps_ftnlen l1);
#define MUMPS_OOC_ALLOC_POINTERS_C \
    F_SYMBOL(ooc_alloc_pointers_c,OOC_ALLOC_POINTERS_C)
void MUMPS_CALL
MUMPS_OOC_ALLOC_POINTERS_C(int* nb_file_type, int* dim, int* ierr);
#define MUMPS_OOC_INIT_VARS_C \
    F_SYMBOL(ooc_init_vars_c,OOC_INIT_VARS_C)
void MUMPS_CALL
MUMPS_OOC_INIT_VARS_C(int* myid_arg, int* size_element, int* async,
                      int* k211, int *ierr);
#define MUMPS_OOC_START_LOW_LEVEL \
    F_SYMBOL(ooc_start_low_level,OOC_START_LOW_LEVEL)
void MUMPS_CALL
MUMPS_OOC_START_LOW_LEVEL(int* ierr);
#define MUMPS_OOC_PRINT_STATS \
    F_SYMBOL(ooc_print_stats,OOC_PRINT_STATS)
void MUMPS_CALL
MUMPS_OOC_PRINT_STATS();
#define MUMPS_OOC_REMOVE_FILE_C \
    F_SYMBOL(ooc_remove_file_c,OOC_REMOVE_FILE_C)
void MUMPS_CALL
MUMPS_OOC_REMOVE_FILE_C(int *ierr, char *name, mumps_ftnlen l1);
#define MUMPS_OOC_END_WRITE_C \
    F_SYMBOL(ooc_end_write_c,OOC_END_WRITE_C)
void MUMPS_CALL
MUMPS_OOC_END_WRITE_C(int *ierr);
#define MUMPS_OOC_IS_ASYNC_AVAIL \
    F_SYMBOL(ooc_is_async_avail,OOC_IS_ASYNC_AVAIL)
void MUMPS_CALL
MUMPS_OOC_IS_ASYNC_AVAIL(int *flag);
#endif /* MUMPS_IO_H */
