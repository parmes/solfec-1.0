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
#include "mumps_io.h"
#include "mumps_io_basic.h"
#include "mumps_io_err.h"
#if ! defined (MUMPS_WIN32) && ! defined (WITHOUT_PTHREAD)
# include "mumps_io_thread.h"
#endif
#if ! defined(MUMPS_WIN32)
double mumps_time_spent_in_sync;
#endif
double read_op_vol,write_op_vol,total_vol;
/**
 * Forward declaration. Definition at the end of the file.
 */
MUMPS_INLINE int
mumps_convert_2fint_to_longlong( int *, int *, long long *);
/* Tests if the request "request_id" has finished. It sets the flag  */
/* argument to 1 if the request has finished (0 otherwise)           */
void MUMPS_CALL
MUMPS_TEST_REQUEST_C(int* request_id,int *flag,int* ierr)
{
  char buf[64]; /* for error message */
#if ! defined(MUMPS_WIN32)
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  switch(mumps_io_flag_async){
  case IO_SYNC:
    /* printf("mumps_test_request_c should not be called with strategy %d\n",mumps_io_flag_async);*/
    /* JY+EA: Allow for this option, since it is similar to wait_request
     * and wait_request is allowed in synchronous mode.
     * We always return TRUE.
     */
    *flag=1;
    break;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  case IO_ASYNC_TH:
    *ierr=mumps_test_request_th(request_id,flag);
    break;
#endif
  default:
    *ierr=-92;
    sprintf(buf,"Error: unknown I/O strategy : %d\n",mumps_io_flag_async);
    mumps_io_error(*ierr,buf);
    return;
  }
#if ! defined(MUMPS_WIN32)
  gettimeofday(&end_time,NULL);
  mumps_time_spent_in_sync=mumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  return;
}
/* Waits for the termination of the request "request_id" */
void MUMPS_CALL
MUMPS_WAIT_REQUEST(int *request_id,int* ierr)
{
  char buf[64]; /* for error message */
#if ! defined(MUMPS_WIN32)
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  if(*request_id==-1)
    return;
  switch(mumps_io_flag_async){
  case IO_SYNC:
    /* printf("mumps_wait_request should not be called with strategy %d\n",mumps_io_flag_async); */
    break;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  case IO_ASYNC_TH:
    *ierr=mumps_wait_request_th(request_id);
    break;
#endif
  default:
    *ierr=-92;
    sprintf(buf,"Error: unknown I/O strategy : %d\n",mumps_io_flag_async);
    mumps_io_error(*ierr,buf);
    return;
    /*    printf("Error: unknown I/O strategy : %d\n",mumps_io_flag_async);
          exit (-3);*/
  }
#if ! defined(MUMPS_WIN32)
  gettimeofday(&end_time,NULL);
  mumps_time_spent_in_sync=mumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  return;
}
/**
 * Inits the I/O OOC mechanism.
 * Because on some computers, file size is limited, the I/O
 * mechanism must be able to handle a multi-file access to data.
 * Hence, we compute mumps_io_nb_file, which is the the number of files
 * we estimate we need.
 * Because of not exact matching between data packets written and size
 * of files, the recoverment may be imperfect. Consequently, we must
 * be able to reallocate if necessary.
 */
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_PREFIX(int * dim, char * str, mumps_ftnlen l1)
{
  int i;
  MUMPS_OOC_STORE_PREFIXLEN = *dim;
  if( *dim > MUMPS_OOC_PREFIX_MAX_LENGTH )
      MUMPS_OOC_STORE_PREFIXLEN = MUMPS_OOC_PREFIX_MAX_LENGTH;
  for(i=0;i<MUMPS_OOC_STORE_PREFIXLEN;i++){
    MUMPS_OOC_STORE_PREFIX[i]=str[i];
  }
  return;
}
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_TMPDIR(int * dim, char * str, mumps_ftnlen l1)
{
  int i;
  MUMPS_OOC_STORE_TMPDIRLEN=*dim;
  if( *dim > MUMPS_OOC_TMPDIR_MAX_LENGTH )
      MUMPS_OOC_STORE_TMPDIRLEN = MUMPS_OOC_TMPDIR_MAX_LENGTH;
  for(i=0;i<MUMPS_OOC_STORE_TMPDIRLEN;i++){
    MUMPS_OOC_STORE_TMPDIR[i]=str[i];
  }
  return;
}
/* Computes the number of files needed. Uses ceil value. */
/*   mumps_io_nb_file=0; */
/*   mumps_io_last_file_opened=-1; */
void MUMPS_CALL
MUMPS_LOW_LEVEL_INIT_OOC_C(int* _myid, int* total_size_io, int* size_element,
                           int* async, int* k211, int * nb_file_type,
                           int * flag_tab, int* ierr)
{
  char buf[64]; /* for error message */
#if defined(MUMPS_WIN32)
  if(*async==IO_ASYNC_AIO||*async==IO_ASYNC_TH){
    mumps_io_is_init_called=0;
    *ierr=-92;
    mumps_io_error(*ierr,"Error: Forbidden value of Async flag with _WIN32\n");
    return;
  }
#endif
#if defined (WITHOUT_PTHREAD)
  if(*async==IO_ASYNC_TH){
    mumps_io_is_init_called=0;
    *ierr=-92;
    mumps_io_error(*ierr,"Error: Forbidden value of Async flag with WITHOUT_PTHREAD\n");
    return;
  }
#endif
  total_vol=0;
  mumps_io_flag_async=*async;
  mumps_io_k211=*k211;
  if (MUMPS_OOC_STORE_PREFIXLEN==-1) {
    *ierr=-92;
    mumps_io_error(*ierr,"Error: prefix not initialized\n");
    return;
  }
  if (MUMPS_OOC_STORE_TMPDIRLEN==-1) {
    *ierr=-92;
    mumps_io_error(*ierr,"Error: tmpdir not initialized\n");
    return;
  }
  *ierr=mumps_init_file_name(MUMPS_OOC_STORE_TMPDIR, MUMPS_OOC_STORE_PREFIX,
		             &MUMPS_OOC_STORE_TMPDIRLEN, &MUMPS_OOC_STORE_PREFIXLEN, _myid);
  if(*ierr<0){
    return;
  }
  /* Re-initialize lenghts to -1 in order to enable the
   * check on initialization next time this routine is called
   */
  MUMPS_OOC_STORE_PREFIXLEN=-1;
  MUMPS_OOC_STORE_TMPDIRLEN=-1;
  *ierr=mumps_init_file_structure(_myid,total_size_io,size_element,*nb_file_type,flag_tab);
  if(*ierr<0){
    return;
  }
#if ! defined(MUMPS_WIN32)
  mumps_time_spent_in_sync=0;
#endif
  if(*async){
    switch(*async){
    case IO_SYNC:
      printf("mumps_low_level_init_ooc_c should not be called with strategy %d\n",mumps_io_flag_async);
      break;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
    case IO_ASYNC_TH:
      mumps_low_level_init_ooc_c_th(async,ierr);
      if(*ierr<0){
        return;
      }
      break;
#endif
    default:
      *ierr=-92;
      sprintf(buf,"Error: unknown I/O strategy : %d\n",*async);
      mumps_io_error(*ierr,buf);
      return;
    }
  }
  mumps_io_is_init_called=1;
  return;
}
/**
 * Writes a contigous block of central memory to the disk.
 */
void MUMPS_CALL
MUMPS_LOW_LEVEL_WRITE_OOC_C(const int * strat_IO,
                            void * address_block,
                            int * block_size_int1,
                            int * block_size_int2,
                            int * inode,
                            int * request_arg,
                            int * type,
                            int * vaddr_int1,
                            int * vaddr_int2,
                            int * ierr)
{
  int ret_code=0;
  long long vaddr,block_size;
  char buf[64]; /* for error message */
#if ! defined(MUMPS_WIN32)
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
/* JY 27/2/08: initialize *request_arg to -1 (null request).
 * There were problems of uninitialized requests in the Fortran
 * code. For example when we use the synchronous version, there are
 * still some tests on *request_arg, which is not initialized.*/
  *request_arg=-1;
  mumps_convert_2fint_to_longlong(vaddr_int1,vaddr_int2,&vaddr);
  mumps_convert_2fint_to_longlong(block_size_int1,block_size_int2,&block_size);
  if(mumps_io_flag_async){
    switch(*strat_IO){
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
    case IO_ASYNC_TH:
      ret_code=mumps_async_write_th(strat_IO, address_block, block_size,inode,request_arg,type,vaddr,ierr);
      if(ret_code<0){
        *ierr=ret_code;
      }
      break;
#endif
    default:
      *ierr=-91;
      sprintf(buf,"Error: unknown I/O strategy : %d\n",*strat_IO);
      mumps_io_error(*ierr,buf);
      return;
    }
  } else {
    ret_code=mumps_io_do_write_block(address_block,block_size,type,vaddr,ierr);
    if(ret_code<0){
      *ierr=ret_code;
    }
  }
#if ! defined(MUMPS_WIN32)
  gettimeofday(&end_time,NULL);
  mumps_time_spent_in_sync=mumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  write_op_vol=write_op_vol+((double)(block_size)*(double)mumps_elementary_data_size);
  return;
}
/**
 * Reads  a contigous block of central memory from the disk.
 */
void MUMPS_CALL
MUMPS_LOW_LEVEL_READ_OOC_C(const int * strat_IO,
                           void * address_block,
                           int * block_size_int1,
                           int * block_size_int2,
                           int * inode,
                           int * request_arg,
                           int * type,
                           int * vaddr_int1,
                           int * vaddr_int2,
                           int * ierr)
{
  char buf[64]; /* for error message */
  long long vaddr,block_size;
#if ! defined(MUMPS_WIN32)
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  mumps_convert_2fint_to_longlong(vaddr_int1,vaddr_int2,&vaddr);
  mumps_convert_2fint_to_longlong(block_size_int1,block_size_int2,&block_size);
  if(mumps_io_flag_async){
      switch(*strat_IO){
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
      case IO_ASYNC_TH:
        mumps_async_read_th(strat_IO,address_block,block_size,inode,request_arg,type,vaddr,ierr);
        break;
#endif
      default:
        *ierr=-91;
        sprintf(buf,"Error: unknown I/O strategy : %d\n",*strat_IO);
        mumps_io_error(*ierr,buf);
        return;
      }
  }else{
    mumps_io_do_read_block(address_block,block_size,type,vaddr,ierr);
    *request_arg=1;
  }
#if ! defined(MUMPS_WIN32)
  gettimeofday(&end_time,NULL);
  mumps_time_spent_in_sync=mumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  read_op_vol=read_op_vol+((double)(block_size)*(double)mumps_elementary_data_size);
  return;
}
/* Emergency read from the MUMPS thread during the solve phase.*/
void MUMPS_CALL
MUMPS_LOW_LEVEL_DIRECT_READ(void * address_block,
                            int * block_size_int1,
                            int * block_size_int2,
                            int * type,
                            int * vaddr_int1,
                            int * vaddr_int2,
                            int * ierr)
{
    /*  int ret_code=0; */
  long long vaddr,block_size;
#if ! defined(MUMPS_WIN32)
  struct timeval start_time,end_time;
  gettimeofday(&start_time,NULL);
#endif
  mumps_convert_2fint_to_longlong(vaddr_int1,vaddr_int2,&vaddr);
  mumps_convert_2fint_to_longlong(block_size_int1,block_size_int2,&block_size);
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
    if(mumps_io_flag_async==IO_ASYNC_TH||mumps_io_flag_async==0)
#else
    if(mumps_io_flag_async==0)
#endif
    {
      *ierr=mumps_io_do_read_block(address_block,block_size,type,vaddr,ierr);
      if(*ierr<0){
	 return;
      }
    } else {
    }
#if ! defined(MUMPS_WIN32)
# if ! defined(WITHOUT_PTHREAD)
# endif
  gettimeofday(&end_time,NULL);
  mumps_time_spent_in_sync=mumps_time_spent_in_sync+((double)end_time.tv_sec+((double)end_time.tv_usec/1000000))-((double)start_time.tv_sec+((double)start_time.tv_usec/1000000));
#endif
  read_op_vol=read_op_vol+((double)(block_size)*(double)mumps_elementary_data_size);
  return;
}
/* Cleans the thread/io management data*/
void MUMPS_CALL
MUMPS_CLEAN_IO_DATA_C(int * myid,int *step,int *ierr)
{
  char buf[64]; /* for error message */
  if(!mumps_io_is_init_called){
    return;
  }
  switch(mumps_io_flag_async){
  case IO_SYNC:
    break;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  case IO_ASYNC_TH:
    *ierr=mumps_clean_io_data_c_th(myid);
    break;
#endif
  default:
    *ierr=-91;
    sprintf(buf,"Error: unknown I/O strategy : %d\n",mumps_io_flag_async);
    mumps_io_error(*ierr,buf);
    return;
  }
  mumps_free_file_pointers(step);
  mumps_io_is_init_called=0;
  return;
}
void MUMPS_CALL
MUMPS_OOC_PRINT_STATS()
{
#if ! defined(MUMPS_WIN32)
  printf("%d: total time spent in i/o mode = %lf\n",mumps_io_myid,mumps_time_spent_in_sync);
#endif
  printf("%d: Volume of read i/o = %lf\n",mumps_io_myid,read_op_vol);
  printf("%d: Volume of write i/o = %lf\n",mumps_io_myid,write_op_vol);
  total_vol=total_vol+read_op_vol+write_op_vol;
  printf("%d: Total i/o volume = %lf\n",mumps_io_myid,total_vol);
  return;
}
void MUMPS_CALL
MUMPS_GET_MAX_NB_REQ_C(int *max,int *ierr)
{
  char buf[64]; /* for error message */
  *ierr=0;
  switch(mumps_io_flag_async){
  case IO_SYNC:
    *max=1;
    break;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  case IO_ASYNC_TH:
    *max=MAX_FINISH_REQ+MAX_IO;
    break;
#endif
  default:
    *ierr=-91;
    sprintf(buf,"Error: unknown I/O strategy : %d\n",mumps_io_flag_async);
    mumps_io_error(*ierr,buf);
    return;
  }
  return;
}
void MUMPS_CALL
MUMPS_GET_MAX_FILE_SIZE_C(double * max_ooc_file_size)
{
  *max_ooc_file_size=(double)(MAX_FILE_SIZE);
  return;
}
void MUMPS_CALL
MUMPS_OOC_GET_NB_FILES_C(const int* type,int* nb_files)
{
  mumps_io_get_nb_files(nb_files,type);
  return;
}
void MUMPS_CALL
MUMPS_OOC_GET_FILE_NAME_C(int* type,int* indice,int* length, char* name, mumps_ftnlen l1)
{
  mumps_io_get_file_name(indice,name,length,type);
  return;
}
void MUMPS_CALL
MUMPS_OOC_SET_FILE_NAME_C(int* type, int* indice, int* length, int *ierr,
                          char* name, mumps_ftnlen l1)
{
  *ierr=mumps_io_set_file_name(indice,name,length,type);
  return;
}
void MUMPS_CALL
MUMPS_OOC_ALLOC_POINTERS_C(int* nb_file_type,int* dim,int* ierr)
{
  int i=0;
  *ierr=mumps_io_alloc_pointers(nb_file_type,dim);
  for(i=0;i<*nb_file_type;i++){
    mumps_io_set_last_file((dim+i),&i);
  }
  return;
}
void MUMPS_CALL
MUMPS_OOC_INIT_VARS_C(int* myid_arg,
                        int* size_element,int* async, int* k211,
                        int *ierr)
{
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
  mumps_time_spent_in_sync=0;
#endif
  mumps_io_k211=*k211;
  *ierr=mumps_io_init_vars(myid_arg,size_element,async);
  return;
}
void MUMPS_CALL
MUMPS_OOC_START_LOW_LEVEL(int *ierr)
{
  char buf[64]; /* for error message */
  read_op_vol=0;
  write_op_vol=0;
  *ierr=mumps_io_open_files_for_read();
  if(*ierr<0){
    return;
  }
  if(mumps_io_flag_async){
    switch(mumps_io_flag_async){
    case IO_SYNC:
      break;
#if ! defined(MUMPS_WIN32) && ! defined(WITHOUT_PTHREAD)
    case IO_ASYNC_TH:
      mumps_low_level_init_ooc_c_th(&mumps_io_flag_async,ierr);
      if(*ierr<0){
        return;
      }
      break;
#endif
    default:
      *ierr=-91;
      sprintf(buf,"Error: unknown I/O strategy : %d\n",mumps_io_flag_async);
      mumps_io_error(*ierr,buf);
      return;
    }
  }
  mumps_io_is_init_called=1;
  return;
}
void MUMPS_CALL
MUMPS_OOC_REMOVE_FILE_C(int *ierr, char *name, mumps_ftnlen l1)
{
  char buf[296]; /* for error message, count 256 chars for name */
  *ierr=remove(name);
  if(*ierr<0){
#if ! defined(MUMPS_WIN32)
    sprintf(buf,"Unable to remove OOC file %s",name);
#else
    sprintf(buf,"Unable to remove OOC file %s with return value %d",name,*ierr);
#endif
    *ierr = -90;
    mumps_io_sys_error(*ierr,buf);
    return;
  }
  return;
}
void MUMPS_CALL
MUMPS_OOC_END_WRITE_C(int *ierr)
{
  return;
}
void MUMPS_CALL
MUMPS_OOC_IS_ASYNC_AVAIL(int *flag)
{
#if ( defined (WITHOUT_PTHREAD) || defined(MUMPS_WIN32) ) && ! defined(WITH_AIO)
  *flag=0;
#else
  *flag=1;
#endif
}
/**
 * IMPORTANT:
 * ==========
 *
 *   After every modification of the code of convert_2fint_to_longlong update
 *   the corresponding fortran subroutines MUMPS_OOC_CONVERT_2INTTOVADDR
 *   and MUMPS_OOC_CONVERT_VADDRTO2INT
 */
MUMPS_INLINE int
mumps_convert_2fint_to_longlong( int * short_int1, int * short_int2,
                                 long long * long_int )
{
  *long_int=((long long)(*short_int1)*((long long)1073741824))+(long long)(*short_int2);
  return 0;
}
