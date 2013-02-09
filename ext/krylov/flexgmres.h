/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision: 2.3 $
 ***********************************************************************EHEADER*/





/******************************************************************************
 *
 * FLEXGMRES flexible gmres
 *
 *****************************************************************************/

#ifndef hypre_KRYLOV_FLEXGMRES_HEADER
#define hypre_KRYLOV_FLEXGMRES_HEADER

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

/**
 * @name Generic FlexGMRES Interface
 *
 * A general description of the interface goes here...
 *
 **/
/*@{*/

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
 * hypre_FlexGMRESData and hypre_FlexGMRESFunctions
 *--------------------------------------------------------------------------*/


/**
 * @name FlexGMRES structs
 *
 * Description...
 **/
/*@{*/

/**
 * The {\tt hypre\_FlexGMRESFunctions} object ...
 **/

typedef struct
{
   char * (*CAlloc)        ( size_t count, size_t elt_size );
   int    (*Free)          ( char *ptr );
   int    (*CommInfo)      ( void  *A, int   *my_id, int   *num_procs );
   void * (*CreateVector)  ( void *vector );
   void * (*CreateVectorArray)  ( int size, void *vectors );
   int    (*DestroyVector) ( void *vector );
   void * (*MatvecCreate)  ( void *A, void *x );
   int    (*Matvec)        ( void *matvec_data, double alpha, void *A,
                             void *x, double beta, void *y );
   int    (*MatvecDestroy) ( void *matvec_data );
   double (*InnerProd)     ( void *x, void *y );
   int    (*CopyVector)    ( void *x, void *y );
   int    (*ClearVector)   ( void *x );
   int    (*ScaleVector)   ( double alpha, void *x );
   int    (*Axpy)          ( double alpha, void *x, void *y );

   int    (*precond)();
   int    (*precond_setup)();

   int    (*modify_pc)();
   

} hypre_FlexGMRESFunctions;

/**
 * The {\tt hypre\_FlexGMRESData} object ...
 **/



typedef struct
{
   int      k_dim;
   int      min_iter;
   int      max_iter;
   int      rel_change;
   int      stop_crit;
   int      converged;
   double   tol;
   double   cf_tol;
   double   a_tol;
   double   rel_residual_norm;

   void   **pre_vecs;

   void  *A;
   void  *r;
   void  *w;
   void  *w_2;
   void  **p;

   void    *matvec_data;
   void    *precond_data;

   hypre_FlexGMRESFunctions * functions;

   /* log info (always logged) */
   int      num_iterations;
 
   int     print_level; /* printing when print_level>0 */
   int     logging;  /* extra computations for logging when logging>0 */
   double  *norms;
   char    *log_file_name;

} hypre_FlexGMRESData;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @name generic FlexGMRES Solver
 *
 * Description...
 **/
/*@{*/

/**
 * Description...
 *
 * @param param [IN] ...
 **/

hypre_FlexGMRESFunctions *
hypre_FlexGMRESFunctionsCreate(
   char * (*CAlloc)        ( size_t count, size_t elt_size ),
   int    (*Free)          ( char *ptr ),
   int    (*CommInfo)      ( void  *A, int   *my_id, int   *num_procs ),
   void * (*CreateVector)  ( void *vector ),
   void * (*CreateVectorArray)  ( int size, void *vectors ),
   int    (*DestroyVector) ( void *vector ),
   void * (*MatvecCreate)  ( void *A, void *x ),
   int    (*Matvec)        ( void *matvec_data, double alpha, void *A,
                             void *x, double beta, void *y ),
   int    (*MatvecDestroy) ( void *matvec_data ),
   double (*InnerProd)     ( void *x, void *y ),
   int    (*CopyVector)    ( void *x, void *y ),
   int    (*ClearVector)   ( void *x ),
   int    (*ScaleVector)   ( double alpha, void *x ),
   int    (*Axpy)          ( double alpha, void *x, void *y ),
   int    (*PrecondSetup)  ( void *vdata, void *A, void *b, void *x ),
   int    (*Precond)       ( void *vdata, void *A, void *b, void *x )
   );

/**
 * Description...
 *
 * @param param [IN] ...
 **/

void *
hypre_FlexGMRESCreate( hypre_FlexGMRESFunctions *fgmres_functions );
int hypre_FlexGMRESDestroy ( void *fgmres_vdata );
int hypre_FlexGMRESGetResidual ( void *fgmres_vdata , void **residual );
int hypre_FlexGMRESSetup ( void *fgmres_vdata , void *A , void *b , void *x );
int hypre_FlexGMRESSolve ( void *fgmres_vdata , void *A , void *b , void *x );
int hypre_FlexGMRESSetKDim ( void *fgmres_vdata , int k_dim );
int hypre_FlexGMRESGetKDim ( void *fgmres_vdata , int *k_dim );
int hypre_FlexGMRESSetTol ( void *fgmres_vdata , double tol );
int hypre_FlexGMRESGetTol ( void *fgmres_vdata , double *tol );
int hypre_FlexGMRESSetAbsoluteTol ( void *fgmres_vdata , double a_tol );
int hypre_FlexGMRESGetAbsoluteTol ( void *fgmres_vdata , double *a_tol );
int hypre_FlexGMRESSetConvergenceFactorTol ( void *fgmres_vdata , double cf_tol );
int hypre_FlexGMRESGetConvergenceFactorTol ( void *fgmres_vdata , double *cf_tol );
int hypre_FlexGMRESSetMinIter ( void *fgmres_vdata , int min_iter );
int hypre_FlexGMRESGetMinIter ( void *fgmres_vdata , int *min_iter );
int hypre_FlexGMRESSetMaxIter ( void *fgmres_vdata , int max_iter );
int hypre_FlexGMRESGetMaxIter ( void *fgmres_vdata , int *max_iter );
int hypre_FlexGMRESSetStopCrit ( void *fgmres_vdata , int stop_crit );
int hypre_FlexGMRESGetStopCrit ( void *fgmres_vdata , int *stop_crit );
int hypre_FlexGMRESSetPrecond ( void *fgmres_vdata , int (*precond )(), int (*precond_setup )(), void *precond_data );
int hypre_FlexGMRESSetPrintLevel ( void *fgmres_vdata , int level );
int hypre_FlexGMRESGetPrintLevel ( void *fgmres_vdata , int *level );
int hypre_FlexGMRESSetLogging ( void *fgmres_vdata , int level );
int hypre_FlexGMRESGetLogging ( void *fgmres_vdata , int *level );
int hypre_FlexGMRESGetNumIterations ( void *fgmres_vdata , int *num_iterations );
int hypre_FlexGMRESGetConverged ( void *fgmres_vdata , int *converged );
int hypre_FlexGMRESGetFinalRelativeResidualNorm ( void *fgmres_vdata , double *relative_residual_norm );
int hypre_FlexGMRESSetModifyPC ( void *fgmres_vdata , int (*modify_pc )());
int hypre_FlexGMRESModifyPCDefault ( void *precond_data , int iteration , double rel_residual_norm );


#ifdef __cplusplus
}
#endif
#endif
