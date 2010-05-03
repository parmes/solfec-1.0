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
 * LGMRES lgmres
 *
 *****************************************************************************/

#ifndef hypre_KRYLOV_LGMRES_HEADER
#define hypre_KRYLOV_LGMRES_HEADER

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

/**
 * @name Generic LGMRES Interface
 *
 * A general description of the interface goes here...
 *
 **/
/*@{*/

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
 * hypre_LGMRESData and hypre_LGMRESFunctions
 *--------------------------------------------------------------------------*/


/**
 * @name LGMRES structs
 *
 * Description...
 **/
/*@{*/

/**
 * The {\tt hypre\_LGMRESFunctions} object ...
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

} hypre_LGMRESFunctions;

/**
 * The {\tt hypre\_LGMRESData} object ...
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

/*lgmres specific stuff */
   int      aug_dim;
   int      approx_constant;
   void   **aug_vecs;
   int     *aug_order;
   void   **a_aug_vecs;
/*---*/

   void  *A;
   void  *r;
   void  *w;
   void  *w_2;
   void  **p;

   void    *matvec_data;
   void    *precond_data;

   hypre_LGMRESFunctions * functions;

   /* log info (always logged) */
   int      num_iterations;
 
   int     print_level; /* printing when print_level>0 */
   int     logging;  /* extra computations for logging when logging>0 */
   double  *norms;
   char    *log_file_name;

} hypre_LGMRESData;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @name generic LGMRES Solver
 *
 * Description...
 **/
/*@{*/

/**
 * Description...
 *
 * @param param [IN] ...
 **/

hypre_LGMRESFunctions *
hypre_LGMRESFunctionsCreate(
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
hypre_LGMRESCreate( hypre_LGMRESFunctions *lgmres_functions );
int hypre_LGMRESDestroy ( void *lgmres_vdata );
int hypre_LGMRESGetResidual ( void *lgmres_vdata , void **residual );
int hypre_LGMRESSetup ( void *lgmres_vdata , void *A , void *b , void *x );
int hypre_LGMRESSolve ( void *lgmres_vdata , void *A , void *b , void *x );
int hypre_LGMRESSetKDim ( void *lgmres_vdata , int k_dim );
int hypre_LGMRESGetKDim ( void *lgmres_vdata , int *k_dim );
int hypre_LGMRESSetAugDim ( void *lgmres_vdata , int aug_dim );
int hypre_LGMRESGetAugDim ( void *lgmres_vdata , int *aug_dim );
int hypre_LGMRESSetTol ( void *lgmres_vdata , double tol );
int hypre_LGMRESGetTol ( void *lgmres_vdata , double *tol );
int hypre_LGMRESSetAbsoluteTol ( void *lgmres_vdata , double a_tol );
int hypre_LGMRESGetAbsoluteTol ( void *lgmres_vdata , double *a_tol );
int hypre_LGMRESSetConvergenceFactorTol ( void *lgmres_vdata , double cf_tol );
int hypre_LGMRESGetConvergenceFactorTol ( void *lgmres_vdata , double *cf_tol );
int hypre_LGMRESSetMinIter ( void *lgmres_vdata , int min_iter );
int hypre_LGMRESGetMinIter ( void *lgmres_vdata , int *min_iter );
int hypre_LGMRESSetMaxIter ( void *lgmres_vdata , int max_iter );
int hypre_LGMRESGetMaxIter ( void *lgmres_vdata , int *max_iter );
int hypre_LGMRESSetStopCrit ( void *lgmres_vdata , int stop_crit );
int hypre_LGMRESGetStopCrit ( void *lgmres_vdata , int *stop_crit );
int hypre_LGMRESSetPrecond ( void *lgmres_vdata , int (*precond )(), int (*precond_setup )(), void *precond_data );
int hypre_LGMRESSetPrintLevel ( void *lgmres_vdata , int level );
int hypre_LGMRESGetPrintLevel ( void *lgmres_vdata , int *level );
int hypre_LGMRESSetLogging ( void *lgmres_vdata , int level );
int hypre_LGMRESGetLogging ( void *lgmres_vdata , int *level );
int hypre_LGMRESGetNumIterations ( void *lgmres_vdata , int *num_iterations );
int hypre_LGMRESGetConverged ( void *lgmres_vdata , int *converged );
int hypre_LGMRESGetFinalRelativeResidualNorm ( void *lgmres_vdata , double *relative_residual_norm );


#ifdef __cplusplus
}
#endif
#endif
