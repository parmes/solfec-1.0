/*
 * solfec.h
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * solfec main module
 */

#include "sol.h"

#ifndef __solfec__
#define __solfec__

#if __cplusplus
extern "C" { /* C */
#endif

/* register new SOLFEC object */
void REGISTER_SOLFEC (SOLFEC *sol);

/* get output directory */
char* OUTPUT_DIR ();

/* get output sub-directory */
char* OUTPUT_SUBDIR ();

/* get input file path */
char* INPUT_FILE ();

/* get write mode flag */
int WRITE_MODE_FLAG ();

/* get continue write flag */
int CONTINUE_WRITE_FLAG ();

/* get wireframe flag */
int WIREFRAME_FLAG ();

/* verbose output time interval */
double VERBOSITY_INTERVAL ();

/* get non-Solfec input arguments */
char** NON_SOLFEC_ARGV (int *argc);

#if __cplusplus
} /* extern C */
#endif

#endif
