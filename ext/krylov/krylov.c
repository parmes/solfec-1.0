#include "krylov.h"

int hypre__global_error = 0;

/* Process the error with code ierr raised in the given line of the
   given source file. */
void hypre_error_handler(char *filename, int line, int ierr)
{
   hypre_error_flag |= ierr;

   fprintf(stderr,
           "hypre error in file \"%s\", line %d, error code = %d\n",
           filename, line, ierr);
}
