/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include "taucs.h"

#ifndef TAUCS_CONFIG_TIMING
double taucs_wtime() { return 0.0; }
double taucs_ctime() { return 0.0; }
#else

#ifdef OSTYPE_WIN32
#define TAUCS_TIMER

#include <windows.h>

double taucs_wtime() {
  double wtime;

  SYSTEMTIME systime;

  GetSystemTime(&systime);

  wtime = 0.0
    +  0.001 * (double) systime.wMilliseconds
    +    1.0 * (double) systime.wSecond
    +   60.0 * (double) systime.wMinute
    + 3600.0 * (double) systime.wHour;

  return wtime; 
}
double taucs_ctime() { 
  double ctime;
  FILETIME creationt,exitt,kernel,user;
  ULARGE_INTEGER ukernel,uuser;
  HANDLE self;

  self = GetCurrentProcess();

  if (!GetProcessTimes(self,
		       &creationt,
		       &exitt,
		       &kernel,
		       &user))
    return 0.0;

  ukernel.LowPart  = kernel.dwLowDateTime;
  ukernel.HighPart = kernel.dwHighDateTime;
  uuser.LowPart  = user.dwLowDateTime;
  uuser.HighPart = user.dwHighDateTime;
  
  ctime = ((double) (signed __int64) ukernel.QuadPart / 10000000.0)
        + ((double) (signed __int64) uuser.QuadPart   / 10000000.0);

  CloseHandle( self );

  return ctime;
}
#endif /* win32 */

/*********************************************************/
/*                                                       */
/*********************************************************/

#ifndef OSTYPE_WIN32
#define TAUCS_TIMER

#include <stdio.h>                                                 
#include <unistd.h>

#ifdef OSTYPE_SOLARIS
#define _XPG4_2
#endif
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>                                                 
#include <sys/timeb.h>                                                 

double taucs_wtime()
{
  struct timeb T;
  static time_t time_diff;
  static time_t mill_diff;
  double dt;
  
  (void) ftime( &T );
  time_diff = T.time;
  mill_diff = T.millitm;

  dt = ((double) time_diff) + (1e-3) * ((double) mill_diff);

  return dt;
}

double taucs_ctime()
{
  struct rusage a;
  
  getrusage(RUSAGE_SELF,&a);

  return (double) 
    (double) ((a.ru_utime).tv_sec +(a.ru_stime).tv_sec ) +
    (double) ((a.ru_utime).tv_usec+(a.ru_stime).tv_usec) * 1.0e-6;
}

#endif /* not win32 */

/*********************************************************/
/*                                                       */
/*********************************************************/

#ifndef TAUCS_TIMER
#define TAUCS_TIMER
double taucs_wtime() { return 0.0; }
double taucs_ctime() { return 0.0; }

#endif

#endif /* TAUCS_CONFIG_TIMING */
