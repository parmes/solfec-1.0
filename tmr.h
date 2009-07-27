/*
 * tmr.h
 * Copyright (C) 2005, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * timer
 */

/* This file is part of Solfec.
 * Solfec is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * Solfec is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Solfec. If not, see <http://www.gnu.org/licenses/>. */

#include <sys/time.h>

#ifndef __tmr__
#define __tmr__

typedef struct timing TIMING;

struct timing
{
  struct timeval time;
  double sec;
};

static inline void timerstart (TIMING *t)
{
  gettimeofday (&t->time, NULL);
}

static inline double timerend (TIMING *t)
{
  struct timeval newtime;
  gettimeofday (&newtime, NULL);
  t->sec = (((double)newtime.tv_sec - (double)t->time.tv_sec) * 1000000. +
    ((double)newtime.tv_usec - (double)t->time.tv_usec)) / 1000000.;
  return t->sec;
}

#endif
