/*
 * lng.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * Python based input language parser
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

#ifndef __lng__
#define __lng__

#if __cplusplus
extern "C" { /* C */
#endif

/* interpret an input file (return 0 on success) */
int lng (const char *path);

/* finalize interpreter */
void lngfinalize ();

/* get positive id of a callback pointer pair;
 * return 0 if the pair was not found  */
int  lngcallback_id (void *data, void *call);

/* set callback pointers for a given id; return 1 when
 * the id was found, or return 0 otherwise */
int  lngcallback_set (int id, void **data, void **call);

/* handle PUT_SPRING spring Python callback */
double springcallback (void *call, double stroke, double velocity);

/* Python object referece increment */
void lng_xincref (void *obj);

/* Python object referece decrement */
void lng_xdecref (void *obj);

/* append Python list with a list storing a vector of double precision numbers */
void lng_append_list_with_list_of_doubles (void *list, double *vector, int length);

#if __cplusplus
} /* extern C */
#endif

#endif
