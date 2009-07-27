/*
 * lis.h
 * Copyright (C) 2008, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------------
 * list merge sort
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

#ifndef __lis__
#define __lis__

#define DOUBLY_LINKED(PREV, NEXT) for (i = list; i; j = i, i = i->NEXT) i->PREV = j;

#define SINGLE_LINKED(PREV, NEXT)

#define IMPLEMENT_LIST_SORT(KIND, CALL, LIST, PREV, NEXT, LE)\
static LIST* CALL (LIST *list)\
{\
  LIST *i, *j, *k, *h, *t;\
  int l, m, n;\
\
  for (l = 1;;l *= 2)\
  {\
    h = t = NULL;\
\
    for (j = list;;)\
    {\
      i = j;\
\
      for (m = 0; m < l && j; j = j->NEXT, m ++);\
      for (n = 0, k = j; n < l && k; k = k->NEXT, n ++);\
\
      if (!j && i == list)\
      {\
	KIND (PREV, NEXT)\
	return list;\
      }\
      else if (!(m+n)) break;\
\
      if (!h) h = (LE (i, j) ? i : j);\
\
      for (; m && n;)\
      {\
	if (LE (i, j))\
	{\
	  if (t) t->NEXT = i;\
	  t = i;\
	  i = i->NEXT;\
	  m --;\
	}\
	else\
	{\
	  if (t) t->NEXT = j;\
	  t = j;\
	  j = j->NEXT;\
	  n --;\
	}\
      }\
\
      while (m)\
      {\
	t->NEXT = i;\
	t = i;\
	i = i->NEXT;\
	m --;\
      }\
\
      while (n)\
      {\
	t->NEXT = j;\
	t = j;\
	j = j->NEXT;\
	n --;\
      }\
    }\
\
    t->NEXT = NULL;\
    list = h;\
  }\
\
  return list;\
}

#endif
