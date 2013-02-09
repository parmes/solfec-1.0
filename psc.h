/*
 * psc.h
 * Copyright (C) 2013, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * parallel self-consistency tests
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

#ifndef __psc__
#define __psc__

/* write body data to file before sending the body via MPI calles;
 * the file is SOLFEC->outpath/bodyID.data; */
void PSC_Write_Body (BODY *bod);

/* after the body has been imported via MPI calls,
 * read body data from file and compare */
void PSC_Test_Body (BODY *bod);

#endif
