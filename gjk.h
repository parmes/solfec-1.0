/*
 * gjk.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------------------
 * distance between a pair of convex polyhedra according to the classical algortihm
 * by Gilbert et al. IEEE J. of Robotics and Automation, 4/2, 1988, pp. 193-203
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

#ifndef __gjk__
#define __gjk__

/* (a,na) and (b,nb) are the two input tables of polyhedrons vertices;
 * 'p' and 'q' are the two outputed closest points, respectively in
 * polyhedron (a,na) and polyhedron (b,nb); the distance is returned */
double gjk (double *a, int na, double *b, int nb, double *p, double *q);

/* (a,na) and (c,r) are the input polyhedron and sphere; 'p' and 'q' are the two outputed
 * closest points, respectively in polyhedron (a,na) and sphere (c,r); the distance is returned */
double gjk_convex_sphere (double *a, int na, double *c, double r, double *p, double *q);

/* (a,na) and p are the input polyhedron and point; 'q' is the outputed
 * closest point on the polyhedron; the distance is returned */
double gjk_convex_point (double *a, int na, double *p, double *q);

/* (a,na) and (b,bsca, brot) are the input polyhedron and ellipsoid; 'p' and 'q' are the two outputed
 * closest points, respectively in polyhedron (a,na) and ellipsoid (b, bsca, brot); the distance is returned */
double gjk_convex_ellip (double *a, int na, double *b, double *bsca, double *brot, double *p, double *q);

/* (a,ra) and (b,rb) are the input spheres; * 'p' and 'q' are the two outputed closest points,
 * respectively in spheres (a,ra) and (b,rb); the distance is returned */
double gjk_sphere_sphere (double *a, double ra, double *b, double rb, double *p, double *q);

/* (a,ra) and (b,bsca,brot)) are the input sphere and ellipsoid; 'p' and 'q' are the two outputed
 * closest points, respectively in sphere (a,ra) and ellipsoid (b,bsca,brot); the distance is returned */
double gjk_sphere_ellip (double *a, double ra, double *b, double *bsca, double *brot, double *p, double *q);

/* (a,asca,arot) and (b,bsca,brot) are the two input ellipsoids; 'p' and 'q' are the two outputed
 * closest points, respectively in (a,asca,arot) and (b,bsca,brot); the distance is returned */
double gjk_ellip_ellip (double *a, double *asca, double *arot, double *b, double *bsca, double *brot, double *p, double *q);

/* (a,asca,arot) and p are the input ellipsoid and point; 'q' is the outputed
 * closest point on the ellipsoid; the distance is returned */
double gjk_ellip_point (double *a, double *asca, double *arot, double *p, double *q);

#endif
