/*
 * svk.c
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * Saint Venant - Kirchhoff material
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

#include "svk.h"

double SVK_Energy_R (double lambda, double mi, double volume, double *F)
{
  double E [9], S [9], trace;

  /* calculate Green tensor: E = (FF - 1) / 2 */
  E [0] = .5*(F [0]*F [0] + F [3]*F [3] + F [6]*F [6] - 1.);
  E [1] = .5*(F [0]*F [1] + F [3]*F [4] + F [6]*F [7]);
  E [2] = .5*(F [0]*F [2] + F [3]*F [5] + F [6]*F [8]);
  E [3] = .5*(F [1]*F [0] + F [4]*F [3] + F [7]*F [6]);
  E [4] = .5*(F [1]*F [1] + F [4]*F [4] + F [7]*F [7] - 1.);
  E [5] = .5*(F [1]*F [2] + F [4]*F [5] + F [7]*F [8]);
  E [6] = .5*(F [2]*F [0] + F [5]*F [3] + F [8]*F [6]);
  E [7] = .5*(F [2]*F [1] + F [5]*F [4] + F [8]*F [7]);
  E [8] = .5*(F [2]*F [2] + F [5]*F [5] + F [8]*F [8] - 1.);

  /* obtain the second PK tensor
   * trough the Saint Venant - Kirchhoff law */
  trace = 2. * mi;
  S [0] = trace * E [0];
  S [1] = trace * E [1];
  S [2] = trace * E [2];
  S [3] = trace * E [3];
  S [4] = trace * E [4];
  S [5] = trace * E [5];
  S [6] = trace * E [6];
  S [7] = trace * E [7];
  S [8] = trace * E [8];

  trace = E [0] + E [4] + E [8];
  S [0] += lambda * trace;
  S [4] += lambda * trace;
  S [8] += lambda * trace;

  /* internal energy: 0.5 * E : C : E */
  return 0.5 * volume * (E[0]*S[0]+E[1]*S[1]+E[2]*S[2]+
    E[3]*S[3]+E[4]*S[4]+E[5]*S[5]+E[6]*S[6]+E[7]*S[7]+E[8]*S[8]);
}

double SVK_Energy_C (double lambda, double mi, double volume, double *F)
{
  double E [9], S [9], trace;

  /* calculate Green tensor: E = (FF - 1) / 2 */
  E [0] = .5 * (F [0]*F [0] + F [1]*F [1] + F [2]*F [2] - 1.);
  E [1] = .5 * (F [3]*F [0] + F [4]*F [1] + F [5]*F [2]);
  E [2] = .5 * (F [6]*F [0] + F [7]*F [1] + F [8]*F [2]);
  E [3] = .5 * (F [0]*F [3] + F [1]*F [4] + F [2]*F [5]);
  E [4] = .5 * (F [3]*F [3] + F [4]*F [4] + F [5]*F [5] - 1.);
  E [5] = .5 * (F [6]*F [3] + F [7]*F [4] + F [8]*F [5]);
  E [6] = .5 * (F [0]*F [6] + F [1]*F [7] + F [2]*F [8]);
  E [7] = .5 * (F [3]*F [6] + F [4]*F [7] + F [5]*F [8]);
  E [8] = .5 * (F [6]*F [6] + F [7]*F [7] + F [8]*F [8] - 1.);

  /* obtain the second PK tensor
   * trough the Saint Venant - Kirchhoff law */
  trace = 2. * mi;
  S [0] = trace * E [0];
  S [1] = trace * E [1];
  S [2] = trace * E [2];
  S [3] = trace * E [3];
  S [4] = trace * E [4];
  S [5] = trace * E [5];
  S [6] = trace * E [6];
  S [7] = trace * E [7];
  S [8] = trace * E [8];

  trace = E [0] + E [4] + E [8];
  S [0] += lambda * trace;
  S [4] += lambda * trace;
  S [8] += lambda * trace;

  /* internal energy: 0.5 * E : C : E */
  return 0.5 * volume * (E[0]*S[0]+E[1]*S[1]+E[2]*S[2]+
    E[3]*S[3]+E[4]*S[4]+E[5]*S[5]+E[6]*S[6]+E[7]*S[7]+E[8]*S[8]);
}

double SVK_Stress_R (double lambda, double mi, double volume, double *F, double *P)
{
  double J, E [9], S [9], trace;

  /* J = det (F) */
  J = F [0]*F [4]*F [8] + F [3]*F [7]*F [2] + F [6]*F [1]*F [5] -
      F [6]*F [4]*F [2] - F [0]*F [7]*F [5] - F [3]*F [1]*F [8];

  /* calculate Green tensor: E = (FF - 1) / 2 */
  E [0] = .5*(F [0]*F [0] + F [3]*F [3] + F [6]*F [6] - 1.);
  E [1] = .5*(F [0]*F [1] + F [3]*F [4] + F [6]*F [7]);
  E [2] = .5*(F [0]*F [2] + F [3]*F [5] + F [6]*F [8]);
  E [3] = .5*(F [1]*F [0] + F [4]*F [3] + F [7]*F [6]);
  E [4] = .5*(F [1]*F [1] + F [4]*F [4] + F [7]*F [7] - 1.);
  E [5] = .5*(F [1]*F [2] + F [4]*F [5] + F [7]*F [8]);
  E [6] = .5*(F [2]*F [0] + F [5]*F [3] + F [8]*F [6]);
  E [7] = .5*(F [2]*F [1] + F [5]*F [4] + F [8]*F [7]);
  E [8] = .5*(F [2]*F [2] + F [5]*F [5] + F [8]*F [8] - 1.);

  /* obtain the second PK tensor
   * trough the Saint Venant - Kirchhoff law */
  trace = 2. * mi;
  S [0] = trace * E [0];
  S [1] = trace * E [1];
  S [2] = trace * E [2];
  S [3] = trace * E [3];
  S [4] = trace * E [4];
  S [5] = trace * E [5];
  S [6] = trace * E [6];
  S [7] = trace * E [7];
  S [8] = trace * E [8];

  trace = E [0] + E [4] + E [8];
  S [0] += lambda * trace;
  S [4] += lambda * trace;
  S [8] += lambda * trace;

  /* now conver S - which is the second PK tensor
   * into the first PK tensor: P = F S */
  P [0] = volume * (F [0]*S [0] + F [1]*S [3] + F [2]*S [6]);
  P [1] = volume * (F [0]*S [1] + F [1]*S [4] + F [2]*S [7]);
  P [2] = volume * (F [0]*S [2] + F [1]*S [5] + F [2]*S [8]);
  P [3] = volume * (F [3]*S [0] + F [4]*S [3] + F [5]*S [6]);
  P [4] = volume * (F [3]*S [1] + F [4]*S [4] + F [5]*S [7]);
  P [5] = volume * (F [3]*S [2] + F [4]*S [5] + F [5]*S [8]);
  P [6] = volume * (F [6]*S [0] + F [7]*S [3] + F [8]*S [6]);
  P [7] = volume * (F [6]*S [1] + F [7]*S [4] + F [8]*S [7]);
  P [8] = volume * (F [6]*S [2] + F [7]*S [5] + F [8]*S [8]);

  return J;
}

double SVK_Stress_C (double lambda, double mi, double volume, double *F, double *P)
{
  double J, E [9], S [9], trace;

  /* J = det (F) */
  J = F [0]*F [4]*F [8] + F [3]*F [7]*F [2] + F [6]*F [1]*F [5] -
      F [6]*F [4]*F [2] - F [0]*F [7]*F [5] - F [3]*F [1]*F [8];

  /* calculate Green tensor: E = (FF - 1) / 2 */
  E [0] = .5 * (F [0]*F [0] + F [1]*F [1] + F [2]*F [2] - 1.);
  E [1] = .5 * (F [3]*F [0] + F [4]*F [1] + F [5]*F [2]);
  E [2] = .5 * (F [6]*F [0] + F [7]*F [1] + F [8]*F [2]);
  E [3] = .5 * (F [0]*F [3] + F [1]*F [4] + F [2]*F [5]);
  E [4] = .5 * (F [3]*F [3] + F [4]*F [4] + F [5]*F [5] - 1.);
  E [5] = .5 * (F [6]*F [3] + F [7]*F [4] + F [8]*F [5]);
  E [6] = .5 * (F [0]*F [6] + F [1]*F [7] + F [2]*F [8]);
  E [7] = .5 * (F [3]*F [6] + F [4]*F [7] + F [5]*F [8]);
  E [8] = .5 * (F [6]*F [6] + F [7]*F [7] + F [8]*F [8] - 1.);

  /* obtain the second PK tensor
   * trough the Saint Venant - Kirchhoff law */
  trace = 2. * mi;
  S [0] = trace * E [0];
  S [1] = trace * E [1];
  S [2] = trace * E [2];
  S [3] = trace * E [3];
  S [4] = trace * E [4];
  S [5] = trace * E [5];
  S [6] = trace * E [6];
  S [7] = trace * E [7];
  S [8] = trace * E [8];

  trace = E [0] + E [4] + E [8];
  S [0] += lambda * trace;
  S [4] += lambda * trace;
  S [8] += lambda * trace;

  /* now conver S - which is the second PK tensor
   * into the first PK tensor: P = F S */
  P [0] = volume * (F [0]*S [0] + F [3]*S [1] + F [6]*S [2]);
  P [1] = volume * (F [1]*S [0] + F [4]*S [1] + F [7]*S [2]);
  P [2] = volume * (F [2]*S [0] + F [5]*S [1] + F [8]*S [2]);
  P [3] = volume * (F [0]*S [3] + F [3]*S [4] + F [6]*S [5]);
  P [4] = volume * (F [1]*S [3] + F [4]*S [4] + F [7]*S [5]);
  P [5] = volume * (F [2]*S [3] + F [5]*S [4] + F [8]*S [5]);
  P [6] = volume * (F [0]*S [6] + F [3]*S [7] + F [6]*S [8]);
  P [7] = volume * (F [1]*S [6] + F [4]*S [7] + F [7]*S [8]);
  P [8] = volume * (F [2]*S [6] + F [5]*S [7] + F [8]*S [8]);

  return J;
}

/* row-wise F => .................................. */

#define SQ(X) ((X)*(X))
static int PERM [] = {0, 3, 6, 1, 4, 7, 2, 5, 8};
#define F(I, J) F[3 * (I-1) + (J-1)]
#define K(I, J) K[PERM[I-1] + dim * PERM[J-1]]
/* things below were generated with Maxima's 'itensor' package,
 * there is not much to comment on apart from that */
void SVK_Tangent_R (double lambda, double mi, double volume, int dim, double *F, double *K)
{
  K (1, 1) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+3*SQ(F(1,1))-3)*lambda + 
             (2*SQ(F(3,1))+2*SQ(F(2,1))+2*SQ(F(1,3))+2*SQ(F(1,2))+6*SQ(F(1,1))-2)*mi)*0.5;
  K (2, 1) = F(1,1)*F(2,1)*lambda+(F(1,3)*F(2,3)+F(1,2)*F(2,2)+2*F(1,1)*F(2,1))*mi;
  K (3, 1) = F(1,1)*F(3,1)*lambda+(F(1,3)*F(3,3)+F(1,2)*F(3,2)+2*F(1,1)*F(3,1))*mi;
  K (4, 1) = F(1,1)*F(1,2)*lambda+(F(3,1)*F(3,2)+F(2,1)*F(2,2)+2*F(1,1)*F(1,2))*mi;
  K (5, 1) = F(1,1)*F(2,2)*lambda+F(1,2)*F(2,1)*mi;
  K (6, 1) = F(1,1)*F(3,2)*lambda+F(1,2)*F(3,1)*mi;
  K (7, 1) = F(1,1)*F(1,3)*lambda+(F(3,1)*F(3,3)+F(2,1)*F(2,3)+2*F(1,1)*F(1,3))*mi;
  K (8, 1) = F(1,1)*F(2,3)*lambda+F(1,3)*F(2,1)*mi;
  K (9, 1) = F(1,1)*F(3,3)*lambda+F(1,3)*F(3,1)*mi;
  K (1, 2) = K(2,1);
  K (2, 2) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+3*SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,1))+2*SQ(F(2,3))+2*SQ(F(2,2))+6*SQ(F(2,1))+2*SQ(F(1,1))-2)*mi)*0.5;
  K (3, 2) = F(2,1)*F(3,1)*lambda+(F(2,3)*F(3,3)+F(2,2)*F(3,2)+2*F(2,1)*F(3,1))*mi;
  K (4, 2) = F(1,2)*F(2,1)*lambda+F(1,1)*F(2,2)*mi;
  K (5, 2) = F(2,1)*F(2,2)*lambda+(F(3,1)*F(3,2)+2*F(2,1)*F(2,2)+F(1,1)*F(1,2))*mi;
  K (6, 2) = F(2,1)*F(3,2)*lambda+F(2,2)*F(3,1)*mi;
  K (7, 2) = F(1,3)*F(2,1)*lambda+F(1,1)*F(2,3)*mi;
  K (8, 2) = F(2,1)*F(2,3)*lambda+(F(3,1)*F(3,3)+2*F(2,1)*F(2,3)+F(1,1)*F(1,3))*mi;
  K (9, 2) = F(2,1)*F(3,3)*lambda+F(2,3)*F(3,1)*mi;
  K (1, 3) = K(3,1);
  K (2, 3) = K(3,2);
  K (3, 3) = ((SQ(F(3,3))+SQ(F(3,2))+3*SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,3))+2*SQ(F(3,2))+6*SQ(F(3,1))+2*SQ(F(2,1))+2*SQ(F(1,1))-2)*mi)*0.5;
  K (4, 3) = F(1,2)*F(3,1)*lambda+F(1,1)*F(3,2)*mi;
  K (5, 3) = F(2,2)*F(3,1)*lambda+F(2,1)*F(3,2)*mi;
  K (6, 3) = F(3,1)*F(3,2)*lambda+(2*F(3,1)*F(3,2)+F(2,1)*F(2,2)+F(1,1)*F(1,2))*mi;
  K (7, 3) = F(1,3)*F(3,1)*lambda+F(1,1)*F(3,3)*mi;
  K (8, 3) = F(2,3)*F(3,1)*lambda+F(2,1)*F(3,3)*mi;
  K (9, 3) = F(3,1)*F(3,3)*lambda+(2*F(3,1)*F(3,3)+F(2,1)*F(2,3)+F(1,1)*F(1,3))*mi;
  K (1, 4) = K(4,1);
  K (2, 4) = K(4,2);
  K (3, 4) = K(4,3);
  K (4, 4) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+3*SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,2))+2*SQ(F(2,2))+2*SQ(F(1,3))+6*SQ(F(1,2))+2*SQ(F(1,1))-2)*mi)*0.5;
  K (5, 4) = F(1,2)*F(2,2)*lambda+(F(1,3)*F(2,3)+2*F(1,2)*F(2,2)+F(1,1)*F(2,1))*mi;
  K (6, 4) = F(1,2)*F(3,2)*lambda+(F(1,3)*F(3,3)+2*F(1,2)*F(3,2)+F(1,1)*F(3,1))*mi;
  K (7, 4) = F(1,2)*F(1,3)*lambda+(F(3,2)*F(3,3)+F(2,2)*F(2,3)+2*F(1,2)*F(1,3))*mi;
  K (8, 4) = F(1,2)*F(2,3)*lambda+F(1,3)*F(2,2)*mi;
  K (9, 4) = F(1,2)*F(3,3)*lambda+F(1,3)*F(3,2)*mi;
  K (1, 5) = K(5,1);
  K (2, 5) = K(5,2);
  K (3, 5) = K(5,3);
  K (4, 5) = K(5,4);
  K (5, 5) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+3*SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,2))+2*SQ(F(2,3))+6*SQ(F(2,2))+2*SQ(F(2,1))+2*SQ(F(1,2))-2)*mi)*0.5;
  K (6, 5) = F(2,2)*F(3,2)*lambda+(F(2,3)*F(3,3)+2*F(2,2)*F(3,2)+F(2,1)*F(3,1))*mi;
  K (7, 5) = F(1,3)*F(2,2)*lambda+F(1,2)*F(2,3)*mi;
  K (8, 5) = F(2,2)*F(2,3)*lambda+(F(3,2)*F(3,3)+2*F(2,2)*F(2,3)+F(1,2)*F(1,3))*mi;
  K (9, 5) = F(2,2)*F(3,3)*lambda+F(2,3)*F(3,2)*mi;
  K (1, 6) = K(6,1);
  K (2, 6) = K(6,2);
  K (3, 6) = K(6,3);
  K (4, 6) = K(6,4);
  K (5, 6) = K(6,5);
  K (6, 6) = ((SQ(F(3,3))+3*SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda + 
             (2*SQ(F(3,3))+6*SQ(F(3,2))+2*SQ(F(3,1))+2*SQ(F(2,2))+2*SQ(F(1,2))-2)*mi)*0.5;
  K (7, 6) = F(1,3)*F(3,2)*lambda+F(1,2)*F(3,3)*mi;
  K (8, 6) = F(2,3)*F(3,2)*lambda+F(2,2)*F(3,3)*mi;
  K (9, 6) = F(3,2)*F(3,3)*lambda+(2*F(3,2)*F(3,3)+F(2,2)*F(2,3)+F(1,2)*F(1,3))*mi;
  K (1, 7) = K(7,1);
  K (2, 7) = K(7,2);
  K (3, 7) = K(7,3);
  K (4, 7) = K(7,4);
  K (5, 7) = K(7,5);
  K (6, 7) = K(7,6);
  K (7, 7) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+3*SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +  
             (2*SQ(F(3,3))+2*SQ(F(2,3))+6*SQ(F(1,3))+2*SQ(F(1,2))+2*SQ(F(1,1))-2)*mi)*0.5;
  K (8, 7) = F(1,3)*F(2,3)*lambda+(2*F(1,3)*F(2,3)+F(1,2)*F(2,2)+F(1,1)*F(2,1))*mi;
  K (9, 7) = F(1,3)*F(3,3)*lambda+(2*F(1,3)*F(3,3)+F(1,2)*F(3,2)+F(1,1)*F(3,1))*mi;
  K (1, 8) = K(8,1);
  K (2, 8) = K(8,2);
  K (3, 8) = K(8,3);
  K (4, 8) = K(8,4);
  K (5, 8) = K(8,5);
  K (6, 8) = K(8,6);
  K (7, 8) = K(8,7);
  K (8, 8) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+3*SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +   
             (2*SQ(F(3,3))+6*SQ(F(2,3))+2*SQ(F(2,2))+2*SQ(F(2,1))+2*SQ(F(1,3))-2)*mi)*0.5;
  K (9, 8) = F(2,3)*F(3,3)*lambda+(2*F(2,3)*F(3,3)+F(2,2)*F(3,2)+F(2,1)*F(3,1))*mi;
  K (1, 9) = K(9,1);
  K (2, 9) = K(9,2);
  K (3, 9) = K(9,3);
  K (4, 9) = K(9,4);
  K (5, 9) = K(9,5);
  K (6, 9) = K(9,6);
  K (7, 9) = K(9,7);
  K (8, 9) = K(9,8);
  K (9, 9) = ((3*SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +    
             (6*SQ(F(3,3))+2*SQ(F(3,2))+2*SQ(F(3,1))+2*SQ(F(2,3))+2*SQ(F(1,3))-2)*mi)*0.5;

  /* volume scaling */
  if (volume != 1.0)
  {
    for (int col = 1; col <= 9; col ++)
    {
      double *val = & K(1, col); /* unroll */
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
    }
  }
}

void SVK_Tangent_Diagonal_R (double lambda, double mi, double volume, int dim, double *F, double *K)
{
  K [0] = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+3*SQ(F(1,1))-3)*lambda + 
             (2*SQ(F(3,1))+2*SQ(F(2,1))+2*SQ(F(1,3))+2*SQ(F(1,2))+6*SQ(F(1,1))-2)*mi)*0.5;
  K [1] = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+3*SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,1))+2*SQ(F(2,3))+2*SQ(F(2,2))+6*SQ(F(2,1))+2*SQ(F(1,1))-2)*mi)*0.5;
  K [2] = ((SQ(F(3,3))+SQ(F(3,2))+3*SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,3))+2*SQ(F(3,2))+6*SQ(F(3,1))+2*SQ(F(2,1))+2*SQ(F(1,1))-2)*mi)*0.5;
  K [3] = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+3*SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,2))+2*SQ(F(2,2))+2*SQ(F(1,3))+6*SQ(F(1,2))+2*SQ(F(1,1))-2)*mi)*0.5;
  K [4] = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+3*SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,2))+2*SQ(F(2,3))+6*SQ(F(2,2))+2*SQ(F(2,1))+2*SQ(F(1,2))-2)*mi)*0.5;
  K [5] = ((SQ(F(3,3))+3*SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda + 
             (2*SQ(F(3,3))+6*SQ(F(3,2))+2*SQ(F(3,1))+2*SQ(F(2,2))+2*SQ(F(1,2))-2)*mi)*0.5;
  K [6] = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+3*SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +  
             (2*SQ(F(3,3))+2*SQ(F(2,3))+6*SQ(F(1,3))+2*SQ(F(1,2))+2*SQ(F(1,1))-2)*mi)*0.5;
  K [7] = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+3*SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +   
             (2*SQ(F(3,3))+6*SQ(F(2,3))+2*SQ(F(2,2))+2*SQ(F(2,1))+2*SQ(F(1,3))-2)*mi)*0.5;
  K [8] = ((3*SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +    
             (6*SQ(F(3,3))+2*SQ(F(3,2))+2*SQ(F(3,1))+2*SQ(F(2,3))+2*SQ(F(1,3))-2)*mi)*0.5;

  if (volume != 1.0)
  {
    K [0] *= volume;
    K [1] *= volume;
    K [2] *= volume;
    K [3] *= volume;
    K [4] *= volume;
    K [5] *= volume;
    K [6] *= volume;
    K [7] *= volume;
    K [8] *= volume;
  }
}

void SVK_Tangent_Derivative_R (int comp,
 double lambda, double mi, double volume, int dim, double *F, double *K)
{
  switch (comp)
  {
  case 0:
    K (1, 1) = 3*F(1,1)*lambda+6*F(1,1)*mi;
    K (2, 1) = F(2,1)*lambda+2*F(2,1)*mi;
    K (3, 1) = F(3,1)*lambda+2*F(3,1)*mi;
    K (4, 1) = F(1,2)*lambda+2*F(1,2)*mi;
    K (5, 1) = F(2,2)*lambda;
    K (6, 1) = F(3,2)*lambda;
    K (7, 1) = F(1,3)*lambda+2*F(1,3)*mi;
    K (8, 1) = F(2,3)*lambda;
    K (9, 1) = F(3,3)*lambda;
    K (2, 2) = F(1,1)*lambda+2*F(1,1)*mi;
    K (3, 2) = 0;
    K (4, 2) = F(2,2)*mi;
    K (5, 2) = F(1,2)*mi;
    K (6, 2) = 0;
    K (7, 2) = F(2,3)*mi;
    K (8, 2) = F(1,3)*mi;
    K (9, 2) = 0;
    K (3, 3) = F(1,1)*lambda+2*F(1,1)*mi;
    K (4, 3) = F(3,2)*mi;
    K (5, 3) = 0;
    K (6, 3) = F(1,2)*mi;
    K (7, 3) = F(3,3)*mi;
    K (8, 3) = 0;
    K (9, 3) = F(1,3)*mi;
    K (4, 4) = F(1,1)*lambda+2*F(1,1)*mi;
    K (5, 4) = F(2,1)*mi;
    K (6, 4) = F(3,1)*mi;
    K (7, 4) = 0;
    K (8, 4) = 0;
    K (9, 4) = 0;
    K (5, 5) = F(1,1)*lambda;
    K (6, 5) = 0;
    K (7, 5) = 0;
    K (8, 5) = 0;
    K (9, 5) = 0;
    K (6, 6) = F(1,1)*lambda;
    K (7, 6) = 0;
    K (8, 6) = 0;
    K (9, 6) = 0;
    K (7, 7) = F(1,1)*lambda+2*F(1,1)*mi;
    K (8, 7) = F(2,1)*mi;
    K (9, 7) = F(3,1)*mi;
    K (8, 8) = F(1,1)*lambda;
    K (9, 8) = 0;
    K (9, 9) = F(1,1)*lambda;
    break;
  case 1:
    K (1, 1) = F(2,1)*lambda+2*F(2,1)*mi;
    K (2, 1) = F(1,1)*lambda+2*F(1,1)*mi;
    K (3, 1) = 0;
    K (4, 1) = F(2,2)*mi;
    K (5, 1) = F(1,2)*mi;
    K (6, 1) = 0;
    K (7, 1) = F(2,3)*mi;
    K (8, 1) = F(1,3)*mi;
    K (9, 1) = 0;
    K (2, 2) = 3*F(2,1)*lambda+6*F(2,1)*mi;
    K (3, 2) = F(3,1)*lambda+2*F(3,1)*mi;
    K (4, 2) = F(1,2)*lambda;
    K (5, 2) = F(2,2)*lambda+2*F(2,2)*mi;
    K (6, 2) = F(3,2)*lambda;
    K (7, 2) = F(1,3)*lambda;
    K (8, 2) = F(2,3)*lambda+2*F(2,3)*mi;
    K (9, 2) = F(3,3)*lambda;
    K (3, 3) = F(2,1)*lambda+2*F(2,1)*mi;
    K (4, 3) = 0 ;
    K (5, 3) = F(3,2)*mi;
    K (6, 3) = F(2,2)*mi;
    K (7, 3) = 0;
    K (8, 3) = F(3,3)*mi;
    K (9, 3) = F(2,3)*mi;
    K (4, 4) = F(2,1)*lambda;
    K (5, 4) = F(1,1)*mi;
    K (6, 4) = 0;
    K (7, 4) = 0;
    K (8, 4) = 0;
    K (9, 4) = 0;
    K (5, 5) = F(2,1)*lambda+2*F(2,1)*mi;
    K (6, 5) = F(3,1)*mi;
    K (7, 5) = 0;
    K (8, 5) = 0;
    K (9, 5) = 0;
    K (6, 6) = F(2,1)*lambda;
    K (7, 6) = 0;
    K (8, 6) = 0;
    K (9, 6) = 0;
    K (7, 7) = F(2,1)*lambda;
    K (8, 7) = F(1,1)*mi;
    K (9, 7) = 0;
    K (8, 8) = F(2,1)*lambda+2*F(2,1)*mi;
    K (9, 8) = F(3,1)*mi ;
    K (9, 9) = F(2,1)*lambda ;
    break;
  case 2:
    K (1, 1) = F(3,1)*lambda+2*F(3,1)*mi;
    K (2, 1) = 0;
    K (3, 1) = F(1,1)*lambda+2*F(1,1)*mi;
    K (4, 1) = F(3,2)*mi;
    K (5, 1) = 0;
    K (6, 1) = F(1,2)*mi;
    K (7, 1) = F(3,3)*mi;
    K (8, 1) = 0;
    K (9, 1) = F(1,3)*mi;
    K (2, 2) = F(3,1)*lambda+2*F(3,1)*mi;
    K (3, 2) = F(2,1)*lambda+2*F(2,1)*mi;
    K (4, 2) = 0;
    K (5, 2) = F(3,2)*mi;
    K (6, 2) = F(2,2)*mi;
    K (7, 2) = 0;
    K (8, 2) = F(3,3)*mi;
    K (9, 2) = F(2,3)*mi;
    K (3, 3) = 3*F(3,1)*lambda+6*F(3,1)*mi;
    K (4, 3) = F(1,2)*lambda;
    K (5, 3) = F(2,2)*lambda;
    K (6, 3) = F(3,2)*lambda+2*F(3,2)*mi;
    K (7, 3) = F(1,3)*lambda;
    K (8, 3) = F(2,3)*lambda;
    K (9, 3) = F(3,3)*lambda+2*F(3,3)*mi;
    K (4, 4) = F(3,1)*lambda;
    K (5, 4) = 0;
    K (6, 4) = F(1,1)*mi;
    K (7, 4) = 0;
    K (8, 4) = 0;
    K (9, 4) = 0;
    K (5, 5) = F(3,1)*lambda;
    K (6, 5) = F(2,1)*mi;
    K (7, 5) = 0;
    K (8, 5) = 0;
    K (9, 5) = 0;
    K (6, 6) = F(3,1)*lambda+2*F(3,1)*mi;
    K (7, 6) = 0;
    K (8, 6) = 0;
    K (9, 6) = 0;
    K (7, 7) = F(3,1)*lambda;
    K (8, 7) = 0;
    K (9, 7) = F(1,1)*mi;
    K (8, 8) = F(3,1)*lambda;
    K (9, 8) = F(2,1)*mi;
    K (9, 9) = F(3,1)*lambda+2*F(3,1)*mi;
    break;
  case 3:
    K (1, 1) = F(1,2)*lambda+2*F(1,2)*mi;
    K (2, 1) = F(2,2)*mi;
    K (3, 1) = F(3,2)*mi;
    K (4, 1) = F(1,1)*lambda+2*F(1,1)*mi;
    K (5, 1) = F(2,1)*mi;
    K (6, 1) = F(3,1)*mi;
    K (7, 1) = 0;
    K (8, 1) = 0;
    K (9, 1) = 0;
    K (2, 2) = F(1,2)*lambda;
    K (3, 2) = 0;
    K (4, 2) = F(2,1)*lambda;
    K (5, 2) = F(1,1)*mi;
    K (6, 2) = 0;
    K (7, 2) = 0;
    K (8, 2) = 0;
    K (9, 2) = 0;
    K (3, 3) = F(1,2)*lambda;
    K (4, 3) = F(3,1)*lambda;
    K (5, 3) = 0;
    K (6, 3) = F(1,1)*mi;
    K (7, 3) = 0;
    K (8, 3) = 0;
    K (9, 3) = 0;
    K (4, 4) = 3*F(1,2)*lambda+6*F(1,2)*mi;
    K (5, 4) = F(2,2)*lambda+2*F(2,2)*mi;
    K (6, 4) = F(3,2)*lambda+2*F(3,2)*mi;
    K (7, 4) = F(1,3)*lambda+2*F(1,3)*mi;
    K (8, 4) = F(2,3)*lambda;
    K (9, 4) = F(3,3)*lambda;
    K (5, 5) = F(1,2)*lambda+2*F(1,2)*mi;
    K (6, 5) = 0;
    K (7, 5) = F(2,3)*mi;
    K (8, 5) = F(1,3)*mi;
    K (9, 5) = 0;
    K (6, 6) = F(1,2)*lambda+2*F(1,2)*mi;
    K (7, 6) = F(3,3)*mi;
    K (8, 6) = 0;
    K (9, 6) = F(1,3)*mi;
    K (7, 7) = F(1,2)*lambda+2*F(1,2)*mi;
    K (8, 7) = F(2,2)*mi;
    K (9, 7) = F(3,2)*mi;
    K (8, 8) = F(1,2)*lambda;
    K (9, 8) = 0;
    K (9, 9) = F(1,2)*lambda;
    break;
  case 4:
    K (1, 1) = F(2,2)*lambda;
    K (2, 1) = F(1,2)*mi;
    K (3, 1) = 0;
    K (4, 1) = F(2,1)*mi;
    K (5, 1) = F(1,1)*lambda;
    K (6, 1) = 0;
    K (7, 1) = 0;
    K (8, 1) = 0;
    K (9, 1) = 0;
    K (2, 2) = F(2,2)*lambda+2*F(2,2)*mi;
    K (3, 2) = F(3,2)*mi;
    K (4, 2) = F(1,1)*mi;
    K (5, 2) = F(2,1)*lambda+2*F(2,1)*mi;
    K (6, 2) = F(3,1)*mi;
    K (7, 2) = 0;
    K (8, 2) = 0;
    K (9, 2) = 0;
    K (3, 3) = F(2,2)*lambda;
    K (4, 3) = 0;
    K (5, 3) = F(3,1)*lambda;
    K (6, 3) = F(2,1)*mi;
    K (7, 3) = 0;
    K (8, 3) = 0;
    K (9, 3) = 0;
    K (4, 4) = F(2,2)*lambda+2*F(2,2)*mi;
    K (5, 4) = F(1,2)*lambda+2*F(1,2)*mi;
    K (6, 4) = 0;
    K (7, 4) = F(2,3)*mi;
    K (8, 4) = F(1,3)*mi;
    K (9, 4) = 0;
    K (5, 5) = 3*F(2,2)*lambda+6*F(2,2)*mi;
    K (6, 5) = F(3,2)*lambda+2*F(3,2)*mi;
    K (7, 5) = F(1,3)*lambda;
    K (8, 5) = F(2,3)*lambda+2*F(2,3)*mi;
    K (9, 5) = F(3,3)*lambda;
    K (6, 6) = F(2,2)*lambda+2*F(2,2)*mi;
    K (7, 6) = 0;
    K (8, 6) = F(3,3)*mi;
    K (9, 6) = F(2,3)*mi;
    K (7, 7) = F(2,2)*lambda;
    K (8, 7) = F(1,2)*mi;
    K (9, 7) = 0;
    K (8, 8) = F(2,2)*lambda+2*F(2,2)*mi;
    K (9, 8) = F(3,2)*mi;
    K (9, 9) = F(2,2)*lambda;
    break;
  case 5:
    K (1, 1) = F(3,2)*lambda;
    K (2, 1) = 0;
    K (3, 1) = F(1,2)*mi;
    K (4, 1) = F(3,1)*mi;
    K (5, 1) = 0;
    K (6, 1) = F(1,1)*lambda;
    K (7, 1) = 0;
    K (8, 1) = 0;
    K (9, 1) = 0;
    K (2, 2) = F(3,2)*lambda;
    K (3, 2) = F(2,2)*mi;
    K (4, 2) = 0;
    K (5, 2) = F(3,1)*mi;
    K (6, 2) = F(2,1)*lambda;
    K (7, 2) = 0;
    K (8, 2) = 0;
    K (9, 2) = 0;
    K (3, 3) = F(3,2)*lambda+2*F(3,2)*mi;
    K (4, 3) = F(1,1)*mi;
    K (5, 3) = F(2,1)*mi;
    K (6, 3) = F(3,1)*lambda+2*F(3,1)*mi;
    K (7, 3) = 0;
    K (8, 3) = 0;
    K (9, 3) = 0;
    K (4, 4) = F(3,2)*lambda+2*F(3,2)*mi;
    K (5, 4) = 0;
    K (6, 4) = F(1,2)*lambda+2*F(1,2)*mi;
    K (7, 4) = F(3,3)*mi;
    K (8, 4) = 0;
    K (9, 4) = F(1,3)*mi;
    K (5, 5) = F(3,2)*lambda+2*F(3,2)*mi;
    K (6, 5) = F(2,2)*lambda+2*F(2,2)*mi;
    K (7, 5) = 0;
    K (8, 5) = F(3,3)*mi;
    K (9, 5) = F(2,3)*mi;
    K (6, 6) = 3*F(3,2)*lambda+6*F(3,2)*mi;
    K (7, 6) = F(1,3)*lambda;
    K (8, 6) = F(2,3)*lambda;
    K (9, 6) = F(3,3)*lambda+2*F(3,3)*mi;
    K (7, 7) = F(3,2)*lambda;
    K (8, 7) = 0;
    K (9, 7) = F(1,2)*mi;
    K (8, 8) = F(3,2)*lambda;
    K (9, 8) = F(2,2)*mi;
    K (9, 9) = F(3,2)*lambda+2*F(3,2)*mi;
    break;
  case 6:
    K (1, 1) = F(1,3)*lambda+2*F(1,3)*mi;
    K (2, 1) = F(2,3)*mi;
    K (3, 1) = F(3,3)*mi;
    K (4, 1) = 0;
    K (5, 1) = 0;
    K (6, 1) = 0;
    K (7, 1) = F(1,1)*lambda+2*F(1,1)*mi;
    K (8, 1) = F(2,1)*mi;
    K (9, 1) = F(3,1)*mi;
    K (2, 2) = F(1,3)*lambda;
    K (3, 2) = 0;
    K (4, 2) = 0;
    K (5, 2) = 0;
    K (6, 2) = 0;
    K (7, 2) = F(2,1)*lambda;
    K (8, 2) = F(1,1)*mi;
    K (9, 2) = 0;
    K (3, 3) = F(1,3)*lambda;
    K (4, 3) = 0;
    K (5, 3) = 0;
    K (6, 3) = 0;
    K (7, 3) = F(3,1)*lambda;
    K (8, 3) = 0;
    K (9, 3) = F(1,1)*mi;
    K (4, 4) = F(1,3)*lambda+2*F(1,3)*mi;
    K (5, 4) = F(2,3)*mi;
    K (6, 4) = F(3,3)*mi;
    K (7, 4) = F(1,2)*lambda+2*F(1,2)*mi;
    K (8, 4) = F(2,2)*mi;
    K (9, 4) = F(3,2)*mi;
    K (5, 5) = F(1,3)*lambda;
    K (6, 5) = 0;
    K (7, 5) = F(2,2)*lambda;
    K (8, 5) = F(1,2)*mi;
    K (9, 5) = 0;
    K (6, 6) = F(1,3)*lambda;
    K (7, 6) = F(3,2)*lambda;
    K (8, 6) = 0;
    K (9, 6) = F(1,2)*mi;
    K (7, 7) = 3*F(1,3)*lambda+6*F(1,3)*mi;
    K (8, 7) = F(2,3)*lambda+2*F(2,3)*mi;
    K (9, 7) = F(3,3)*lambda+2*F(3,3)*mi;
    K (8, 8) = F(1,3)*lambda+2*F(1,3)*mi;
    K (9, 8) = 0;
    K (9, 9) = F(1,3)*lambda+2*F(1,3)*mi;
    break;
  case 7:
    K (1, 1) = F(2,3)*lambda;
    K (2, 1) = F(1,3)*mi;
    K (3, 1) = 0;
    K (4, 1) = 0;
    K (5, 1) = 0;
    K (6, 1) = 0;
    K (7, 1) = F(2,1)*mi;
    K (8, 1) = F(1,1)*lambda;
    K (9, 1) = 0;
    K (2, 2) = F(2,3)*lambda+2*F(2,3)*mi;
    K (3, 2) = F(3,3)*mi;
    K (4, 2) = 0;
    K (5, 2) = 0;
    K (6, 2) = 0;
    K (7, 2) = F(1,1)*mi;
    K (8, 2) = F(2,1)*lambda+2*F(2,1)*mi;
    K (9, 2) = F(3,1)*mi;
    K (3, 3) = F(2,3)*lambda;
    K (4, 3) = 0;
    K (5, 3) = 0;
    K (6, 3) = 0;
    K (7, 3) = 0;
    K (8, 3) = F(3,1)*lambda;
    K (9, 3) = F(2,1)*mi;
    K (4, 4) = F(2,3)*lambda;
    K (5, 4) = F(1,3)*mi;
    K (6, 4) = 0;
    K (7, 4) = F(2,2)*mi;
    K (8, 4) = F(1,2)*lambda;
    K (9, 4) = 0;
    K (5, 5) = F(2,3)*lambda+2*F(2,3)*mi;
    K (6, 5) = F(3,3)*mi;
    K (7, 5) = F(1,2)*mi;
    K (8, 5) = F(2,2)*lambda+2*F(2,2)*mi;
    K (9, 5) = F(3,2)*mi;
    K (6, 6) = F(2,3)*lambda;
    K (7, 6) = 0;
    K (8, 6) = F(3,2)*lambda;
    K (9, 6) = F(2,2)*mi;
    K (7, 7) = F(2,3)*lambda+2*F(2,3)*mi;
    K (8, 7) = F(1,3)*lambda+2*F(1,3)*mi;
    K (9, 7) = 0;
    K (8, 8) = 3*F(2,3)*lambda+6*F(2,3)*mi;
    K (9, 8) = F(3,3)*lambda+2*F(3,3)*mi;
    K (9, 9) = F(2,3)*lambda+2*F(2,3)*mi;
    break;
  case 8:
    K (1, 1) = F(3,3)*lambda;
    K (2, 1) = 0;
    K (3, 1) = F(1,3)*mi;
    K (4, 1) = 0;
    K (5, 1) = 0;
    K (6, 1) = 0;
    K (7, 1) = F(3,1)*mi;
    K (8, 1) = 0;
    K (9, 1) = F(1,1)*lambda;
    K (2, 2) = F(3,3)*lambda;
    K (3, 2) = F(2,3)*mi;
    K (4, 2) = 0;
    K (5, 2) = 0;
    K (6, 2) = 0;
    K (7, 2) = 0;
    K (8, 2) = F(3,1)*mi;
    K (9, 2) = F(2,1)*lambda;
    K (3, 3) = F(3,3)*lambda+2*F(3,3)*mi;
    K (4, 3) = 0;
    K (5, 3) = 0;
    K (6, 3) = 0;
    K (7, 3) = F(1,1)*mi;
    K (8, 3) = F(2,1)*mi;
    K (9, 3) = F(3,1)*lambda+2*F(3,1)*mi;
    K (4, 4) = F(3,3)*lambda;
    K (5, 4) = 0;
    K (6, 4) = F(1,3)*mi;
    K (7, 4) = F(3,2)*mi;
    K (8, 4) = 0;
    K (9, 4) = F(1,2)*lambda;
    K (5, 5) = F(3,3)*lambda;
    K (6, 5) = F(2,3)*mi;
    K (7, 5) = 0;
    K (8, 5) = F(3,2)*mi;
    K (9, 5) = F(2,2)*lambda;
    K (6, 6) = F(3,3)*lambda+2*F(3,3)*mi;
    K (7, 6) = F(1,2)*mi;
    K (8, 6) = F(2,2)*mi;
    K (9, 6) = F(3,2)*lambda+2*F(3,2)*mi;
    K (7, 7) = F(3,3)*lambda+2*F(3,3)*mi;
    K (8, 7) = 0;
    K (9, 7) = F(1,3)*lambda+2*F(1,3)*mi;
    K (8, 8) = F(3,3)*lambda+2*F(3,3)*mi;
    K (9, 8) = F(2,3)*lambda+2*F(2,3)*mi;
    K (9, 9) = 3*F(3,3)*lambda+6*F(3,3)*mi;
    break;
  }

  K (1, 2) = K(2,1);
  K (1, 3) = K(3,1);
  K (2, 3) = K(3,2);
  K (1, 4) = K(4,1);
  K (2, 4) = K(4,2);
  K (3, 4) = K(4,3);
  K (1, 5) = K(5,1);
  K (2, 5) = K(5,2);
  K (3, 5) = K(5,3);
  K (4, 5) = K(5,4);
  K (1, 6) = K(6,1);
  K (2, 6) = K(6,2);
  K (3, 6) = K(6,3);
  K (4, 6) = K(6,4);
  K (5, 6) = K(6,5);
  K (1, 7) = K(7,1);
  K (2, 7) = K(7,2);
  K (3, 7) = K(7,3);
  K (4, 7) = K(7,4);
  K (5, 7) = K(7,5);
  K (6, 7) = K(7,6);
  K (1, 8) = K(8,1);
  K (2, 8) = K(8,2);
  K (3, 8) = K(8,3);
  K (4, 8) = K(8,4);
  K (5, 8) = K(8,5);
  K (6, 8) = K(8,6);
  K (7, 8) = K(8,7);
  K (1, 9) = K(9,1);
  K (2, 9) = K(9,2);
  K (3, 9) = K(9,3);
  K (4, 9) = K(9,4);
  K (5, 9) = K(9,5);
  K (6, 9) = K(9,6);
  K (7, 9) = K(9,7);
  K (8, 9) = K(9,8);

  /* volume scaling */
  if (volume != 1.0)
  {
    for (int col = 1; col <= 9; col ++)
    {
      double *val = & K(1, col); /* unroll */
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
    }
  }
}

/* column-wise F => .................................. */

#undef F
#undef K
#define F(I, J) F[(I-1) + 3 * (J-1)]
#define K(I, J) K[(I-1) + dim * (J-1)]
/* things below were generated with Maxima's 'itensor' package,
 * there is not much to comment on apart from that */
void SVK_Tangent_C (double lambda, double mi, double volume, int dim, double *F, double *K)
{
  K (1, 1) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+3*SQ(F(1,1))-3)*lambda + 
             (2*SQ(F(3,1))+2*SQ(F(2,1))+2*SQ(F(1,3))+2*SQ(F(1,2))+6*SQ(F(1,1))-2)*mi)*0.5;
  K (2, 1) = F(1,1)*F(2,1)*lambda+(F(1,3)*F(2,3)+F(1,2)*F(2,2)+2*F(1,1)*F(2,1))*mi;
  K (3, 1) = F(1,1)*F(3,1)*lambda+(F(1,3)*F(3,3)+F(1,2)*F(3,2)+2*F(1,1)*F(3,1))*mi;
  K (4, 1) = F(1,1)*F(1,2)*lambda+(F(3,1)*F(3,2)+F(2,1)*F(2,2)+2*F(1,1)*F(1,2))*mi;
  K (5, 1) = F(1,1)*F(2,2)*lambda+F(1,2)*F(2,1)*mi;
  K (6, 1) = F(1,1)*F(3,2)*lambda+F(1,2)*F(3,1)*mi;
  K (7, 1) = F(1,1)*F(1,3)*lambda+(F(3,1)*F(3,3)+F(2,1)*F(2,3)+2*F(1,1)*F(1,3))*mi;
  K (8, 1) = F(1,1)*F(2,3)*lambda+F(1,3)*F(2,1)*mi;
  K (9, 1) = F(1,1)*F(3,3)*lambda+F(1,3)*F(3,1)*mi;
  K (1, 2) = K(2,1);
  K (2, 2) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+3*SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,1))+2*SQ(F(2,3))+2*SQ(F(2,2))+6*SQ(F(2,1))+2*SQ(F(1,1))-2)*mi)*0.5;
  K (3, 2) = F(2,1)*F(3,1)*lambda+(F(2,3)*F(3,3)+F(2,2)*F(3,2)+2*F(2,1)*F(3,1))*mi;
  K (4, 2) = F(1,2)*F(2,1)*lambda+F(1,1)*F(2,2)*mi;
  K (5, 2) = F(2,1)*F(2,2)*lambda+(F(3,1)*F(3,2)+2*F(2,1)*F(2,2)+F(1,1)*F(1,2))*mi;
  K (6, 2) = F(2,1)*F(3,2)*lambda+F(2,2)*F(3,1)*mi;
  K (7, 2) = F(1,3)*F(2,1)*lambda+F(1,1)*F(2,3)*mi;
  K (8, 2) = F(2,1)*F(2,3)*lambda+(F(3,1)*F(3,3)+2*F(2,1)*F(2,3)+F(1,1)*F(1,3))*mi;
  K (9, 2) = F(2,1)*F(3,3)*lambda+F(2,3)*F(3,1)*mi;
  K (1, 3) = K(3,1);
  K (2, 3) = K(3,2);
  K (3, 3) = ((SQ(F(3,3))+SQ(F(3,2))+3*SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,3))+2*SQ(F(3,2))+6*SQ(F(3,1))+2*SQ(F(2,1))+2*SQ(F(1,1))-2)*mi)*0.5;
  K (4, 3) = F(1,2)*F(3,1)*lambda+F(1,1)*F(3,2)*mi;
  K (5, 3) = F(2,2)*F(3,1)*lambda+F(2,1)*F(3,2)*mi;
  K (6, 3) = F(3,1)*F(3,2)*lambda+(2*F(3,1)*F(3,2)+F(2,1)*F(2,2)+F(1,1)*F(1,2))*mi;
  K (7, 3) = F(1,3)*F(3,1)*lambda+F(1,1)*F(3,3)*mi;
  K (8, 3) = F(2,3)*F(3,1)*lambda+F(2,1)*F(3,3)*mi;
  K (9, 3) = F(3,1)*F(3,3)*lambda+(2*F(3,1)*F(3,3)+F(2,1)*F(2,3)+F(1,1)*F(1,3))*mi;
  K (1, 4) = K(4,1);
  K (2, 4) = K(4,2);
  K (3, 4) = K(4,3);
  K (4, 4) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+3*SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,2))+2*SQ(F(2,2))+2*SQ(F(1,3))+6*SQ(F(1,2))+2*SQ(F(1,1))-2)*mi)*0.5;
  K (5, 4) = F(1,2)*F(2,2)*lambda+(F(1,3)*F(2,3)+2*F(1,2)*F(2,2)+F(1,1)*F(2,1))*mi;
  K (6, 4) = F(1,2)*F(3,2)*lambda+(F(1,3)*F(3,3)+2*F(1,2)*F(3,2)+F(1,1)*F(3,1))*mi;
  K (7, 4) = F(1,2)*F(1,3)*lambda+(F(3,2)*F(3,3)+F(2,2)*F(2,3)+2*F(1,2)*F(1,3))*mi;
  K (8, 4) = F(1,2)*F(2,3)*lambda+F(1,3)*F(2,2)*mi;
  K (9, 4) = F(1,2)*F(3,3)*lambda+F(1,3)*F(3,2)*mi;
  K (1, 5) = K(5,1);
  K (2, 5) = K(5,2);
  K (3, 5) = K(5,3);
  K (4, 5) = K(5,4);
  K (5, 5) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+3*SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,2))+2*SQ(F(2,3))+6*SQ(F(2,2))+2*SQ(F(2,1))+2*SQ(F(1,2))-2)*mi)*0.5;
  K (6, 5) = F(2,2)*F(3,2)*lambda+(F(2,3)*F(3,3)+2*F(2,2)*F(3,2)+F(2,1)*F(3,1))*mi;
  K (7, 5) = F(1,3)*F(2,2)*lambda+F(1,2)*F(2,3)*mi;
  K (8, 5) = F(2,2)*F(2,3)*lambda+(F(3,2)*F(3,3)+2*F(2,2)*F(2,3)+F(1,2)*F(1,3))*mi;
  K (9, 5) = F(2,2)*F(3,3)*lambda+F(2,3)*F(3,2)*mi;
  K (1, 6) = K(6,1);
  K (2, 6) = K(6,2);
  K (3, 6) = K(6,3);
  K (4, 6) = K(6,4);
  K (5, 6) = K(6,5);
  K (6, 6) = ((SQ(F(3,3))+3*SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda + 
             (2*SQ(F(3,3))+6*SQ(F(3,2))+2*SQ(F(3,1))+2*SQ(F(2,2))+2*SQ(F(1,2))-2)*mi)*0.5;
  K (7, 6) = F(1,3)*F(3,2)*lambda+F(1,2)*F(3,3)*mi;
  K (8, 6) = F(2,3)*F(3,2)*lambda+F(2,2)*F(3,3)*mi;
  K (9, 6) = F(3,2)*F(3,3)*lambda+(2*F(3,2)*F(3,3)+F(2,2)*F(2,3)+F(1,2)*F(1,3))*mi;
  K (1, 7) = K(7,1);
  K (2, 7) = K(7,2);
  K (3, 7) = K(7,3);
  K (4, 7) = K(7,4);
  K (5, 7) = K(7,5);
  K (6, 7) = K(7,6);
  K (7, 7) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+3*SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +  
             (2*SQ(F(3,3))+2*SQ(F(2,3))+6*SQ(F(1,3))+2*SQ(F(1,2))+2*SQ(F(1,1))-2)*mi)*0.5;
  K (8, 7) = F(1,3)*F(2,3)*lambda+(2*F(1,3)*F(2,3)+F(1,2)*F(2,2)+F(1,1)*F(2,1))*mi;
  K (9, 7) = F(1,3)*F(3,3)*lambda+(2*F(1,3)*F(3,3)+F(1,2)*F(3,2)+F(1,1)*F(3,1))*mi;
  K (1, 8) = K(8,1);
  K (2, 8) = K(8,2);
  K (3, 8) = K(8,3);
  K (4, 8) = K(8,4);
  K (5, 8) = K(8,5);
  K (6, 8) = K(8,6);
  K (7, 8) = K(8,7);
  K (8, 8) = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+3*SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +   
             (2*SQ(F(3,3))+6*SQ(F(2,3))+2*SQ(F(2,2))+2*SQ(F(2,1))+2*SQ(F(1,3))-2)*mi)*0.5;
  K (9, 8) = F(2,3)*F(3,3)*lambda+(2*F(2,3)*F(3,3)+F(2,2)*F(3,2)+F(2,1)*F(3,1))*mi;
  K (1, 9) = K(9,1);
  K (2, 9) = K(9,2);
  K (3, 9) = K(9,3);
  K (4, 9) = K(9,4);
  K (5, 9) = K(9,5);
  K (6, 9) = K(9,6);
  K (7, 9) = K(9,7);
  K (8, 9) = K(9,8);
  K (9, 9) = ((3*SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +    
             (6*SQ(F(3,3))+2*SQ(F(3,2))+2*SQ(F(3,1))+2*SQ(F(2,3))+2*SQ(F(1,3))-2)*mi)*0.5;

  /* volume scaling */
  if (volume != 1.0)
  {
    for (int col = 1; col <= 9; col ++)
    {
      double *val = & K(1, col); /* unroll */
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
    }
  }
}

void SVK_Tangent_Diagonal_C (double lambda, double mi, double volume, int dim, double *F, double *K)
{
  K [0] = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+3*SQ(F(1,1))-3)*lambda + 
             (2*SQ(F(3,1))+2*SQ(F(2,1))+2*SQ(F(1,3))+2*SQ(F(1,2))+6*SQ(F(1,1))-2)*mi)*0.5;
  K [1] = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+3*SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,1))+2*SQ(F(2,3))+2*SQ(F(2,2))+6*SQ(F(2,1))+2*SQ(F(1,1))-2)*mi)*0.5;
  K [2] = ((SQ(F(3,3))+SQ(F(3,2))+3*SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,3))+2*SQ(F(3,2))+6*SQ(F(3,1))+2*SQ(F(2,1))+2*SQ(F(1,1))-2)*mi)*0.5;
  K [3] = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+3*SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,2))+2*SQ(F(2,2))+2*SQ(F(1,3))+6*SQ(F(1,2))+2*SQ(F(1,1))-2)*mi)*0.5;
  K [4] = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+3*SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +
             (2*SQ(F(3,2))+2*SQ(F(2,3))+6*SQ(F(2,2))+2*SQ(F(2,1))+2*SQ(F(1,2))-2)*mi)*0.5;
  K [5] = ((SQ(F(3,3))+3*SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda + 
             (2*SQ(F(3,3))+6*SQ(F(3,2))+2*SQ(F(3,1))+2*SQ(F(2,2))+2*SQ(F(1,2))-2)*mi)*0.5;
  K [6] = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+3*SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +  
             (2*SQ(F(3,3))+2*SQ(F(2,3))+6*SQ(F(1,3))+2*SQ(F(1,2))+2*SQ(F(1,1))-2)*mi)*0.5;
  K [7] = ((SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+3*SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +   
             (2*SQ(F(3,3))+6*SQ(F(2,3))+2*SQ(F(2,2))+2*SQ(F(2,1))+2*SQ(F(1,3))-2)*mi)*0.5;
  K [8] = ((3*SQ(F(3,3))+SQ(F(3,2))+SQ(F(3,1))+SQ(F(2,3))+SQ(F(2,2))+SQ(F(2,1))+SQ(F(1,3))+SQ(F(1,2))+SQ(F(1,1))-3)*lambda +    
             (6*SQ(F(3,3))+2*SQ(F(3,2))+2*SQ(F(3,1))+2*SQ(F(2,3))+2*SQ(F(1,3))-2)*mi)*0.5;

  if (volume != 1.0)
  {
    K [0] *= volume;
    K [1] *= volume;
    K [2] *= volume;
    K [3] *= volume;
    K [4] *= volume;
    K [5] *= volume;
    K [6] *= volume;
    K [7] *= volume;
    K [8] *= volume;
  }
}

void SVK_Tangent_Derivative_C (int comp,
 double lambda, double mi, double volume, int dim, double *F, double *K)
{
  switch (comp)
  {
  case 0:
    K (1, 1) = 3*F(1,1)*lambda+6*F(1,1)*mi;
    K (2, 1) = F(2,1)*lambda+2*F(2,1)*mi;
    K (3, 1) = F(3,1)*lambda+2*F(3,1)*mi;
    K (4, 1) = F(1,2)*lambda+2*F(1,2)*mi;
    K (5, 1) = F(2,2)*lambda;
    K (6, 1) = F(3,2)*lambda;
    K (7, 1) = F(1,3)*lambda+2*F(1,3)*mi;
    K (8, 1) = F(2,3)*lambda;
    K (9, 1) = F(3,3)*lambda;
    K (2, 2) = F(1,1)*lambda+2*F(1,1)*mi;
    K (3, 2) = 0;
    K (4, 2) = F(2,2)*mi;
    K (5, 2) = F(1,2)*mi;
    K (6, 2) = 0;
    K (7, 2) = F(2,3)*mi;
    K (8, 2) = F(1,3)*mi;
    K (9, 2) = 0;
    K (3, 3) = F(1,1)*lambda+2*F(1,1)*mi;
    K (4, 3) = F(3,2)*mi;
    K (5, 3) = 0;
    K (6, 3) = F(1,2)*mi;
    K (7, 3) = F(3,3)*mi;
    K (8, 3) = 0;
    K (9, 3) = F(1,3)*mi;
    K (4, 4) = F(1,1)*lambda+2*F(1,1)*mi;
    K (5, 4) = F(2,1)*mi;
    K (6, 4) = F(3,1)*mi;
    K (7, 4) = 0;
    K (8, 4) = 0;
    K (9, 4) = 0;
    K (5, 5) = F(1,1)*lambda;
    K (6, 5) = 0;
    K (7, 5) = 0;
    K (8, 5) = 0;
    K (9, 5) = 0;
    K (6, 6) = F(1,1)*lambda;
    K (7, 6) = 0;
    K (8, 6) = 0;
    K (9, 6) = 0;
    K (7, 7) = F(1,1)*lambda+2*F(1,1)*mi;
    K (8, 7) = F(2,1)*mi;
    K (9, 7) = F(3,1)*mi;
    K (8, 8) = F(1,1)*lambda;
    K (9, 8) = 0;
    K (9, 9) = F(1,1)*lambda;
    break;
  case 1:
    K (1, 1) = F(2,1)*lambda+2*F(2,1)*mi;
    K (2, 1) = F(1,1)*lambda+2*F(1,1)*mi;
    K (3, 1) = 0;
    K (4, 1) = F(2,2)*mi;
    K (5, 1) = F(1,2)*mi;
    K (6, 1) = 0;
    K (7, 1) = F(2,3)*mi;
    K (8, 1) = F(1,3)*mi;
    K (9, 1) = 0;
    K (2, 2) = 3*F(2,1)*lambda+6*F(2,1)*mi;
    K (3, 2) = F(3,1)*lambda+2*F(3,1)*mi;
    K (4, 2) = F(1,2)*lambda;
    K (5, 2) = F(2,2)*lambda+2*F(2,2)*mi;
    K (6, 2) = F(3,2)*lambda;
    K (7, 2) = F(1,3)*lambda;
    K (8, 2) = F(2,3)*lambda+2*F(2,3)*mi;
    K (9, 2) = F(3,3)*lambda;
    K (3, 3) = F(2,1)*lambda+2*F(2,1)*mi;
    K (4, 3) = 0 ;
    K (5, 3) = F(3,2)*mi;
    K (6, 3) = F(2,2)*mi;
    K (7, 3) = 0;
    K (8, 3) = F(3,3)*mi;
    K (9, 3) = F(2,3)*mi;
    K (4, 4) = F(2,1)*lambda;
    K (5, 4) = F(1,1)*mi;
    K (6, 4) = 0;
    K (7, 4) = 0;
    K (8, 4) = 0;
    K (9, 4) = 0;
    K (5, 5) = F(2,1)*lambda+2*F(2,1)*mi;
    K (6, 5) = F(3,1)*mi;
    K (7, 5) = 0;
    K (8, 5) = 0;
    K (9, 5) = 0;
    K (6, 6) = F(2,1)*lambda;
    K (7, 6) = 0;
    K (8, 6) = 0;
    K (9, 6) = 0;
    K (7, 7) = F(2,1)*lambda;
    K (8, 7) = F(1,1)*mi;
    K (9, 7) = 0;
    K (8, 8) = F(2,1)*lambda+2*F(2,1)*mi;
    K (9, 8) = F(3,1)*mi ;
    K (9, 9) = F(2,1)*lambda ;
    break;
  case 2:
    K (1, 1) = F(3,1)*lambda+2*F(3,1)*mi;
    K (2, 1) = 0;
    K (3, 1) = F(1,1)*lambda+2*F(1,1)*mi;
    K (4, 1) = F(3,2)*mi;
    K (5, 1) = 0;
    K (6, 1) = F(1,2)*mi;
    K (7, 1) = F(3,3)*mi;
    K (8, 1) = 0;
    K (9, 1) = F(1,3)*mi;
    K (2, 2) = F(3,1)*lambda+2*F(3,1)*mi;
    K (3, 2) = F(2,1)*lambda+2*F(2,1)*mi;
    K (4, 2) = 0;
    K (5, 2) = F(3,2)*mi;
    K (6, 2) = F(2,2)*mi;
    K (7, 2) = 0;
    K (8, 2) = F(3,3)*mi;
    K (9, 2) = F(2,3)*mi;
    K (3, 3) = 3*F(3,1)*lambda+6*F(3,1)*mi;
    K (4, 3) = F(1,2)*lambda;
    K (5, 3) = F(2,2)*lambda;
    K (6, 3) = F(3,2)*lambda+2*F(3,2)*mi;
    K (7, 3) = F(1,3)*lambda;
    K (8, 3) = F(2,3)*lambda;
    K (9, 3) = F(3,3)*lambda+2*F(3,3)*mi;
    K (4, 4) = F(3,1)*lambda;
    K (5, 4) = 0;
    K (6, 4) = F(1,1)*mi;
    K (7, 4) = 0;
    K (8, 4) = 0;
    K (9, 4) = 0;
    K (5, 5) = F(3,1)*lambda;
    K (6, 5) = F(2,1)*mi;
    K (7, 5) = 0;
    K (8, 5) = 0;
    K (9, 5) = 0;
    K (6, 6) = F(3,1)*lambda+2*F(3,1)*mi;
    K (7, 6) = 0;
    K (8, 6) = 0;
    K (9, 6) = 0;
    K (7, 7) = F(3,1)*lambda;
    K (8, 7) = 0;
    K (9, 7) = F(1,1)*mi;
    K (8, 8) = F(3,1)*lambda;
    K (9, 8) = F(2,1)*mi;
    K (9, 9) = F(3,1)*lambda+2*F(3,1)*mi;
    break;
  case 3:
    K (1, 1) = F(1,2)*lambda+2*F(1,2)*mi;
    K (2, 1) = F(2,2)*mi;
    K (3, 1) = F(3,2)*mi;
    K (4, 1) = F(1,1)*lambda+2*F(1,1)*mi;
    K (5, 1) = F(2,1)*mi;
    K (6, 1) = F(3,1)*mi;
    K (7, 1) = 0;
    K (8, 1) = 0;
    K (9, 1) = 0;
    K (2, 2) = F(1,2)*lambda;
    K (3, 2) = 0;
    K (4, 2) = F(2,1)*lambda;
    K (5, 2) = F(1,1)*mi;
    K (6, 2) = 0;
    K (7, 2) = 0;
    K (8, 2) = 0;
    K (9, 2) = 0;
    K (3, 3) = F(1,2)*lambda;
    K (4, 3) = F(3,1)*lambda;
    K (5, 3) = 0;
    K (6, 3) = F(1,1)*mi;
    K (7, 3) = 0;
    K (8, 3) = 0;
    K (9, 3) = 0;
    K (4, 4) = 3*F(1,2)*lambda+6*F(1,2)*mi;
    K (5, 4) = F(2,2)*lambda+2*F(2,2)*mi;
    K (6, 4) = F(3,2)*lambda+2*F(3,2)*mi;
    K (7, 4) = F(1,3)*lambda+2*F(1,3)*mi;
    K (8, 4) = F(2,3)*lambda;
    K (9, 4) = F(3,3)*lambda;
    K (5, 5) = F(1,2)*lambda+2*F(1,2)*mi;
    K (6, 5) = 0;
    K (7, 5) = F(2,3)*mi;
    K (8, 5) = F(1,3)*mi;
    K (9, 5) = 0;
    K (6, 6) = F(1,2)*lambda+2*F(1,2)*mi;
    K (7, 6) = F(3,3)*mi;
    K (8, 6) = 0;
    K (9, 6) = F(1,3)*mi;
    K (7, 7) = F(1,2)*lambda+2*F(1,2)*mi;
    K (8, 7) = F(2,2)*mi;
    K (9, 7) = F(3,2)*mi;
    K (8, 8) = F(1,2)*lambda;
    K (9, 8) = 0;
    K (9, 9) = F(1,2)*lambda;
    break;
  case 4:
    K (1, 1) = F(2,2)*lambda;
    K (2, 1) = F(1,2)*mi;
    K (3, 1) = 0;
    K (4, 1) = F(2,1)*mi;
    K (5, 1) = F(1,1)*lambda;
    K (6, 1) = 0;
    K (7, 1) = 0;
    K (8, 1) = 0;
    K (9, 1) = 0;
    K (2, 2) = F(2,2)*lambda+2*F(2,2)*mi;
    K (3, 2) = F(3,2)*mi;
    K (4, 2) = F(1,1)*mi;
    K (5, 2) = F(2,1)*lambda+2*F(2,1)*mi;
    K (6, 2) = F(3,1)*mi;
    K (7, 2) = 0;
    K (8, 2) = 0;
    K (9, 2) = 0;
    K (3, 3) = F(2,2)*lambda;
    K (4, 3) = 0;
    K (5, 3) = F(3,1)*lambda;
    K (6, 3) = F(2,1)*mi;
    K (7, 3) = 0;
    K (8, 3) = 0;
    K (9, 3) = 0;
    K (4, 4) = F(2,2)*lambda+2*F(2,2)*mi;
    K (5, 4) = F(1,2)*lambda+2*F(1,2)*mi;
    K (6, 4) = 0;
    K (7, 4) = F(2,3)*mi;
    K (8, 4) = F(1,3)*mi;
    K (9, 4) = 0;
    K (5, 5) = 3*F(2,2)*lambda+6*F(2,2)*mi;
    K (6, 5) = F(3,2)*lambda+2*F(3,2)*mi;
    K (7, 5) = F(1,3)*lambda;
    K (8, 5) = F(2,3)*lambda+2*F(2,3)*mi;
    K (9, 5) = F(3,3)*lambda;
    K (6, 6) = F(2,2)*lambda+2*F(2,2)*mi;
    K (7, 6) = 0;
    K (8, 6) = F(3,3)*mi;
    K (9, 6) = F(2,3)*mi;
    K (7, 7) = F(2,2)*lambda;
    K (8, 7) = F(1,2)*mi;
    K (9, 7) = 0;
    K (8, 8) = F(2,2)*lambda+2*F(2,2)*mi;
    K (9, 8) = F(3,2)*mi;
    K (9, 9) = F(2,2)*lambda;
    break;
  case 5:
    K (1, 1) = F(3,2)*lambda;
    K (2, 1) = 0;
    K (3, 1) = F(1,2)*mi;
    K (4, 1) = F(3,1)*mi;
    K (5, 1) = 0;
    K (6, 1) = F(1,1)*lambda;
    K (7, 1) = 0;
    K (8, 1) = 0;
    K (9, 1) = 0;
    K (2, 2) = F(3,2)*lambda;
    K (3, 2) = F(2,2)*mi;
    K (4, 2) = 0;
    K (5, 2) = F(3,1)*mi;
    K (6, 2) = F(2,1)*lambda;
    K (7, 2) = 0;
    K (8, 2) = 0;
    K (9, 2) = 0;
    K (3, 3) = F(3,2)*lambda+2*F(3,2)*mi;
    K (4, 3) = F(1,1)*mi;
    K (5, 3) = F(2,1)*mi;
    K (6, 3) = F(3,1)*lambda+2*F(3,1)*mi;
    K (7, 3) = 0;
    K (8, 3) = 0;
    K (9, 3) = 0;
    K (4, 4) = F(3,2)*lambda+2*F(3,2)*mi;
    K (5, 4) = 0;
    K (6, 4) = F(1,2)*lambda+2*F(1,2)*mi;
    K (7, 4) = F(3,3)*mi;
    K (8, 4) = 0;
    K (9, 4) = F(1,3)*mi;
    K (5, 5) = F(3,2)*lambda+2*F(3,2)*mi;
    K (6, 5) = F(2,2)*lambda+2*F(2,2)*mi;
    K (7, 5) = 0;
    K (8, 5) = F(3,3)*mi;
    K (9, 5) = F(2,3)*mi;
    K (6, 6) = 3*F(3,2)*lambda+6*F(3,2)*mi;
    K (7, 6) = F(1,3)*lambda;
    K (8, 6) = F(2,3)*lambda;
    K (9, 6) = F(3,3)*lambda+2*F(3,3)*mi;
    K (7, 7) = F(3,2)*lambda;
    K (8, 7) = 0;
    K (9, 7) = F(1,2)*mi;
    K (8, 8) = F(3,2)*lambda;
    K (9, 8) = F(2,2)*mi;
    K (9, 9) = F(3,2)*lambda+2*F(3,2)*mi;
    break;
  case 6:
    K (1, 1) = F(1,3)*lambda+2*F(1,3)*mi;
    K (2, 1) = F(2,3)*mi;
    K (3, 1) = F(3,3)*mi;
    K (4, 1) = 0;
    K (5, 1) = 0;
    K (6, 1) = 0;
    K (7, 1) = F(1,1)*lambda+2*F(1,1)*mi;
    K (8, 1) = F(2,1)*mi;
    K (9, 1) = F(3,1)*mi;
    K (2, 2) = F(1,3)*lambda;
    K (3, 2) = 0;
    K (4, 2) = 0;
    K (5, 2) = 0;
    K (6, 2) = 0;
    K (7, 2) = F(2,1)*lambda;
    K (8, 2) = F(1,1)*mi;
    K (9, 2) = 0;
    K (3, 3) = F(1,3)*lambda;
    K (4, 3) = 0;
    K (5, 3) = 0;
    K (6, 3) = 0;
    K (7, 3) = F(3,1)*lambda;
    K (8, 3) = 0;
    K (9, 3) = F(1,1)*mi;
    K (4, 4) = F(1,3)*lambda+2*F(1,3)*mi;
    K (5, 4) = F(2,3)*mi;
    K (6, 4) = F(3,3)*mi;
    K (7, 4) = F(1,2)*lambda+2*F(1,2)*mi;
    K (8, 4) = F(2,2)*mi;
    K (9, 4) = F(3,2)*mi;
    K (5, 5) = F(1,3)*lambda;
    K (6, 5) = 0;
    K (7, 5) = F(2,2)*lambda;
    K (8, 5) = F(1,2)*mi;
    K (9, 5) = 0;
    K (6, 6) = F(1,3)*lambda;
    K (7, 6) = F(3,2)*lambda;
    K (8, 6) = 0;
    K (9, 6) = F(1,2)*mi;
    K (7, 7) = 3*F(1,3)*lambda+6*F(1,3)*mi;
    K (8, 7) = F(2,3)*lambda+2*F(2,3)*mi;
    K (9, 7) = F(3,3)*lambda+2*F(3,3)*mi;
    K (8, 8) = F(1,3)*lambda+2*F(1,3)*mi;
    K (9, 8) = 0;
    K (9, 9) = F(1,3)*lambda+2*F(1,3)*mi;
    break;
  case 7:
    K (1, 1) = F(2,3)*lambda;
    K (2, 1) = F(1,3)*mi;
    K (3, 1) = 0;
    K (4, 1) = 0;
    K (5, 1) = 0;
    K (6, 1) = 0;
    K (7, 1) = F(2,1)*mi;
    K (8, 1) = F(1,1)*lambda;
    K (9, 1) = 0;
    K (2, 2) = F(2,3)*lambda+2*F(2,3)*mi;
    K (3, 2) = F(3,3)*mi;
    K (4, 2) = 0;
    K (5, 2) = 0;
    K (6, 2) = 0;
    K (7, 2) = F(1,1)*mi;
    K (8, 2) = F(2,1)*lambda+2*F(2,1)*mi;
    K (9, 2) = F(3,1)*mi;
    K (3, 3) = F(2,3)*lambda;
    K (4, 3) = 0;
    K (5, 3) = 0;
    K (6, 3) = 0;
    K (7, 3) = 0;
    K (8, 3) = F(3,1)*lambda;
    K (9, 3) = F(2,1)*mi;
    K (4, 4) = F(2,3)*lambda;
    K (5, 4) = F(1,3)*mi;
    K (6, 4) = 0;
    K (7, 4) = F(2,2)*mi;
    K (8, 4) = F(1,2)*lambda;
    K (9, 4) = 0;
    K (5, 5) = F(2,3)*lambda+2*F(2,3)*mi;
    K (6, 5) = F(3,3)*mi;
    K (7, 5) = F(1,2)*mi;
    K (8, 5) = F(2,2)*lambda+2*F(2,2)*mi;
    K (9, 5) = F(3,2)*mi;
    K (6, 6) = F(2,3)*lambda;
    K (7, 6) = 0;
    K (8, 6) = F(3,2)*lambda;
    K (9, 6) = F(2,2)*mi;
    K (7, 7) = F(2,3)*lambda+2*F(2,3)*mi;
    K (8, 7) = F(1,3)*lambda+2*F(1,3)*mi;
    K (9, 7) = 0;
    K (8, 8) = 3*F(2,3)*lambda+6*F(2,3)*mi;
    K (9, 8) = F(3,3)*lambda+2*F(3,3)*mi;
    K (9, 9) = F(2,3)*lambda+2*F(2,3)*mi;
    break;
  case 8:
    K (1, 1) = F(3,3)*lambda;
    K (2, 1) = 0;
    K (3, 1) = F(1,3)*mi;
    K (4, 1) = 0;
    K (5, 1) = 0;
    K (6, 1) = 0;
    K (7, 1) = F(3,1)*mi;
    K (8, 1) = 0;
    K (9, 1) = F(1,1)*lambda;
    K (2, 2) = F(3,3)*lambda;
    K (3, 2) = F(2,3)*mi;
    K (4, 2) = 0;
    K (5, 2) = 0;
    K (6, 2) = 0;
    K (7, 2) = 0;
    K (8, 2) = F(3,1)*mi;
    K (9, 2) = F(2,1)*lambda;
    K (3, 3) = F(3,3)*lambda+2*F(3,3)*mi;
    K (4, 3) = 0;
    K (5, 3) = 0;
    K (6, 3) = 0;
    K (7, 3) = F(1,1)*mi;
    K (8, 3) = F(2,1)*mi;
    K (9, 3) = F(3,1)*lambda+2*F(3,1)*mi;
    K (4, 4) = F(3,3)*lambda;
    K (5, 4) = 0;
    K (6, 4) = F(1,3)*mi;
    K (7, 4) = F(3,2)*mi;
    K (8, 4) = 0;
    K (9, 4) = F(1,2)*lambda;
    K (5, 5) = F(3,3)*lambda;
    K (6, 5) = F(2,3)*mi;
    K (7, 5) = 0;
    K (8, 5) = F(3,2)*mi;
    K (9, 5) = F(2,2)*lambda;
    K (6, 6) = F(3,3)*lambda+2*F(3,3)*mi;
    K (7, 6) = F(1,2)*mi;
    K (8, 6) = F(2,2)*mi;
    K (9, 6) = F(3,2)*lambda+2*F(3,2)*mi;
    K (7, 7) = F(3,3)*lambda+2*F(3,3)*mi;
    K (8, 7) = 0;
    K (9, 7) = F(1,3)*lambda+2*F(1,3)*mi;
    K (8, 8) = F(3,3)*lambda+2*F(3,3)*mi;
    K (9, 8) = F(2,3)*lambda+2*F(2,3)*mi;
    K (9, 9) = 3*F(3,3)*lambda+6*F(3,3)*mi;
    break;
  }

  K (1, 2) = K(2,1);
  K (1, 3) = K(3,1);
  K (2, 3) = K(3,2);
  K (1, 4) = K(4,1);
  K (2, 4) = K(4,2);
  K (3, 4) = K(4,3);
  K (1, 5) = K(5,1);
  K (2, 5) = K(5,2);
  K (3, 5) = K(5,3);
  K (4, 5) = K(5,4);
  K (1, 6) = K(6,1);
  K (2, 6) = K(6,2);
  K (3, 6) = K(6,3);
  K (4, 6) = K(6,4);
  K (5, 6) = K(6,5);
  K (1, 7) = K(7,1);
  K (2, 7) = K(7,2);
  K (3, 7) = K(7,3);
  K (4, 7) = K(7,4);
  K (5, 7) = K(7,5);
  K (6, 7) = K(7,6);
  K (1, 8) = K(8,1);
  K (2, 8) = K(8,2);
  K (3, 8) = K(8,3);
  K (4, 8) = K(8,4);
  K (5, 8) = K(8,5);
  K (6, 8) = K(8,6);
  K (7, 8) = K(8,7);
  K (1, 9) = K(9,1);
  K (2, 9) = K(9,2);
  K (3, 9) = K(9,3);
  K (4, 9) = K(9,4);
  K (5, 9) = K(9,5);
  K (6, 9) = K(9,6);
  K (7, 9) = K(9,7);
  K (8, 9) = K(9,8);

  /* volume scaling */
  if (volume != 1.0)
  {
    for (int col = 1; col <= 9; col ++)
    {
      double *val = & K(1, col); /* unroll */
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
      (*val) *= volume; val ++;
    }
  }
}
