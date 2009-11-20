/*
 * cli.h
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * constraints clique
 */

/*
   Underlying ideas are simple:

   1. In a constraint graph, constraints adjacent to a single body form a clique.

   2. Move all local dynamics data to cliques and store cliques at bodies.

   3. In parallel, balance bodies together with cliques (let them affect weights).

   4. Implement clique based, decoupled (parallelizable) solver.

   5. Develop other solvers based on cliques.
*/

#ifndef __cli__
#define __cli__

typedef struct cliadj CLIADJ;
typedef struct clique CLIQUE;

struct cliadj
{
  double W [9]; /* off-diagonal block (&W[0] can be casted at CLIADJ*) */

  double *R; /* adjacent constraint reaction (can be casted at CLIQUE*) */

  BODY *bod;

  double *upper; /* if lower diagonal, points to the transposed upper diagonal block; NULL otherwise */
}

struct clique
{
  double R [3]; /* reaction (&R[0] can be casted at CLIQUE*) */

  double W [9]; /* diagonal block */

  CON *con;

  CLIADJ *adj; /* adjacent constraints (for every pair of constraints between same pair of bodies adj->R repeats twice with different body pointer) */

  int nadj, sadj;

  CLIQUE *p, *n;
};

/* insert constraint into body pair cliques */
void CLIQUE_Insert (CON *con, BODY *one, BODY *two);

/* remove constraint from its cliques */
void CLIQUE_Remove (CON *con);

/* update body clique */
void CLIQUE_Update (BODY *bod);

/* destroy body clique */
void CLIQUE_Destroy (BODY *bod);

#endif
