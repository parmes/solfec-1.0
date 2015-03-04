#ifndef __costy__
#define __costy__

// refine convex edges in mesh (looking from outside) when generating TetGen tetrahedrel dense mesh
// returns the size of mtrout (allocated inside)
int refineEdgesTetgen (MESH *msh, int gg, double edgeLength, double **mtrout);

#endif
