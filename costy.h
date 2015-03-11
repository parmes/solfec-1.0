#ifndef __costy__
#define __costy__

// refine convex edges in mesh (looking from outside) when generating TetGen tetrahedrel refined mesh
// returns the size of mtroutlist (allocated inside)
int refineEdgesTetgen (MESH *msh, double MIN_ANGLE, double MAX_ANGLE, double edgeLength, double* mtroutlist);

#endif