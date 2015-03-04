#include <cmath>

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

extern "C"
{
#include "sol.h"
#include "dom.h"
#include "msh.h"
}

#define PI 3.14159265358979323846

void dummy()
{
}

//______________________________________________________________________

class Vector_3 {

public:

	Vector_3() : x(0.0), y(0.0), z(0.0) {};
	Vector_3(double _x, double _y, double _z) : x(_x), y (_y), z(_z) {};

	void normalize() {
		double norm = std::sqrt(x*x + y*y + z*z);
		x = x / norm;
		y = y / norm;
		z = z / norm;
	}

	void neg() {
		x = -1.0 * x;
		y = -1.0 * y;
		z = -1.0 * z;
	}

	void zero() {
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}

	double x, y, z;
};


double inner(const Vector_3& a, const Vector_3& b) {

	return (a.x*b.x + a.y*b.y + a.z*b.z);

}


double norm(const Vector_3& a) {

	return std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);

}


void cross(const Vector_3& a, const Vector_3& b, Vector_3& c) {

	c.x = a.y*b.z - a.z*b.y;
	c.y = a.z*b.x - a.x*b.z;
	c.z = a.x*b.y - a.y*b.x;

}


// calculate angle between two vectors
double angle(const Vector_3& a, const Vector_3& b) {

  Vector_3 c;
  cross(a,b,c);

  return std::atan2( norm(c), inner(a,b) );

}


//______________________________________________________________________

typedef struct face Face;

// edge class
class Edge {

public:

	Edge() {};
	Edge(int VertsIndx_1, int VertsIndx_2, Face* f_1, Face* f_2) {

		VertsIndx[0] = VertsIndx_1;
		VertsIndx[1] = VertsIndx_2;
		f1 = f_1;
		f2 = f_2;

	}

	int VertsIndx[2];
	Face* f1;
	Face* f2;

};


// return whether face is on the left of an edge looking from outside
bool onLeft(Face* fac, Edge& edg) {

	int i = 0;

	for (; i < fac->type; ++i)
		if (edg.VertsIndx[0] == fac->nodes[i])
			break;

	if (i < (fac->type)-1)
		if (edg.VertsIndx[1] == fac->nodes[i+1]) 
			return true;

	if (i == (fac->type)-1)
		if (edg.VertsIndx[1] == fac->nodes[0]) 
			return true;

	return false;

}


// calculate dihedral angle, +ve for convex looking from outside
double dihedralAngle(Face* f_0, Face* f_1, Vector_3& p, bool left) {

	Vector_3 n0, n1, n;

	n0.x = f_0->normal[0]; n0.y = f_0->normal[1]; n0.z = f_0->normal[2];
	n1.x = f_1->normal[0]; n1.y = f_1->normal[1]; n1.z = f_1->normal[2];

	double theta = PI - angle(n0,n1);

	if (left)
	{
		cross(n1, n0, n);
		if (inner(n,p) > 0) return ( theta );
		else return ( -1*theta );
	} else {
		cross(n0, n1, n);
		if (inner(n,p) > 0) return ( theta );
		else return ( -1*theta );
	}

}


//______________________________________________________________________

// establish connectivity of an edge; searching for second face incident on an edge
void connectivity(std::vector<Edge>& edges, std::vector<Edge>& uqEdges) {

	std::map<int,int> vtx_edge; // map to store <vertex index,edge number in vector>
	std::pair<std::map<int,int>::iterator,bool> ret1, ret2;

	// loop through surface faces
	for (int i = 0; i < (signed)edges.size(); ++i)
	{
		// insert vertex indices into map
		ret1 = vtx_edge.insert( std::pair<int,int>(edges[i].VertsIndx[0],i) );
		ret2 = vtx_edge.insert( std::pair<int,int>(edges[i].VertsIndx[1],i) );

		// if both already exist & are therefore not inserted
		if (ret1.second == false && ret2.second == false)
		{
			// compare the edge numbers in the map
			if (ret1.first->second == ret2.first->second)
			{
				// store face on edge object
				edges[i].f2 = edges[ret1.first->second].f1;
				//edges[ret1.first->second].f2 = edges[i].f1;

				// store unique edge vector reference
				uqEdges.push_back( edges[i] );
			}
		}
	}

}


// find the convex edges in the mesh (g = 0,1,2)
void findEdgesToRefine(MESH* msh, std::vector<Edge>& uqEdges, int g, std::vector<Edge>& convexEdges) {

	double theta = 0.0;
	bool on_left = false;
	Vector_3 vec1, vec2, axis;

	// loop through edges vector
	for (std::vector<Edge>::iterator it = uqEdges.begin(); it != uqEdges.end(); ++it)
	{
		// determine axis direction
		vec1.x = msh->ref_nodes[ (*it).VertsIndx[0] ][0]; vec1.y = msh->ref_nodes[ (*it).VertsIndx[0] ][1]; vec1.z = msh->ref_nodes[ (*it).VertsIndx[0] ][2];
		vec2.x = msh->ref_nodes[ (*it).VertsIndx[1] ][0]; vec2.y = msh->ref_nodes[ (*it).VertsIndx[1] ][1]; vec2.z = msh->ref_nodes[ (*it).VertsIndx[1] ][2];

		if (msh->ref_nodes[ (*it).VertsIndx[0] ][g] > msh->ref_nodes[ (*it).VertsIndx[1] ][g])
		{
			axis.x = vec1.x - vec2.x;
			axis.y = vec1.y - vec2.y;
			axis.z = vec1.z - vec2.z;
		}
		else
		{
			axis.x = vec2.x - vec1.x;
			axis.y = vec2.y - vec1.y;
			axis.z = vec2.z - vec1.z;
		}
		axis.normalize();

		// establish if f1 is on left
		on_left = onLeft( (*it).f1, *it );
		theta = dihedralAngle( (*it).f1, (*it).f2, axis, on_left );

		if (theta >= 0 && theta < (179*(PI/180)) )
			convexEdges.push_back( *it );
	}

}


// output mtr file with mesh edge length constraint on vertices
void outputMTR(MESH* msh, std::vector<Edge>& convexEdges, double length, std::vector<double> mtrout) {

	std::set<int> uqVerts;

	// find unique vertex indices from vector
	for (std::vector<Edge>::iterator it = convexEdges.begin(); it != convexEdges.end(); ++it)
	{
		uqVerts.insert( (*it).VertsIndx[0] );
		uqVerts.insert( (*it).VertsIndx[1] );
	}

	//____________________________________
	// output file that is readable by TetGen

	for (int i = 0; i < msh->nodes_count; ++i)
	{
		if (uqVerts.count(i) > 0)
		        mtrout.push_back (length); 
		else 
		        mtrout.push_back (1.0); 
	}
}

//______________________________________________________________________

extern "C"
{
// refine convex edges in mesh (looking from outside) when generating TetGen tetrahedrel dense mesh
int refineEdgesTetgen (MESH *msh, int gg, double edgeLength, double **mtrout) {
	int k = 0;
	FACE* _face;
	Edge _edge;
	std::vector<Edge> _edges, _uqEdges, _convexEdges;

	MESH* _msh = MESH_Copy (msh);
        MESH_Update (_msh, NULL, NULL, NULL); /* reference configuration */

	// loop through surface faces
	for (_face = _msh->faces; _face; _face = _face->n)
	{
		// loop through edges of face & store in edge list
		k = _face->type;
		for (int i = 0; i < k-1; ++i)
			_edges.push_back( Edge(_face->nodes[i], _face->nodes[i+1], _face, NULL) );

		_edges.push_back( Edge(_face->nodes[k-1], _face->nodes[0], _face, NULL) );
	}

	connectivity(_edges, _uqEdges);
	findEdgesToRefine(_msh, _uqEdges, gg, _convexEdges);

	std::vector<double> _mtrout;
	outputMTR(_msh, _convexEdges, edgeLength, _mtrout);

	(*mtrout) = (double*) malloc (sizeof (double) *  _mtrout.size());
	memcpy ((*mtrout), &(_mtrout[0]), sizeof (double) * _mtrout.size());

	//std::string edgeMeshFname = "convex_edges.vtk";
	//edges2VTK(edgeMeshFname, _msh, _convexEdges);

	_edges.clear();
	_uqEdges.clear();
	_convexEdges.clear();
	MESH_Destroy (_msh);

	return _mtrout.size();
}
}

//______________________________________________________________________

// output edges in VTK file format
void edges2VTK(std::string& fname, MESH* msh, std::vector<Edge>& convexEdges) {

	int i = -1;
	std::vector<Vector_3> uqVerts;
	std::map<int,int> verts_indx;
	std::pair<std::map<int,int>::iterator,bool> ret;

	for (std::vector<Edge>::iterator it = convexEdges.begin(); it != convexEdges.end(); ++it)
	{
		ret = verts_indx.insert( std::pair<int,int>((*it).VertsIndx[0], ++i) );
		if (ret.second == true) 
			uqVerts.push_back( Vector_3(msh->ref_nodes[ (*it).VertsIndx[0] ][0], msh->ref_nodes[ (*it).VertsIndx[0] ][1], msh->ref_nodes[ (*it).VertsIndx[0] ][2]) );
		else --i;

		ret = verts_indx.insert( std::pair<int,int>((*it).VertsIndx[1], ++i) );
		if (ret.second == true) 
			uqVerts.push_back( Vector_3(msh->ref_nodes[ (*it).VertsIndx[1] ][0], msh->ref_nodes[ (*it).VertsIndx[1] ][1], msh->ref_nodes[ (*it).VertsIndx[1] ][2]) );
		else --i;
	}

	//_______________________________________________________

	std::ofstream fin;
	fin.open( fname.c_str() );

	// print out the header data
	fin << "# vtk DataFile Version 2.0\n";
	fin << "Convex edges\n";
	fin << "ASCII\n";

	//print out the coordinates of the vertices
	fin << "DATASET UNSTRUCTURED_GRID\n";
	fin << "POINTS " << uqVerts.size() << " float\n";

	for (std::vector<Vector_3>::iterator it= uqVerts.begin(); it != uqVerts.end(); ++it)  
		fin << (*it).x  << " " << (*it).y << " " << (*it).z << std::endl;

  //_______________________________________________________

	// print out the element index list
	fin << std::endl;
	fin << "CELLS " << convexEdges.size() << " " << 3*convexEdges.size() << std::endl; //(2+1)*no_of_edges

	// loop through edges vector
	for (std::vector<Edge>::iterator it = convexEdges.begin(); it != convexEdges.end(); ++it)
		fin << "2 " << verts_indx[ (*it).VertsIndx[0] ] << " " <<  verts_indx[ (*it).VertsIndx[1] ] << std::endl;

  //_______________________________________________________

	//print out VTK cell type
	fin << std::endl;
	fin << "CELL_TYPES " << convexEdges.size() << std::endl;
	
	for(int j = 0; j < (int)convexEdges.size(); ++j) 
		fin << "3\n";

	fin.close();
}
