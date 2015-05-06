#include <cmath>

#include <assert.h> 
#include <algorithm>
#include <vector>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <string.h>

extern "C"
{
#include "../sol.h"
#include "../dom.h"
#include "../msh.h"
}

#include "costy.h"

#define PI 3.14159265358979323846

typedef struct face Face;

//______________________________________________________________________

// Edge class
class Edge {

public:

	Edge() {};
	Edge(int VertIndx_1, int VertIndx_2, int _f1_VertIndx, Face* _fac1, Face* _fac2) {

		VertsIndx[0] = VertIndx_1;
		VertsIndx[1] = VertIndx_2;
		f1_VertIndx = _f1_VertIndx;
		fac1 = _fac1;
		fac2 = _fac2;
	}

	int VertsIndx[2];
	int f1_VertIndx, f2_VertIndx;
	Face* fac1;
	Face* fac2;

};


// find vertex index of second incident face on give edge
int findIndx(Face* _fac, Edge& edg) {

	bool found = false;

	// loop through vertex indices of face 
	int k = _fac->type;
	for (int i = 0; i < k; ++i)
	{
		// check if vertex node of face equals that belonging to edg.VertsIndx[1]
		if (_fac->nodes[i] == edg.VertsIndx[1])
		{
			if (i >= 0  && i < k-2) return _fac->nodes[i+2];
			if (i == k-2) return _fac->nodes[0];
			if (i == k-1) return _fac->nodes[1];

			found = true;
			break;
		}
	}

	assert(found);

	return 0;

}


// find convex edges and return in std::vector convexEdges; angles given in degrees
void findConvexEdges(MESH *msh, std::vector<Edge>& edges, double ANGLE_MIN, double ANGLE_MAX, std::vector<Edge>& convexEdges) {

	double vol = 0.0, theta = 0.0;
	Vector_3 dvw, dvx, dvy, cd, n1, n2;

	for (std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); ++it)
	{
		dvw.x = msh->ref_nodes[ (*it).VertsIndx[1] ][0] - msh->ref_nodes[ (*it).VertsIndx[0] ][0];
		dvw.y = msh->ref_nodes[ (*it).VertsIndx[1] ][1] - msh->ref_nodes[ (*it).VertsIndx[0] ][1];
		dvw.z = msh->ref_nodes[ (*it).VertsIndx[1] ][2] - msh->ref_nodes[ (*it).VertsIndx[0] ][2];

		dvx.x = msh->ref_nodes[ (*it).f1_VertIndx ][0] - msh->ref_nodes[ (*it).VertsIndx[0] ][0];
		dvx.y = msh->ref_nodes[ (*it).f1_VertIndx ][1] - msh->ref_nodes[ (*it).VertsIndx[0] ][1];
		dvx.z = msh->ref_nodes[ (*it).f1_VertIndx ][2] - msh->ref_nodes[ (*it).VertsIndx[0] ][2];

		dvy.x = msh->ref_nodes[ (*it).f2_VertIndx ][0] - msh->ref_nodes[ (*it).VertsIndx[0] ][0];
		dvy.y = msh->ref_nodes[ (*it).f2_VertIndx ][1] - msh->ref_nodes[ (*it).VertsIndx[0] ][1];
		dvy.z = msh->ref_nodes[ (*it).f2_VertIndx ][2] - msh->ref_nodes[ (*it).VertsIndx[0] ][2];

		cross(dvx,dvy,cd);
		vol = (1.0/6.0)*inner(dvw,cd);

		// if vol positive, the edge is convex
		if (vol > 0) 
		{
			// check angle between faces 
			n1.x = (*it).fac1->normal[0]; n1.y = (*it).fac1->normal[1]; n1.z = (*it).fac1->normal[2];
			n2.x = (*it).fac2->normal[0]; n2.y = (*it).fac2->normal[1]; n2.z = (*it).fac2->normal[2];
			theta = PI - angle(n1,n2);
			theta = theta*(180.0/PI); // convert to degrees

			// check if theta is within limits
			if (theta >= ANGLE_MIN && theta < ANGLE_MAX)
				convexEdges.push_back( *it );
		}
	}

}

//______________________________________________________________________

// output mtr file with mesh edge length constraint on vertices
void outputMTR(MESH* msh, std::vector<Edge>& convexEdges, double length, std::vector<double>& mtrout) {

	std::set<int> uqVerts;

	// find unique vertex indices from convexEdges
	for (std::vector<Edge>::iterator it = convexEdges.begin(); it != convexEdges.end(); ++it)
	{
		uqVerts.insert( (*it).VertsIndx[0] );
		uqVerts.insert( (*it).VertsIndx[1] );
	}

	// output file that is readable by TetGen
	for (int i = 0; i < msh->nodes_count; ++i)
	{
		if (uqVerts.count(i) > 0)
		    mtrout.push_back(length); 
		else 
			mtrout.push_back(0);
	}

}

//______________________________________________________________________

// output edges in VTK file format
void edges2VTK(std::string& fname, MESH* msh, std::vector<Edge>& edges) {

	int i = -1;
	std::vector<Vector_3> uqVerts;
	std::map<int,int> verts_indx;
	std::pair<std::map<int,int>::iterator,bool> ret;

	for (std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); ++it)
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
	fin << "Edges\n";
	fin << "ASCII\n";
	fin << "DATASET UNSTRUCTURED_GRID\n";
	fin << "POINTS " << uqVerts.size() << " float\n";

	for (std::vector<Vector_3>::iterator it= uqVerts.begin(); it != uqVerts.end(); ++it)  
		fin << (*it).x  << " " << (*it).y << " " << (*it).z << std::endl;

  //_______________________________________________________

	// print out the element index list
	fin << std::endl;
	fin << "CELLS " << edges.size() << " " << 3*edges.size() << std::endl; // (2+1)*no_of_edges

	// loop through edges std::vector
	for (std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); ++it)
		fin << "2 " << verts_indx[ (*it).VertsIndx[0] ] << " " <<  verts_indx[ (*it).VertsIndx[1] ] << std::endl;

  //_______________________________________________________

	//print out VTK cell type
	fin << std::endl;
	fin << "CELL_TYPES " << edges.size() << std::endl;
	
	for(i = 0; i < (int)edges.size(); ++i) 
		fin << "3\n";

	fin.close();

}

//______________________________________________________________________

// Export body reference mesh into VTK
void EXPORT_SOLFEC_MESH (std::string& fname, MESH* msh) {

	int n = 0, elno = 0, el_size = 0;
	ELEMENT *ele;

	std::ofstream fin;
	fin.open( fname.c_str() );

	// file output start
	fin << "# vtk DataFile Version 2.0\n";
	fin << "SOLFEC ORIGINAL BODY MESH\n";
	fin << "ASCII\n";
	fin << std::endl;
	fin << "DATASET UNSTRUCTURED_GRID\n";
	fin << "Points " << msh->nodes_count << " float\n";


	for (; n < msh->nodes_count; n++) 
		fin << msh->ref_nodes[n][0] << " " <<  msh->ref_nodes[n][1] << " " <<  msh->ref_nodes[n][2] << std::endl;

	//______________________________________________________
	elno = msh->surfeles_count + msh->bulkeles_count; // element no

	// calculate size
	for (ele = msh->surfeles; ele; ele = ele->next)
	{
		if (ele->type == 4) el_size += 5;
		else if (ele->type == 8) el_size += 9;
	}

	for (ele = msh->bulkeles; ele; ele = ele->next)
	{
		if (ele->type == 4) el_size += 5;
		else if (ele->type == 8) el_size += 9;
	}

	fin << std::endl;
	fin << "CELLS " << elno << " " << el_size << std::endl;

	for (ele = msh->surfeles; ele; ele = ele->next)
	{
		if (ele->type == 4) fin << "4 " << ele->nodes[0] << " " << ele->nodes[1] << " " << ele->nodes[2] << " " <<  ele->nodes[3] << std::endl;
		else if (ele->type == 8) fin << "8 " << ele->nodes[0] << " " << ele->nodes[1] << " " << ele->nodes[2] << " " <<  ele->nodes[3] << " " <<  ele->nodes[4] << " " <<  ele->nodes[5] << " " <<  ele->nodes[6] << " " <<  ele->nodes[7] << std::endl;
	}

	for (ele = msh->bulkeles; ele; ele = ele->next)
	{
		if (ele->type == 4) fin << "4 " << ele->nodes[0] << " " << ele->nodes[1] << " " << ele->nodes[2] << " " <<  ele->nodes[3] << std::endl;
		else if (ele->type == 8) fin << "8 " << ele->nodes[0] << " " << ele->nodes[1] << " " << ele->nodes[2] << " " <<  ele->nodes[3] << " " <<  ele->nodes[4] << " " <<  ele->nodes[5] << " " <<  ele->nodes[6] << " " <<  ele->nodes[7] << std::endl;
	}

	//______________________________________________________
	fin << std::endl;
	fin << "CELL_TYPES " << elno << std::endl;

	for (ele = msh->surfeles; ele; ele = ele->next)
	{
		if (ele->type == 4) fin << "10\n";
		else if (ele->type == 8) fin << "12\n";
	}

	for (ele = msh->bulkeles; ele; ele = ele->next)
	{
		if (ele->type == 4) fin << "10\n";
		else if (ele->type == 8) fin << "12\n";
	}
	// file output complete

	fin.close();

}


//______________________________________________________________________

extern "C"
{

// refine convex edges in mesh when generating TetGen tetrahedrel dense mesh
int refineEdgesTetgen (MESH *msh, double MIN_ANGLE, double MAX_ANGLE, double edgeLength, double* mtroutlist) {
	
	int j = -1, k = 0, z = 0, e;
	bool found;
	FACE* _face;

	std::vector<int> vecEdge;
	std::vector<int>::iterator fit;

	std::multimap<int,int> vtx_edge; // map to store <vertex index, insertion index> instead of using hash
	std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> ret;

	std::vector<Edge> _edges, _convexEdges;

	MESH* _msh = MESH_Copy (msh);
	MESH_Update (_msh, NULL, NULL, NULL); // reference configuration

	// loop through surface faces
	for (_face = _msh->faces; _face; _face = _face->n)
	{
		// loop through edges of face & store in edge list
		k = _face->type;
		for (int i = 0; i < k; ++i)
		{
			// introduce z variable to compensate for edge at (2,0)
			if (i != k-1) z = i+1;
			else z = 0;

			// bool to indicate if edge is already in vector
			e = 0;
			found = false;

			// check if _face->nodes[i] exists in vtx_edge
			if (vtx_edge.count(_face->nodes[i]) > 0)
			{
				// loop through values
				ret = vtx_edge.equal_range(_face->nodes[i]);
				for (std::multimap<int,int>::iterator it = ret.first; it != ret.second; ++it, ++e)
				{	
					if ((int)vecEdge.size() > e) 
						vecEdge[e] = it->second;
					else 
						vecEdge.push_back(it->second);
				}
			}

			// check if vecEdge is not empty
			if (e > 0)
			{
				// check if _face->nodes[z] exists in vtx_edge
				if (vtx_edge.count(_face->nodes[z]) > 0)
				{
					// loop through values
					ret = vtx_edge.equal_range(_face->nodes[z]);
					for (std::multimap<int,int>::iterator it = ret.first; it != ret.second; ++it)
					{
						fit = find (vecEdge.begin(), vecEdge.begin()+e, it->second);
						if ( fit != (vecEdge.begin()+e) )
						{
							found = true;
							break;
						}
					}
				}
			}

			// if edge is not found
			if (!found)
			{
				if (z != k-1) _edges.push_back( Edge(_face->nodes[i], _face->nodes[z], _face->nodes[z+1], _face, NULL) );
				else _edges.push_back( Edge(_face->nodes[i], _face->nodes[z], _face->nodes[0], _face, NULL) );

				++j;
				vtx_edge.insert( std::pair<int,int>(_face->nodes[i], j) );
				vtx_edge.insert( std::pair<int,int>(_face->nodes[z], j) );
			} else {
				_edges[*fit].f2_VertIndx = findIndx(_face, _edges[*fit]);
				_edges[*fit].fac2 = _face;
			}

		}
	}

	// find the convex edges
	// std::cout << "-->Number of edges in model: " << _edges.size() << std::endl;
	findConvexEdges(_msh, _edges, MIN_ANGLE, MAX_ANGLE, _convexEdges);
	// std::cout << "Number of convex edges in model: " << _convexEdges.size() << std::endl;

	// output mtrpointlist to be read in TetGen
	std::vector<double> _mtrout;
	outputMTR(_msh, _convexEdges, edgeLength, _mtrout);
	//std::cout << "mtrout: " << _mtrout.size() << std::endl;

	memcpy (mtroutlist, &_mtrout[0], sizeof(double)*_mtrout.size());
#if 0
	// output files check
	std::string edgeMeshFname_convex = "convex_edges.vtk";
	edges2VTK(edgeMeshFname_convex, _msh, _convexEdges);

	std::string edgeMeshFname_all = "all_edges.vtk";
	edges2VTK(edgeMeshFname_all, _msh, _edges);

	std::string VTKFname = "mesh_original.vtk";
	EXPORT_SOLFEC_MESH(VTKFname, _msh);
#endif

	_edges.clear();
	_convexEdges.clear();
	MESH_Destroy (_msh);

	return _mtrout.size();

}

}
