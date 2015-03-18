#ifndef __costy__
#define __costy__

// refine convex edges in mesh (looking from outside) when generating TetGen tetrahedrel refined mesh
// returns the size of mtroutlist (allocated inside)
#if __cplusplus
extern "C"
{
#endif

int refineEdgesTetgen (MESH *msh, double MIN_ANGLE, double MAX_ANGLE, double edgeLength, double* mtroutlist);

int rbmm_main(BODY* _bod, double volume, double quality, char* _output);

#if __cplusplus
}
#endif

#if __cplusplus
//_____________________________________________________________________________________________________________________________________

// Vector class
class Vector_3 {

public:

	Vector_3() : x(0.0), y(0.0), z(0.0) {};
	Vector_3(double _x, double _y, double _z) : x(_x), y (_y), z(_z) {};

	void normalize() {
		double norm = sqrt(x*x + y*y + z*z);
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

static double inner(const Vector_3& a, const Vector_3& b) {

	return (a.x*b.x + a.y*b.y + a.z*b.z);

}

static double norm(const Vector_3& a) {

	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);

}

static void cross(const Vector_3& a, const Vector_3& b, Vector_3& c) {

	c.x = a.y*b.z - a.z*b.y;
	c.y = a.z*b.x - a.x*b.z;
	c.z = a.x*b.y - a.y*b.x;

}

// calculate angle between two vectors
static double angle(const Vector_3& a, const Vector_3& b) {

  Vector_3 c;
  cross(a,b,c);

  return atan2( norm(c), inner(a,b) );

}
//_____________________________________________________________________________________________________________________________________
#endif

#endif
