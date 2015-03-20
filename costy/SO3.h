//
//  S03.h
//  maths matrix
//

#include <stdio.h>      
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <matrix.h> 

#define PI 3.14159265358

// class & function declarations
class SO3;
double norm(const matrix<double>& x);
void spin(matrix<double>& A, const matrix<double>& B, const int i);


// calc norm of vector
double norm(const matrix<double>& x)
{
    double temp = 0.0;
    for (int i = 1; i < ( 1 + x.rowno() ); ++i) temp += ( x.val(i,1) * x.val(i,1) );
    
    return sqrt(temp);
}


// spin function
void spin(matrix<double>& A, const matrix<double>& B, const int i)
{
    A(1,2) = -1 * B.val(i,3);
    A(1,3) = B.val(i,2);
    A(2,1) = B.val(i,3);
    A(2,3) = -1 * B.val(i,1);
    A(3,1) = -1 * B.val(i,2);
    A(3,2) = B.val(i,1);
}



class SO3 {
    
public:
    
    SO3(std::vector<Facet>& fac)
    {
        
        initialise();

        nofaces = fac.size();
        normals(nofaces, 3);
        integ(nofaces, 3);

        for (int i = 0; i < fac.size(); ++i)
        {
            normals(i, 1) = fac[i].normal.x;
            normals(i, 2) = fac[i].normal.y;
            normals(i, 3) = fac[i].normal.z;

            integ(i, 1) = fac[i].integ.x;
            integ(i, 2) = fac[i].integ.y;
            integ(i, 3) = fac[i].integ.z;
        }

    }
    
    
    //destructor
    ~SO3()
    {
    }

        
    //initialise matrix R
    matrix<double> initR(const matrix<double>& w)
    {
        matrix<double> Jn, R;
        double theta = norm(w);
        Jn =  (  (w.val(1,1) * Jwx) + (w.val(2,1) * Jwy) + (w.val(3,1) * Jwz)  ) / theta;
        R = I3 + (sin(theta) * Jn) + (  2 * sin(theta/2) * sin(theta/2) * (Jn * Jn)  );
        return R;
    }


    //update matrix R using new values of w
    void updateR0(matrix<double>& R0, const matrix<double>& w)
    {
        matrix<double> Jn;
        double theta = norm(w);
        Jn =  (  (w.val(1,1) * Jwx) + (w.val(2,1) * Jwy) + (w.val(3,1) * Jwz)  ) / theta;
        R0 = R0 * (   I3 + (sin(theta) * Jn) + (  2 * sin(theta/2) * sin(theta/2) * (Jn * Jn)  )   );
    }

    
    //calc objective function val
    matrix<double> objective(const matrix<double>& R0, const matrix<double>& R)
    {
        matrix<double> A(3,1), tmp_normal(3,3), tmp_integ(3,1), S, Q, obj;

		S = tp(R);
        Q = tp(R0);

		A.zero();
		tmp_normal.zero();
        tmp_integ.zero();
        
        for (int i = 1; i < (1+nofaces); ++i)
        {     
            spin(tmp_normal, normals, i);
			
            tmp_integ(1,1) = integ.val(i,1);
            tmp_integ(2,1) = integ.val(i,2);
            tmp_integ(3,1) = integ.val(i,3);

			A += (tmp_normal * S * Q * tmp_integ);
			
			//tmp_normal.zero();
            //tmp_integ.zero();
        }
        
        obj = 0.5 * tp(A) * A;
        
        return obj;
    }

   
    //calc gradient
    matrix<double> jacobian(const matrix<double>& R0, const matrix<double>& R)
    {
        matrix<double> tmp_normal(3,3), tmp_integ(3,1);
        matrix<double> S, Q, jacob(3,1), A(3,1), B(3,1), C(3,1), D(3,1), tmp(1,1);
        
        S = tp(R);
        Q = tp(R0);
        
        A.zero();
        B.zero();
        C.zero();
        D.zero();
		tmp_normal.zero();
        tmp_integ.zero();
        tmp.zero();
        
        for (int i = 1; i < (1+nofaces); ++i)
        {     
            spin(tmp_normal, normals, i);
			
            tmp_integ(1,1) = integ.val(i,1);
            tmp_integ(2,1) = integ.val(i,2);
            tmp_integ(3,1) = integ.val(i,3);
            
			A += (tmp_normal * S * Q * tmp_integ);
            B += (tmp_normal * -1 * Jwx * Q * tmp_integ);
            C += (tmp_normal * -1 * Jwy * Q * tmp_integ);
            D += (tmp_normal * -1 * Jwz * Q * tmp_integ);
			
        }

        tmp = tp(A) * B;
        jacob(1,1) = tmp.val(1,1);
        tmp = tp(A) * C;
        jacob(2,1) = tmp.val(1,1);
        tmp = tp(A) * D;
        jacob(3,1) = tmp.val(1,1);
        
        /*
        
        for (int i = 0; i < 3; ++i)
        {
            temp.zero();
			A.zero();
			B.zero();
			
            for (int j = 1; j < (1+nofaces); ++j)
            {           
				tmp_normal.zero();
                tmp_integ.zero();
      
                spin(tmp_normal, normals, j);
				
                tmp_integ(1,1) = integ.val(j,1);
                tmp_integ(2,1) = integ.val(j,2);
                tmp_integ(3,1) = integ.val(j,3);
                
				A += (tmp_normal * tp(R) * tmp_integ);
				B += (tmp_normal * (-1 * (*pJ[i])) * tp(Ro)  * tmp_integ);
            }

			temp = tp(A) * B;
            jacob(i+1,1) = temp.val(1,1);
        }
        */
        return jacob;
    }
    
    
    //calc hessian
    matrix<double> hessian(const matrix<double>& R0, const matrix<double>& R)
    {
        matrix<double> tmp_normal(3,3), tp_tmp_normal(3,3), tmp_integ(3,1), tp_tmp_integ(1,3);
        matrix<double> S, Q, hess(3,3), A(1,3), B(3,1), C(3,1), D(3,1), tmp(1,1);
        
        S = tp(R);
        Q = tp(R0);
        
        A.zero();
        B.zero();
        C.zero();
        D.zero();
		tmp_normal.zero();
        tp_tmp_normal.zero();
        tmp_integ.zero();
        tp_tmp_integ.zero();
        tmp.zero();
        
        for (int i = 0; i < 3; ++i)
        {
            for (int k = 0; k < 3; ++k)
            {
                
                for (int j = 1; j < (1+nofaces); ++j)
                {              
                    spin(tmp_normal, normals, j);
                    tp_tmp_normal = tp(tmp_normal);
                    
                    tmp_integ(1,1) = integ.val(j,1);
                    tmp_integ(2,1) = integ.val(j,2);
                    tmp_integ(3,1) = integ.val(j,3);
                    
                    tp_tmp_integ(1,1) = integ.val(j,1);
                    tp_tmp_integ(1,2) = integ.val(j,2);
                    tp_tmp_integ(1,3) = integ.val(j,3);
                    
					A += (tp_tmp_integ * R0 * (*pJ[k]) * tp_tmp_normal);   
					B += (tmp_normal * -1 * (*pJ[i]) * Q * tmp_integ);
                    
					C += (tmp_normal * S * Q * tmp_integ);
                    D += (tmp_normal * -1 * (*pH[i][k]) * Q * tmp_integ);
                    
                    //tmp_normal.zero();
                    //tmp_integ.zero();
                }
                
                tmp = (A * B) + (tp(C) * D);
                hess(i+1,k+1) = tmp.val(1,1);

				A.zero();
				B.zero();
				C.zero();
                D.zero();
                tmp.zero();
            }
        }
        
        return hess;
         
    }
    
    
private:
    
    void initialise()
    {
        matrix<double> id3(3,3), Jx(3,3), Jy(3,3), Jz(3,3);
        
        id3.identity();
        
        Jx(1,1) = 0;
        Jx(1,2) = 0;
        Jx(1,3) = 0;
        Jx(2,1) = 0;
        Jx(2,2) = 0;
        Jx(2,3) = -1;
        Jx(3,1) = 0;
        Jx(3,2) = 1;
        Jx(3,3) = 0;
        
        Jy(1,1) = 0;
        Jy(1,2) = 0;
        Jy(1,3) = 1;
        Jy(2,1) = 0;
        Jy(2,2) = 0;
        Jy(2,3) = 0;
        Jy(3,1) = -1;
        Jy(3,2) = 0;
        Jy(3,3) = 0;
        
        Jz(1,1) = 0;
        Jz(1,2) = -1;
        Jz(1,3) = 0;
        Jz(2,1) = 1;
        Jz(2,2) = 0;
        Jz(2,3) = 0;
        Jz(3,1) = 0;
        Jz(3,2) = 0;
        Jz(3,3) = 0;
        
        I3 = id3;
        Jwx = Jx;
        Jwy = Jy;
        Jwz = Jz;
        
        pJ[0] = &Jwx;
        pJ[1] = &Jwy;
        pJ[2] = &Jwz;
        
        dRxx = Jwx * Jwx;
        dRyy = Jwy * Jwy;
        dRzz = Jwz * Jwz;
        dRxy = 0.5 * ( (Jwx*Jwy) + (Jwy*Jwx) );
        dRxz = 0.5 * ( (Jwx*Jwz) + (Jwz*Jwx) );
        dRyz = 0.5 * ( (Jwy*Jwz) + (Jwz*Jwy) );
        
        pH[0][0] = &dRxx;
        pH[0][1] = &dRxy;
        pH[0][2] = &dRxz;
        
        pH[1][0] = &dRxy;
        pH[1][1] = &dRyy;
        pH[1][2] = &dRyz;
        
        pH[2][0] = &dRxz;
        pH[2][1] = &dRyz;
        pH[2][2] = &dRzz;
        
    }
    
    int nofaces;
    matrix<double> normals, integ;
    matrix<double> I3, Jwx, Jwy, Jwz;
    matrix<double> dRxx, dRyy, dRzz, dRxy, dRxz, dRyz, *pJ[3], *pH[3][3];

};