//
//  opt.h
//  otimisation class
//


#include <SO3.h>  
#include <matrix.h>
#include <cmath>
#include <stdio.h>      
#include <stdlib.h>

#define PI 3.14159265358

 
class OPT; // class declaration


class OPT {
    
public:

    // Newton optimisation with local parameterisation of the manifold SO(3)
    matrix<double> newtSO3(matrix<double>& w, const int nTrial, SO3& func, const double tol)
    {
        double err = 1;
        matrix<double> h(3,1), R(3,3), R0(3,3), jacob, hess, err_tmp(1,1);
		
		R.identity();
        R0.identity();

        h.zero();
        err_tmp.zero();
        
        // carry out optimisation
        for (int k = 0; k < nTrial; ++k)
        {
            if (err < tol) return h;
            else {
                //compute gradient & hessian
                jacob = func.jacobian(R0, R);
                hess = func.hessian(R0, R);

                // std::cout << "Jacobian:\n";
                // jacob.output();
                // std::cout << "Hessian:\n";
                // hess.output();
                
                w = -1 * inv(hess) * jacob;
                func.updateR0(R0,w);

                // std::cout << "w:\n";
                // w.output();
                
                err_tmp = (tp(jacob) * -1 * w) / 2;
                err = err_tmp.val(1,1);
                h = w + h;
                
                std::cout << "Iteration: " << k+1 << "; w: " << h.val(1,1) << ", " <<  h.val(2,1) << ", " << h.val(3,1) << "; Error: " << err << std::endl; 
            }   
        }

        // if optimisation does not converge
        if (err >= tol)
        {
            std::cout << "Error: Convergence not achieved.\n";;
            std::cout << "Try increasing trial no and/or increasing tolerance.\n";
            std::cout << "Please always remember to check input!\n";
            exit(EXIT_FAILURE);
        }
            
        return h;
    }

};