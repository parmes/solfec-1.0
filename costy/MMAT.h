//
//  MMAT.h
//  Maths Matrix Class v1
//
//  Created by Costy Kodsi on 03/04/2013.
//  Copyright (c) 2013
//  
//  Remaining Work:
//  (i) Make (class) into template format
//  (ii) Make code more efficient:
//      a. inv function (pass result matrix as argument saving on memory copying)
//      b. tp function (pass result matrix as argument saving on memory copying)
//  (iii) Add cross product & dot product
//  (iv) Set inv values to 0 where required
//  (v) determinant could become too large for double, store as power 10
//  (vi) add const to member functions
//  (vii) add 1d and 2d input (i.e. read) array functionality
//
//  Changes To Be Made:
//  (i) Change class name to MAT from MMAT
//
//
//  Usage Guide:
//
//  USE LIKE MATLAB
//  JUST REMEMBER: INITIALISE TYPE, only need to specify dimensions if you will specify each value
//  e.g. MMAT A, B, C(2,2);
//  e.g. C(1,1) = 5; C(1,2) = 10; ... 
//  e.g. A = C + C;
//  
//  -C not possible, instead use -1*A
//  Access component value, C.val(1,1) ...
//
//  Access matrix dimension: C.rowno(), C.colno()
//  Available functions:
//  A = inv(C);
//  A = tp(C);
//  A = det(C);
//  (note: be careful with inv and det, no warning for poorly conditioned matrices)
//
//  View Result: C.output();
//


#ifndef maths_matrix_MMAT_h
#define maths_matrix_MMAT_h
/*
#include <stdlib.h>
#include <cstdio>
#include <math.h>
#include <cmath>
#include <Accelerate/Accelerate.h>
*/
using namespace std;

class MMAT;
MMAT tp(const MMAT& );
MMAT inv(const MMAT& );
double det(const MMAT& );


class MMAT {
    
public:
    
    MMAT(); //default constructor
    MMAT(const MMAT& p); //copy constructor
    MMAT(const int RowSize, const int ColSize); //assign size 
    ~MMAT(); //destructor
    
    MMAT& operator= (const MMAT& p); //equal operator, e.g. A=B
    
    //addition operators
    MMAT operator+ (const double a);
    friend MMAT operator+ (const double a, const MMAT& p);
    MMAT operator+ (const MMAT& p);
    MMAT& operator+= (const MMAT& p);
    friend MMAT operator+ (const MMAT& p, const MMAT& q);
    
    //subtraction operators
    MMAT operator- (const double a); 
    friend MMAT operator- (const double a, const MMAT& p);
    MMAT operator- (const MMAT& p);
    MMAT& operator-= (const MMAT& p);
    friend MMAT operator- (const MMAT& p, const MMAT& q);
   
    //multiplication operators
    MMAT operator* (const double a);
    friend MMAT operator* (const double a, const MMAT& p);
    MMAT operator* (const MMAT& p);
    friend MMAT operator* (const MMAT& p, const MMAT& q);
    
    //division operators
    MMAT operator/ (const double a);
    friend MMAT operator/ (const double a, const MMAT& p);
    MMAT operator/ (const MMAT& p);
    friend MMAT operator/ (const MMAT& p, const MMAT& q);
    
    MMAT& tp(); //transpose matrix
    MMAT& zero(); //set all elements of matrix to zero
    MMAT& identity(); //set all A[i][i] elements to 1
    
    //access matrix component values
    double& operator() (const int row, const int col)
    {
        if (row>0 && row<=nrows && col>0 && col<=ncols)
        {
            return m[row-1][col-1];  
        }
        else
        {
            cout << "Error: Object Not Defined Or Data Request Outside Bounds. \n";
            exit (EXIT_FAILURE);
        }
    }
    
    //access matrix component values
    double val(const int row, const int col) const
    {
        if (row>0 && row<=nrows && col>0 && col<=ncols)
        {
            return m[row-1][col-1];  
        }
        else
        {
            cout << "Error: Object Not Defined Or Data Request Outside Bounds. \n";
            exit (EXIT_FAILURE);
        }
    }
    
    //access number of rows in matrix
    int rowno() const
    {return nrows;}
    
    //access number of columns in matrix
    int colno() const
    {return ncols;}
    
    //print output in matrix form
    void output() const;
    
private:
    
    //basic maths functions (in order): 
    //A[i][j] + b, A[i][j] - b, b - A[i][j], A[i][j] * b, A[i][j] / b, b / A[i][j]
    MMAT& add(const double value);
    MMAT& subt1(const double value);
    MMAT& subt2(const double value);
    MMAT& mult(const double value);
    MMAT& div1(const double value);
    MMAT& div2(const double value);
    
    //allocate matrix
    void initMat()
    {
        m = new double*[nrows];
        for (int i = 0; i < nrows; ++i)
        {
            m[i] = new double[ncols];
        }
    }
    
    //declaring private variables
    int nrows, ncols;
    double **m;
    
};
    


//-----------------------------------------------------------------------//
    
//constructor
MMAT::MMAT()
{
    m = NULL;
    nrows = 0;
    ncols = 0;
}

//dynamically allocate matrix give no of rows and cols
MMAT::MMAT(const int RowSize, const int ColSize) : nrows(RowSize),ncols(ColSize)
{
    initMat();
    
    for (int i = 0; i < nrows; ++i)
    {
        for (int j = 0; j < ncols; ++j)
        {
            m[i][j] = 0;
        }
    }
}


//copy constructor
MMAT::MMAT(const MMAT& p) : nrows(p.nrows),ncols(p.ncols)
{
    initMat();
    
    for (int i = 0; i < nrows; ++i)
    {
        for(int j = 0; j < ncols; ++j)
        {
            m[i][j] = p.m[i][j];
        }
    }
}


//destructor 
MMAT::~MMAT()
{
    for (int i = 0; i < nrows; ++i)
    {
        delete [] m[i];
    }
    delete [] m;
    m = NULL;
}


//-----------------------------------------------------------------------//

//matrix equal operator
MMAT& MMAT::operator= (const MMAT& p)
{
    if (this == &p)
    {
        return *this;
    }
    else
    {
        if (nrows != p.nrows || ncols != p.ncols)
        {
            this->~MMAT();
            nrows = p.nrows;
            ncols = p.ncols;
            initMat();
        }
        
        for (int i = 0; i < nrows; ++i)
        {
            for (int j = 0; j < ncols; ++j)
            {
                m[i][j] = p.m[i][j];
            }
        }
        
        return *this;
    }
    
}


//-----------------------------------------------------------------------//

//add matrix to double (RHS), e.g. A + 1.467
MMAT MMAT::operator+ (double a)
{
    MMAT matrix(*this);
    matrix.add(a);
    return matrix;
    
}


//add matrix to double, e.g. 1.46789 + A
MMAT operator+ (const double a, const MMAT& p)
{
    MMAT matrix(p);
    matrix.add(a);
    return matrix;
}


//add matrix to matrix, e.g. A + B, (RHS)
MMAT MMAT::operator+ (const MMAT& p)
{
    // check dimensions
    if (nrows == p.nrows && ncols == p.ncols)
    {
        MMAT matrix(*this);
        
        for (int i = 0; i < nrows; ++i)
        {
            for (int j =0; j < ncols; ++j)
            {
                matrix.m[i][j] += p.m[i][j];
            }
        }
        
        return matrix;
    }
    else
    {
        cout << "Error Matrix Dimensions for Addition Operation Do Not Agree. \n";
        exit (EXIT_FAILURE);
    }
}


//addition operation += (i.e. A = A + B)
MMAT& MMAT::operator+= (const MMAT& p)
{
    MMAT matrix(*this);
    
    if (matrix.nrows == p.nrows && matrix.ncols == p.ncols)
    {
        *this = *this + p;
        return *this;   
    }
    else
    {
        cout << "Error Matrix Dimensions for Addition Operation Do Not Agree. \n";
        exit (EXIT_FAILURE);
    }
}


//add matrix to matrix, e.g. B + A
MMAT operator+ (const MMAT& p, const MMAT& q)
{
    // check dimensions
    if (p.nrows == q.nrows && p.ncols == q.ncols)
    {
        MMAT matrix(p.nrows,p.ncols);
        for (int i = 0; i < p.nrows; ++i)
        {
            for (int j =0; j < p.ncols; ++j)
            {
                matrix.m[i][j] = p.m[i][j] + q.m[i][j];
            }
        }
        
        return matrix;
    }
    else
    {
        cout << "Error Matrix Dimensions for Addition Operation Do Not Agree. \n";
        exit (EXIT_FAILURE);
    }
}


//-----------------------------------------------------------------------//

//subtract double from matrix, e.g. A - 1.426
MMAT MMAT::operator- (const double a)
{
    MMAT matrix(*this);
    matrix.subt1(a);
    return matrix;
}


// subtract matrix from double, e.g. 1.46789 - A
MMAT operator- (const double a, const MMAT& p)
{
    MMAT matrix(p);
    matrix.subt2(a);
    return matrix;
}


//subtract matrix from matrix (RHS)
MMAT MMAT::operator- (const MMAT& p)
{
    // check dimensions
    if (nrows == p.nrows && ncols == p.ncols)
    {
        MMAT matrix(*this);
        
        for (int i = 0; i < nrows; ++i)
        {
            for (int j =0; j < ncols; ++j)
            {
                matrix.m[i][j] -= p.m[i][j];
            }
        }
        
        return matrix;
    }
    else
    {
        cout << "Error Matrix Dimensions for Subtraction Operation Do Not Agree. \n";
        exit (EXIT_FAILURE);
    }
}


//subtraction operation (i.e. A = A - B)
MMAT& MMAT::operator-= (const MMAT& p)
{
    MMAT matrix(*this);
    
    if (matrix.nrows == p.nrows && matrix.ncols == p.ncols)
    {
        *this = *this - p;
        return *this;   
    }
    else
    {
        cout << "Error Matrix Dimensions for Subtraction Operation Do Not Agree. \n";
        exit (EXIT_FAILURE);
    }
}


// subtract matrix from matrix, e.g. B - A
MMAT operator- (const MMAT& p, const MMAT& q)
{
    // check dimensions
    if (p.nrows == q.nrows && p.ncols == q.ncols)
    {
        MMAT matrix(p.nrows,p.ncols);
        for (int i = 0; i < p.nrows; ++i)
        {
            for (int j =0; j < p.ncols; ++j)
            {
                matrix.m[i][j] = p.m[i][j] - q.m[i][j];
            }
        }
        
        return matrix;
    }
    else
    {
        cout << "Error Matrix Dimensions for Subtraction Operation Do Not Agree. \n";
        exit (EXIT_FAILURE);
    }
}


//-----------------------------------------------------------------------//

//multiply matix by double, e.g. A * a
MMAT MMAT::operator* (const double a)
{
    MMAT matrix(*this);
    matrix.mult(a);
    return matrix;
}

//multiply matrix with double, e.g. 1.46789 * A
MMAT operator* (const double a, const MMAT& p)
{
    MMAT matrix(p);
    matrix.mult(a);
    return matrix;
}


//multiply matrix with matrix, e.g. A * B
MMAT MMAT::operator* (const MMAT& p)
{
    // check dimensions
    if (ncols == p.nrows)
    {
        MMAT matrix(nrows,p.ncols);
        
        for (int i = 0; i < nrows; ++i)
        {
            for (int j =0; j < p.ncols; ++j)
            {
                matrix.m[i][j] = 0.0;
                
                for (int k = 0; k < ncols; ++k)
                {
                    matrix.m[i][j] += m[i][k] * p.m[k][j];                      
                }
            }
        }
        
        return matrix;
    }
    else
    {
        cout << "Error Matrix Dimensions for Multiplication Operation Do Not Agree. \n";
        exit (EXIT_FAILURE);
    }
}


//multiply matrix with matrix, e.g. B * A
MMAT operator* (const MMAT& p, const MMAT& q)
{
    // check dimensions
    if (p.ncols == q.nrows)
    {
        MMAT matrix(p.nrows,q.ncols);
        
        for (int i = 0; i < p.nrows; ++i)
        {
            for (int j =0; j < q.ncols; ++j)
            {
                matrix.m[i][j] = 0.0;
                
                for (int k = 0; k < p.ncols; ++k)
                {
                    matrix.m[i][j] += p.m[i][k] * q.m[k][j];                      
                }
            }
        }
        
        return matrix;
    }
    else
    {
        cout << "Error Matrix Dimensions for Multiplication Operation Do Not Agree. \n";
        exit (EXIT_FAILURE);
    }
}


//-----------------------------------------------------------------------//

//divide matrix by double, e.g. A / 1.426
MMAT MMAT::operator/ (const double a)
{
    MMAT matrix(*this);
    matrix.div1(a);
    return matrix;
}

//divide matrix by double, e.g. 1.46789 / A
MMAT operator/ (const double a, const MMAT& p)
{
    MMAT matrix(p);
    matrix.div2(a);
    return matrix;
}


//divide matrix by matrix,e.g. A / B (RHS)
MMAT MMAT::operator/ (const MMAT& p)
{
    // check dimensions
    if (nrows == ncols && p.nrows == p.ncols && nrows == p.nrows && ncols == p.ncols)
    {
        MMAT temp(*this);
        MMAT matrix(nrows,ncols);
        matrix = temp * inv(p);
        return matrix;
    }
    else
    {
        cout << "Error Matrix Dimensions for Division Operation Do Not Agree. \n";
        exit (EXIT_FAILURE);
    }
}


//divide matrix by matrix, e.g. B / A
MMAT operator/ (const MMAT& p, const MMAT& q)
{
    // check dimensions
    if (p.nrows == p.ncols && q.nrows == q.ncols && p.nrows == q.nrows && p.ncols == q.ncols)
    {
        MMAT matrix(p.nrows,p.ncols);
        matrix = p * inv(q);
        return matrix;
    }
    else
    {
        cout << "Error Matrix Dimensions for Addition Operation Do Not Agree. \n";
        exit (EXIT_FAILURE);
    }
}


//-----------------------------------------------------------------------//

//transpose matrix
MMAT& MMAT::tp()
{
    // check if square matrix
    if (nrows == ncols)
    {
        double temp;
        
        for (int i = 0; i < nrows; ++i)
        {
            for (int j = i + 1; j < ncols; ++j)
            {
                temp = m[i][j];
                m[i][j] = m[j][i];
                m[j][i] = temp;
            }
        }
    
    }
    else
    {
        MMAT temp(ncols,nrows);
        
        for (int i = 0; i < nrows; ++i)
        {
            for (int j = 0; j < ncols; ++j)
            {
                temp.m[j][i] = m[i][j];
            }
        }
        
        *this = temp;
    }
    
    return *this;
}


//set all matrix elements to zero
MMAT& MMAT::zero()
{
    for (int i = 0; i < nrows; ++i)
    {
        for (int j = 0; j < ncols; ++j)
        {
            m[i][j] = 0;

        }
    }
    
    return *this;
}


//set A[i][i] matrix elements to 1
MMAT& MMAT::identity()
{
    if (nrows == ncols)
    {
        for (int i = 0; i < nrows; ++i)
        {
            for (int j = 0; j < ncols; ++j)
            {
                if (i == j)
                {
                    m[i][j] = 1;
                }
            }
        }
        
    }
    else
    {
        cout << "Error: Not Square Matrix! Therefore Cannot Exectue Request.\n";
    }
    
    return *this;
}

//-----------------------------------------------------------------------//

//output matrix
void MMAT::output() const
{
    int i, j;
    
    if (nrows>0 && ncols>0)
    {
        cout << "[";
        for (i = 0; i < nrows; ++i)
        {
            if (i > 0)
            {
                cout << " ";
            }
            for (j = 0; j < ncols-1; ++j)
            {
                cout << m[i][j] << ", ";
            }
            if (i < nrows-1)
            {
                cout << m[i][ncols-1] << ";\n";
            }
            else
            {
                cout << m[i][ncols-1] << "]\n";
            }   
        }
    }
    
}    


//-----------------------------------------------------------------------//

//addition function
MMAT& MMAT::add(const double value)
{
    for (int i = 0; i < nrows; ++i)
    {
        for (int j = 0; j < ncols; ++j)
        {
            m[i][j] += value; 
        }
    }
    
    return *this;
}


//subtraction function 1 (m[i][j] - val)
MMAT& MMAT::subt1(const double value)
{
    for (int i = 0; i < nrows; ++i)
    {
        for (int j = 0; j < ncols; ++j)
        {
            m[i][j] -= value; 
        }
    }
    
    return *this;
}


//subtraction function 2 (val - m[i][j])
MMAT& MMAT::subt2(const double value)
{
    for (int i = 0; i < nrows; ++i)
    {
        for (int j = 0; j < ncols; ++j)
        {
            m[i][j] = value - m[i][j]; 
        }
    }
    
    return *this;
}


// multiplication function
MMAT& MMAT::mult(const double value)
{
    for (int i = 0; i < nrows; ++i)
    {
        for (int j = 0; j < ncols; ++j)
        {
            m[i][j] *= value; 
        }
    }
    
    return *this;
}


//division function 1 (A[i][j] / val)
MMAT& MMAT::div1(const double value)
{
    for (int i = 0; i < nrows; ++i)
    {
        for (int j = 0; j < ncols; ++j)
        {
            m[i][j] /= value; 
        }
    }
    
    return *this;
}


//division function 2 (val / A[i][j])
MMAT& MMAT::div2(const double value)
{
    for (int i = 0; i < nrows; ++i)
    {
        for (int j = 0; j < ncols; ++j)
        {
            m[i][j] = value / m[i][j]; 
        }
    }
    
    return *this;
}

//-----------------------------------------------------------------------//

//transpose function, A = tp(B)
MMAT tp(const MMAT& p)
{
    MMAT temp(p.colno(),p.rowno());
    
    for (int i = 1; i < (1 + p.rowno()); ++i)
    {
        for (int j = 1; j < (1 + p.colno()); ++j)
        {
            temp(j,i) = p.val(i,j);
        }
    }
    
    return temp; 
        
}


//matrix inverse
MMAT inv(const MMAT& p)
{
    if (p.rowno() == p.colno())
    {
        int m, n, i, j; 
        m = p.rowno();
        n = p.colno();
        
        MMAT temp(m,n);
        
        //define single array for input into LAPACK function
        double* a;
        a = new double[m*n];
        for (i = 0; i < m; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                a[i*m+j] = p.val(i+1,j+1);
            }
        }
        
        //IPIV = min(rowno,colno), redundant for square matrix
        int* ipiv;
        if (m <= n)
        {
            ipiv = new int[m];
        }else
        {
            ipiv = new int[n];
        }
        
        int lwork = m*m;
        double *work = new double[lwork];
        int info;
        
        //LU decomoposition of a general matrix
        dgetrf_(&m,&n,a,&m,ipiv,&info);
        
        //Generate inverse of a matrix given its LU decomposition
        dgetri_(&m,a,&m,ipiv,work,&lwork,&info);
        
        //define output function
        for (i = 0; i < m; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                temp(i+1,j+1) = a[i*m+j];
            }
        }
        
        delete ipiv;
        delete work;
        
        return temp;
    }
    else
    {
        cout << "Error: Not Square Matrix! Therefore Cannot Inverse.\n";
        exit (EXIT_FAILURE);
    }
    
}


//determinant of a matrix, function gives correct value but not sign sometimes!!!
double det(const MMAT& p)
{
    if (p.rowno() == p.colno())
    {
        int m, n, i, j, negative = 0;
        double determinant;
        m = p.rowno();
        n = p.colno();
        
        MMAT temp(m,n);
        
        //define single array for input into LAPACK function
        double* a;
        a = new double[m*n];
        for (i = 0; i < m; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                a[i*m+j] = p.val(i+1,j+1);
            }
        }
        
        //IPIV = min(rowno,colno), redundant for square matrix
        int* ipiv;
        if (m <= n)
        {
            ipiv = new int[m];
        }else
        {
            ipiv = new int[n];
        }
        
        int info;
        
        dgetrf_(&m,&n,a,&m,ipiv,&info);
        
        //define output function
        for (i = 0; i < m; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                temp(i+1,j+1) = a[i*m+j];
            }
        }
        
        //calc determinant
        determinant = 1;
        for (i = 0; i < m; ++i)
        {
            determinant *= temp(i+1,i+1);
        }
        
        //calc sign
        for (i = 0; i < m; ++i)
        {
            if (ipiv[i] != (i+1))
            {
                negative = !negative; 
            }
        }
        
        
        delete ipiv;
        
        return negative?-determinant:determinant;
    }
    else
    {
        cout << "Error: Not Square Matrix! Therefore Cannot Create Determinant.\n";
        exit (EXIT_FAILURE);
    }

}

#endif