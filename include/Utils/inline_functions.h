/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/

#ifndef SMARTMATH_INLINEFUNCTIONS_H
#define SMARTMATH_INLINEFUNCTIONS_H

#include <vector>
#include <cmath>
#include "../LinearAlgebra/Eigen/Eigen"

namespace smartmath{

const double ZERO = 1.0e-15;

template <class T>
T inverse(T x){
    if(fabs(x)<=ZERO){
        std::cout<<"ERROR: Division by zero."<<std::endl;
        throw std::exception();
    }
    return 1.0/x;
}

//MATH STUFFS
inline int factorial(int n)
{
    if(n<0)
        smart_throw("FACT: factorial of non positive integer does not exist");

    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

inline int combination(int n, int k)
{
    int max = std::max(n,k);
    int min = std::min(n,k);
    int res = 1;
    int j = 1;

    while(j<=min){
        res *= max+j;
        j++;
    }

    return res/factorial(min);
}


inline void rep(std::vector<std::vector<int> > &res, const std::vector<int> &values, std::vector<int> &item, unsigned int count){
    if (count < item.size()){
        for (unsigned int i = 0; i < values.size(); i++) {
            item[count] = values[i];
            unsigned int tmp_count = count + 1;
            rep(res, values, item, tmp_count);
        }
    }else{
        res.push_back(item);
    }
}


inline void variations(const std::vector<int> values, const int k, std::vector<std::vector<int> > &res){
    res.clear();

    std::vector<int> item(k);
    rep(res, values, item, 0);
}


typedef double (*fun)(double);
inline double bisection_method(fun f, double lb, double ub, double prec){
    double f_low = f(lb);
    double f_up  = f(ub);
    if( (f_low*f_up) < 0 && (ub-lb) <= prec )
    {
        return lb;
    }

    double temp   = (ub+lb)/2;
    double f_temp = f(temp);

    if( (f_temp*f_low) > 0 )
        lb = temp;
    
    if( (f_temp*f_up) > 0 )
        ub = temp;
     
return bisection_method(f, lb, ub, prec);
}

inline double Legendre(int l, int m, double x)
{
    if(x*x>1.0)
        smart_throw("LEGENDRE: real number must be in [-1,1]");

    if(l<0)
        return Legendre(-l-1,m,x);

    if(m<0)
        return pow(-1.0,-m)*double(factorial(l+m))*Legendre(l,-m,x)/double(factorial(l-m));

    if(m>l)
        return 0.0;        

    /* l>=m>=0 */
    double out=1.0; // default value (l=0)

    if(l==1){
        if(m==1){
            out=-sqrt(1.0-x*x);
        }
        if(m==0){
            out=x;
        }       
    }

    if(l>1){
        if(l==m){
            out=-double(2*l-1)*sqrt(1.0-x*x)*Legendre(l-1,l-1,x);
        }
        else if(m==l-1){
            out=double(2*m+1)*x*Legendre(m,m,x);
        }        
        else{
            out=(double(2*l-1)*x*Legendre(l-1,m,x)-double(l-1+m)*Legendre(l-2,m,x))/double(l-m); 
        }        
    }

    return out;
}


inline double Legendre_derivative(int l, int m, double x)
{
    if(x*x>=1.0)
        smart_throw("LEGENDRE: real number must be in ]-1,1[");  

    return (double(l)*x*Legendre(l,m,x)-double(l+m)*Legendre(l-1,m,x))/(x*x-1.0); 
}



}

#endif // SMARTMATH_INLINEFUNCTIONS_H
