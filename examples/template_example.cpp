/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
--------- Author: Annalisa Riccardi ----------------------------------
*/


#include <iostream>

#include "../include/smartmath.h"

using namespace std;
using namespace smartmath;

int main(){

    cout << "Welcome to SMART-MATH!" << endl;

    //int n = 5;
    //int k = 3;
    //std::cout << "\nCombination: " << smartmath::combination(n,k) << std::endl;

    //double x=-0.5276;
    //int l=4;
    //int m=-3;
    //std::cout << Legendre(l,m,x) << std::endl;
    //std::cout << 105.0*x*pow(1.0-x*x,1.5)/5040.0 << std::endl; // P(4,-3,x)
    //std::cout << Legendre_derivative(l,m,x) << std::endl;

    /* Return as it is
    std::vector<double> times(4), func1 = times;
    std::vector<double> func2 = func1;
    std::vector< std::vector<double> > funcs(times.size());
    for(unsigned int k=0;k<times.size();k++){
    	times[k] = double(k+1);
		func1[k] = pow(times[k],times.size()-1);
		func2[k] = -pow(times[k],times.size()-1);
		funcs[k].push_back(func1[k]);
		funcs[k].push_back(func2[k]);
    }   
    double t = -2.0;
    std::cout << Lagrange1d(times, func1, t) << std::endl;
    std::vector<double> interpolated_vector = LagrangeNd(times, funcs, t);
   	std::cout << interpolated_vector[0] << ", " << interpolated_vector[1] << std::endl;
     */

    std::vector<double> roots = {1,2};
    std::vector<double> coeffs ;

    std::cout << "\nRoots are: \n" ;
    for( double &number : roots)
        cout << number << " " ;
    cout << "\n";

    coeffs = vieta_root2coef(roots);


    std::cout << "\nCoefficients are: \n" ;
    for( double &number : coeffs)
        cout << number << " " ;
    cout << "\n";

}


