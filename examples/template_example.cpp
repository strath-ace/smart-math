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

    std::vector<double> times(4);
    std::vector<double> func1 = times;
    for(unsigned int k=0;k<times.size();k++){
    	times[k] = double(k+1);
		func1[k] = pow(times[k],times.size()-1);
    }
    double t = -2.0;
   	std::cout << Lagrange1d(times, func1, t) << std::endl;

}


