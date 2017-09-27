/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
--------- Author: Annalisa Riccardi ----------------------------------
*/


#include <iostream>
#include <Eigen/Dense>

#include "../include/smartmath.h"

using namespace std;

double H(const std::vector<double> vec){

	return 0.5 * vec[1] * vec[1] - cos(vec[0]);
}

int main(){

    cout << "Welcome to SMART-MATH!" << endl;
    
}
