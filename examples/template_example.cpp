/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors---------
------------ e-mail: annalisa.riccardi@strath.ac.uk -------------------------
------------- Author: Annalisa Riccardi -------------------------------------
*/


#include <iostream>

#include "../include/smartmath.h"

using namespace std;
using namespace smartmath;

int main(){

    cout << "Welcome to SMART-MATH!" << endl;

    int n = 5;
    int k = 3;
    cout << "\nCombination between " << n << " and " << k << ": " << combination(n,k) << endl;

}
