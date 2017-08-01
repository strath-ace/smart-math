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
    cout << endl;

    cout << "Following is a simple example (the mathematical pendulum) on the advantage of using symplectic integrators for Hamiltonian systems" << endl;

    /* Creating the dynamics */
    smartmath::dynamics::pendulum<double> *dyn = new smartmath::dynamics::pendulum<double>();

	/* Creating integrators */
    smartmath::integrator::ABM<double> prop1(dyn, 6); // predictor-corrector of order 6 
    smartmath::integrator::yoshida6<double> prop2(dyn); // symplectic integrator of order 6 (note that the dynamics has to be coded in an Hamiltonian form)

    /* Setting initial conditions */
    std::vector<double> x(2), x1 = x, x2 = x;
    x[0] = 0.1;
    x[1] = 0.01;

    cout << "Initial conditions are: (" << x[0] << ", " << x[1] << "). The corresponding value for the Hamiltonian H is: " << H(x) << endl;

    double dt = 1.0e5; // time of propagation
    int steps = 1e6; // number of integration steps

    prop1.integrate(0.0, dt, steps, x, x1); // propagation with non-symplectic integrator
    prop2.integrate(0.0, dt, steps, x, x2); // propagation with symplectic integrator

	cout << "Values for H after an elapsed time of " << dt << " are: " << H(x2) << " with symplectic integration and " << H(x1) << " without" << endl;
    cout << endl;

    delete dyn;

}
