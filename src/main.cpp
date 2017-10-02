/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
--------- Author: Annalisa Riccardi ----------------------------------
*/

#include "../include/smartmath.h"

using namespace std;

int main(){

	cout << "This is an example to illustrate integration with mixed variables" << endl;

	/* Creating the dynamics */
	smartmath::dynamics::spring<double> *dyn = new smartmath::dynamics::spring<double>();

	/* Creating integrator */
	smartmath::integrator::leapfrog<double> prop(dyn, true);

	/* Setting initial conditions */
	std::vector<double> x(2), xf;
	x[0] = 0.1; // angle
	x[1] = 0.01; // angular velocity

	/* Integration */
	double t_0 = 5.0;
	double t_f = 10.0; // projected final time
	prop.integrate(t_0, t_f, 1e2, x, xf); 

	cout << "The analytical solution to the harmonic oscillator is compared to the numerical one:" << endl;

	cout << sqrt(x[0] * x[0] + x[1] * x[1]) * sin(t_f - (t_0 - atan2(x[0], x[1]))) << " VS " << xf[0] << endl;

	delete dyn;
}
