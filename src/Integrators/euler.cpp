/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#include "../../include/Integrators/euler.h"


using namespace smartmath;
using namespace smartmath::integrator;

template < class T >
euler<T>::euler(const dynamics::base_dynamics<T> *dyn) : base_integrator<T>("Forward Euler integration scheme", dyn)
{
}

template < class T >
euler<T>::~euler(){

}

template < class T >
int euler<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

	xfinal.clear();

	std::vector<T> dx(x0);
	std::vector<T> x(x0);

	double h = (tend-ti)/nsteps;

    for(int i=0; i<nsteps; i++){
		m_dyn->evaluate(ti+i*h, x, dx);
		for(size_t j=0; j<x.size(); j++){
			x[j] += h*dx[j];
		}
	}

	for(int i=0; i<x0.size(); i++)
	    xfinal.push_back(x[i]);

	return 0;
}


template class euler<double>;
template class euler<float>;
template class euler<long double>;
#ifdef ENABLE_SMARTUQ
template class euler<smartuq::polynomial::chebyshev_polynomial<double> >;
template class euler<smartuq::polynomial::chebyshev_polynomial<float> >;
template class euler<smartuq::polynomial::chebyshev_polynomial<long double> >;
template class euler<smartuq::polynomial::taylor_polynomial<double> >;
template class euler<smartuq::polynomial::taylor_polynomial<float> >;
template class euler<smartuq::polynomial::taylor_polynomial<long double> >;
#endif
