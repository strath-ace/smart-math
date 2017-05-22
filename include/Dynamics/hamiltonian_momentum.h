/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: romain.serra@strath.ac.uk -----------------------
----------------- Author: Romain Serra -------------------------------
*/


#ifndef SMARTMATH_HAMILTONIAN_MOMENTUM_H
#define SMARTMATH_HAMILTONIAN_MOMENTUM_H

#include "base_hamiltonian.h"
#include "../exception.h"

namespace smartmath
{
    namespace dynamics {
        /**
         * @brief The hamiltonian_momentum class is a template abstract class. Any Hamiltonian system using position and momenta as canonical variables needs to inherit from it
         *
         * The hamiltonian_momentum class is a template abstract class. Any Hamiltonian system using position and momenta as canonical variables needs to inherit from it
         */
        template < class T >
        class hamiltonian_momentum: public base_hamiltonian<T>
        {

        protected:
        	using smartmath::dynamics::base_hamiltonian<T>::m_dim;

        public:
             /**
             * @brief hamiltonian_momentum constructor
             *
             * The constructor initialize the name of the Hamiltonian dynamics, its half-dimension and a flag about its separability
             * @param name integrator name
             * @param dim half-order of the Hamiltonian system
             * @param separable boolean precising whether the system is separable or not
             */
            hamiltonian_momentum(const std::string &name, const int &dim, const bool &separable): base_hamiltonian<T>(name, dim, separable){}

            /**
             * @brief ~hamiltonian_momentum deconstructor
             */
            ~hamiltonian_momentum(){}

            /**
             * @brief DHp computes the partial derivative of the Hamiltonian with respect to the momenta
             *
             * The method computes the partial derivative of the Hamiltonian with respect to p the vector of momentum
             * @param[in] time in scaled units
             * @param[in] q position vector in scaled units
             * @param[in] p momenta vector in scaled units
             * @param[out] dH vector of partial derivatives of H w.r.t. the vector q
             * @return exit flag (0=success)
             */
			int DHp(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const{

			    if(p.size() != m_dim)
			        smartmath_exception("DHP: the momentum must have the correct dimension");

			    dH.clear();
			    for(int i = 0; i < m_dim; i++)
			        dH.push_back(p[i]);

			    return 0;
			};

        };

    }
}

#endif // SMARTMATH_HAMILTONIAN_MOMENTUM_H
