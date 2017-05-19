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

        template < class T >
        class hamiltonian_momentum: public base_hamiltonian<T>
        {

        protected:
        	using smartmath::dynamics::base_hamiltonian<T>::m_dim;

        public:
            hamiltonian_momentum(const std::string &name, const int &dim, const bool &separable): base_hamiltonian<T>(name, dim, separable){}

            ~hamiltonian_momentum(){}

            using smartmath::dynamics::base_hamiltonian<T>::evaluate;

            virtual int DHq(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const = 0;//{

			//     /* sanity checks */
			//     if(q.size() != m_dim)
			//         smartmath_exception("DHQ: the position must have dimension 3");
			//     if(p.size() != m_dim)
			//         smartmath_exception("DHQ: the momentum must have dimension 3");   

			//     dH.clear();
			//     for(int i = 0; i < m_dim; i++)
			//         dH.push_back(0.0 * q[i]); // dummy force

			//     return 0;

			// };

			int DHp(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const{

			    /* sanity checks */
			    if(q.size() != m_dim)
			        smartmath_exception("DHP: the position must have dimension 3");
			    if(p.size() != m_dim)
			        smartmath_exception("DHP: the momentum must have dimension 3");   

			    dH.clear();
			    for(int i = 0; i < m_dim; i++)
			        dH.push_back(p[i]);

			    return 0;

			};

        };

    }
}

#endif // SMARTMATH_HAMILTONIAN_MOMENTUM_H
