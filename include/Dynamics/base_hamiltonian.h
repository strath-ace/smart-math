/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: romain.serra@strath.ac.uk -----------------------
----------------- Author: Romain Serra -------------------------------
*/


#ifndef SMARTMATH_BASE_HAMILTONIAN_H
#define SMARTMATH_BASE_HAMILTONIAN_H

#include "base_dynamics.h"
#include "../exception.h"

namespace smartmath
{
    namespace dynamics {

        template < class T >
        class base_hamiltonian: public base_dynamics<T>
        {

        protected:
            using base_dynamics<T>::m_name;
			int m_dim;
			bool m_separable;

        public:
            base_hamiltonian(const std::string &name, const int &dim, const bool &separable): base_dynamics<T>(name), m_dim(dim), m_separable(separable){}

            ~base_hamiltonian(){}

            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const{

            	/* sanity checks */
				if(state.size() != 2 * m_dim)
					smartmath_throw("EVALUATE: the Hamiltonian state must have a consistent dimension");
				if(state.size() != dstate.size())
					smartmath_throw("EVALUATE: the derivative of the state must have same length");				

				dstate.clear();

				std::vector<T> dHp, dHq, q, p;
				for(int k = 0; k < m_dim; k++)
					{
						q.push_back(state[k]);
						p.push_back(state[k + m_dim]);
					}

				DHq(t, q, p, dHq);
				DHp(t, q, p, dHp);

				for(int k = 0; k < m_dim; k++)
					{
						dstate.push_back(dHp[k]);
					}
				for(int k = 0; k < m_dim; k++)
					{
						dstate.push_back(-dHq[k]);
					}

            	return 0;
            }

			virtual int Hamiltonian(const double &t, const std::vector<T> &q, const std::vector<T> &p, T &H) const = 0;

			virtual int DHq(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const = 0;

			virtual int DHp(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const = 0;

        };

    }
}

#endif // SMARTMATH_HAMILTONIAN_H
