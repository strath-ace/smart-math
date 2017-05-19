/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: romain.serra@strath.ac.uk -----------------------
----------------- Author: Romain Serra -------------------------------
*/


#ifndef SMARTMATH_HAMILTONIAN_MIXEDVAR_H
#define SMARTMATH_HAMILTONIAN_MIXEDVAR_H

#include "base_hamiltonian.h"
#include "../exception.h"

namespace smartmath
{
    namespace dynamics {

        template < class T >
        class hamiltonian_mixedvar: public base_hamiltonian<T>
        {

        protected:
        	using smartmath::dynamics::base_hamiltonian<T>::m_dim;

        public:
            hamiltonian_mixedvar(const std::string &name, const int &dim, const bool &separable): base_hamiltonian<T>(name, dim, separable){}

            ~hamiltonian_mixedvar(){}

            using smartmath::dynamics::base_hamiltonian<T>::evaluate;

            virtual int DHq(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const = 0;

            virtual int DHp(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const = 0;

            virtual int DHq2(const double &t, const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &dH2) const = 0;

            virtual int DHp2(const double &t, const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &dH2) const = 0;            

            int evaluate2(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const{

            	/* sanity checks */
				if(state.size() != 2 * m_dim)
					smartmath_throw("EVALUATE: the Hamiltonian state must have a consistent dimension");
				if(state.size() != dstate.size())
					smartmath_throw("EVALUATE: the derivative of the state must have same length");				

				dstate.clear();

				std::vector<T> dHp2, dHq2, q2, p2;
				for(int k = 0; k < m_dim; k++)
					{
						q2.push_back(state[k]);
						p2.push_back(state[k + m_dim]);
					}

				DHq2(t, q2, p2, dHq2);
				DHp2(t, q2, p2, dHp2);

				for(int k = 0; k < m_dim; k++)
					{
						dstate.push_back(dHp2[k]);
					}
				for(int k = 0; k < m_dim; k++)
					{
						dstate.push_back(-dHq2[k]);
					}

            	return 0;
            }

            virtual int conversion(const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &q2, std::vector<T> &p2) const = 0; 

            virtual int conversion2(const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &q, std::vector<T> &p) const = 0; 

        };

    }
}

#endif // SMARTMATH_HAMILTONIAN_MIXEDVAR_H
