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

        /**
         * @brief The base_hamiltonian class is a template abstract class. Any Hamiltonian system added to the toolbox needs to inherit from it
         *
         * The base_hamiltonian class is a template abstract class. Any Hamiltonian system added to the toolbox needs to inherit from it
         * The canonical variables are called q and p. 
         */
        template < class T >
        class base_hamiltonian: public base_dynamics<T>
        {

        protected:
            using base_dynamics<T>::m_name;
            /**
             * @brief m_dim half-dimension of Hamiltonian system
             */             
			unsigned int m_dim;
            /**
             * @brief m_separable flag equal to 1 if the system is separable i.e. if the Hamiltonian is of the type V(p) + W(q), 0 otherwise
             */             
			bool m_separable; 

        public:
            /**
             * @brief base_hamiltonian constructor
             *
             * The constructor initializes the name of the Hamiltonian dynamics, its half-dimension and a flag about its separability
             * @param name integrator name
             * @param dim half-order of the Hamiltonian system
             * @param separable boolean precising whether the system is separable or not
             */
            base_hamiltonian(const std::string &name, const unsigned int &dim, const bool &separable = false): base_dynamics<T>(name), m_dim(dim), m_separable(separable){}

            /**
             * @brief ~base_hamiltonian deconstructor
             */
            ~base_hamiltonian(){}

            /**
             * @brief get_dim performs get the dimension of the Hamiltonian system
             *
             * The method returns the dimension of the implemented Hamiltonian system
             * @return m_dim half-order of the system
             */
            unsigned int get_dim() const{
            	return m_dim;
            }

            /**
             * @brief is_separable performs get the dimension of the Hamiltonian system
             *
             * The method returns a boolean that is true if the system is seperable, false otherwise
             * @return m_separable boolean about the separability of the Hamiltonian 
             */
            bool is_separable() const{
            	return m_separable;
            }            

            /**
             * @brief evaluate differential equations of the implemented Hamiltonian system
             * @param[in] t time in scaled units
             * @param[in] state vector in scaled units
             * @param[out] dstate derivative in scaled units
             * @return exit flag (0=success)
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const{

            	/* sanity checks */
				if(state.size() != 2 * m_dim)
					smartmath_throw("EVALUATE: the Hamiltonian state must have a consistent dimension");

				dstate.clear();
				/* reconstituting the canonical variables q and p from the state vector */
				std::vector<T> q, p;
				for(unsigned int k = 0; k < m_dim; k++)
				{
					q.push_back(state[k]);
					p.push_back(state[k + m_dim]);
				}
                
				/* computing the partial derivatives of the Hamiltonian w.r.t. q and p */ 
                std::vector<T> dHp = p, dHq = q;
				DHq(t, q, p, dHq);       
				DHp(t, q, p, dHp);

				/* reconstituting the state derivative */
				for(unsigned int k = 0; k < m_dim; k++)
					dstate.push_back(dHp[k]);
				for(unsigned int k = 0; k < m_dim; k++)
					dstate.push_back(-dHq[k]);

            	return 0;
            }

            /**
             * @brief DHq computes the partial derivative of the Hamiltonian with respect to q
             *
             * The method computes the partial derivative of the Hamiltonian with respect to q
             * @param[in] t time in scaled units
             * @param[in] q vector in scaled units
             * @param[in] p vector in scaled units
             * @param[out] dH vector of partial derivatives of H w.r.t. the vector q
             * @return exit flag (0=success)
             */
			virtual int DHq(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const = 0;

            /**
             * @brief DHq computes the partial derivative of the Hamiltonian with respect to p
             *
             * The method computes the partial derivative of the Hamiltonian with respect to p
             * @param[in] t time in scaled units
             * @param[in] q vector in scaled units
             * @param[in] p vector in scaled units
             * @param[out] dH vector of partial derivatives of H w.r.t. the vector p
             * @return exit flag (0=success)
             */
			virtual int DHp(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const = 0;

        };

    }
}

#endif // SMARTMATH_HAMILTONIAN_H
