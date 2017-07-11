/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_EULER_SYMPLECTIC_H
#define SMARTMATH_EULER_SYMPLECTIC_H

#include "base_symplectic.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The euler_symplectic class is a instantiation of symplectic integrators with order 1.
         *
         * The euler_symplectic class is a 1st order instantiation of symplectic integrators. It has two versions: kick-drift and drift-kick.
         */
        template < class T >
        class euler_symplectic: public base_symplectic<T>
        {

        protected:
            using base_symplectic<T>::m_dyn;
            using base_symplectic<T>::m_stages;
            using base_symplectic<T>::m_c; // drift coefficients
            using base_symplectic<T>::m_d; // kick coefficients
            /**
             * @brief m_flag boolean to chose between drift-kick and kick-drift versions
             */ 
            bool m_flag;

        public:

            /**
             * @brief euler_symplectic constructor
             *
             * The constructor initialize a pointer to the dynamics to integrate and a flag to decide on the integration scheme to use
             * @param dyn Hamiltonian system to integrate
             * @param flag boolean to know what algorithm to use (drift-kick or kick-drift)
             */
            euler_symplectic(const dynamics::base_hamiltonian<T> *dyn, const bool &flag) : base_symplectic<T>("symplectic version of explicit Euler integrator", dyn, 1), m_flag(flag){
                
                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("EULER_SYMPLECTIC: symplectic integrator cannot operate on non-separable Hamiltonian");

                /* computation integration coefficients depending on chosen algorithm */
                std::vector<double> c(m_stages), d(m_stages);
                if(flag){ //drift-kick
                    c[0] = 1.0;
                    d[0] = 1.0;
                }
                else{ //kick-drift
                    ++m_stages;
                    c[0] = 0.0;
                    d[0] = 1.0;
                    c.push_back(1.0);
                    d.push_back(0.0);
                }
                m_c = c;
                m_d = d;
            }

            /**
             * @brief ~euler_symplectic deconstructor
             */
            ~euler_symplectic(){}

        };

    }
}

#endif // SMARTMATH_EULER_SYMPLECTIC_H
