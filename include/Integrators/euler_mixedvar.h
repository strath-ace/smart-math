/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_EULER_MIXEDVAR_H
#define SMARTMATH_EULER_MIXEDVAR_H

#include "symplectic_mixedvar.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The euler_mixedvar class is an adaptation of the leapfrog algorithm for mixed variables.
         *
         * The euler_mixedvar class is an adaptation of the symplectic Euler algorithm for mixed variables. It has two versions: kick-drift and drift-kick.
         */
        template < class T >
        class euler_mixedvar: public symplectic_mixedvar<T>
        {

        protected:
            using symplectic_mixedvar<T>::m_ham;
            using symplectic_mixedvar<T>::m_stages;
            using symplectic_mixedvar<T>::m_c; 
            using symplectic_mixedvar<T>::m_d;             
            /**
             * @brief m_flag boolean to chose between drift-kick and kick-drift versions
             */             
            bool m_flag;

        public:

            /**
             * @brief euler_mixedvar constructor
             *
             * The constructor initialize a pointer to the dynamics to integrate and a flag to decide on the integration scheme to use
             * @param dyn Hamiltonian system to integrate
             * @param flag boolean to know what algorithm to use (kick-drift-kick or drift-kick-drift)
             */
            euler_mixedvar(const dynamics::hamiltonian_mixedvar<T> *dyn, const bool &flag) : symplectic_mixedvar<T>("symplectic Euler integrator with mixed variables", dyn, 2), m_flag(flag){
                
                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("EULER_MIXEDVAR: symplectic integrator cannot operate on non-separable Hamiltonian");

                std::vector<double> c(m_stages), d(m_stages);
                if(flag)
                { //drift-kick
                    c[0] = 1.0;
                    d[0] = 0.0;
                    c[1] = 0.0;
                    d[1] = 1.0;                  
                }
                else
                { //kick-drift                 
                    c[0] = 0.0;
                    d[0] = 1.0;
                    c[1] = 1.0;
                    d[1] = 0.0;
                }
                m_c = c;
                m_d = d;
            }

            /**
             * @brief ~euler_mixedvar deconstructor
             */
            ~euler_mixedvar(){}

        };

    }
}

#endif // SMARTMATH_EULER_MIXEDVAR_H
