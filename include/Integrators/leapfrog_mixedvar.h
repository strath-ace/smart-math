/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2017 University of Strathclyde and Authors-------
-------- e-mail: romain.serra@strath.ac.uk ---------------------------
--------- Author: Romain Serra ---------------------------------------
*/

#ifndef SMARTMATH_LEAPFROG_MIXEDVAR_H
#define SMARTMATH_LEAPFROG_MIXEDVAR_H

#include "symplectic_mixedvar.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The leapfrog_mixedvar class is an adaptation of the leapfrog algorithm for mixed variables.
         *
         * The leapfrog class is an adaptation of the leapfrog algorithm for mixed variables. It has two versions: kick-drift-kick and drift-kick-drift.
         */
        template < class T >
        class leapfrog_mixedvar: public symplectic_mixedvar<T>
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
             * @brief leapfrog_mixedvar constructor
             *
             * The constructor initializes a pointer to the dynamics to integrate and a flag to decide on the integration scheme to use
             * @param dyn Hamiltonian system to integrate
             * @param flag boolean to know what algorithm to use (kick-drift-kick or drift-kick-drift)
             */
            leapfrog_mixedvar(const dynamics::hamiltonian_mixedvar<T> *dyn, const bool &flag) : symplectic_mixedvar<T>("leapfrog integrator with mixed variables", dyn, 3), m_flag(flag){
                
                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("LEAPFROG_MIXEDVAR: symplectic integrator cannot operate on non-separable Hamiltonian");

                std::vector<double> c(m_stages), d(m_stages);
                if(flag)
                { //drift-kick-drift
                    c[0] = 1.0 / 2.0;
                    d[0] = 0.0;
                    c[1] = 0.0;
                    d[1] = 1.0;
                    c[2] = 1.0 / 2.0;
                    d[2] = 0.0;                                       
                }
                else
                { //kick-drift-kick
                    c[0] = 0.0;
                    d[0] = 1.0 / 2.0;
                    c[1] = 1.0;
                    d[1] = 0.0;
                    c[2] = 0.0;
                    d[2] = 1.0 / 2.0; 
                }
                m_c = c;
                m_d = d;
                
            }

            /**
             * @brief ~leapfrog_mixedvar deconstructor
             */
            ~leapfrog_mixedvar(){}

        };

    }
}

#endif // SMARTMATH_LEAPFROG_MIXEDVAR_H
