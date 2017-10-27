/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2017 University of Strathclyde and Authors-------
-------- e-mail: romain.serra@strath.ac.uk ---------------------------
--------- Author: Romain Serra ---------------------------------------
*/

#ifndef SMARTMATH_LEAPFROG_H
#define SMARTMATH_LEAPFROG_H

#include "base_symplectic.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The leapfrog class is an instantiation of symplectic integrators with order 2.
         *
         * The leapfrog class is a 2nd order instantiation of symplectic integrators. It has two versions: kick-drift-kick and drift-kick-drift.
         */
        template < class T >
        class leapfrog: public base_symplectic<T>
        {

        protected:
            using base_symplectic<T>::m_dyn;
            using base_symplectic<T>::m_stages;
            using base_symplectic<T>::m_c; 
            using base_symplectic<T>::m_d; 
            /**
             * @brief m_flag boolean to chose between drift-kick and kick-drift versions
             */          
            bool m_flag;

        public:

            /**
             * @brief leapfrog constructor
             *
             * The constructor initializes a pointer to the dynamics to integrate and a flag to decide on the integration scheme to use
             * @param dyn Hamiltonian system to integrate
             * @param flag boolean to know what algorithm to use (kick-drift-kick or drift-kick-drift)
             */
            leapfrog(const dynamics::base_hamiltonian<T> *dyn, const bool &flag) : base_symplectic<T>("leapfrog integrator", dyn, 2), m_flag(flag){
                
                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("LEAPFROG: symplectic integrator cannot operate on non-separable Hamiltonian");

                /* computation integration coefficients depending on chosen algorithm */
                std::vector<double> c(m_stages), d(m_stages);
                if(flag)
                { //drift-kick-drift
                    c[0] = 1.0 / 2.0;
                    d[0] = 1.0;
                    c[1] = 1.0 / 2.0;
                    d[1] = 0.0;
                }
                else
                { //kick-drift-kick
                    c[0] = 0.0;
                    d[0] = 1.0 / 2.0;
                    c[1] = 1.0;
                    d[1] = 1.0 / 2.0; 
                }
                m_c = c;
                m_d = d;
            }

            /**
             * @brief ~leapfrog deconstructor
             */
            ~leapfrog(){}

        };

    }
}

#endif // SMARTMATH_LEAPFROG_H
