/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2017 University of Strathclyde and Authors-------
-------- e-mail: romain.serra@strath.ac.uk ---------------------------
--------- Author: Romain Serra ---------------------------------------
*/

#ifndef SMARTMATH_YOSHIDA6_H
#define SMARTMATH_YOSHIDA6_H

#include "base_symplectic.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The %yoshida6 class is an instantiation of symplectic integrators with order 6.
         *
         * The %yoshida6 class is a 6th order instantiation of symplectic integrators from Yoshida (1990). 
         */
        template < class T >
        class yoshida6: public base_symplectic<T>
        {

        protected:
            using base_symplectic<T>::m_dyn;
            using base_symplectic<T>::m_stages;
            using base_symplectic<T>::m_c;
            using base_symplectic<T>::m_d;

        public:

            /**
             * @brief yoshida6 constructor
             *
             * The constructor initializes a pointer to the dynamics to integrate and precomputes the integration coefficients
             * @param dyn Hamiltonian system to integrate
             */
            yoshida6(const dynamics::base_hamiltonian<T> *dyn) : base_symplectic<T>("Forest scheme", dyn, 8){

                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("YOSHIDA6: symplectic integrator cannot operate on non-separable Hamiltonian");

                std::vector<double> c(m_stages), d(m_stages);

                c[0] = 0.392256805238780;
                c[7] = c[0];
                c[1] = 0.510043411918458;
                c[6] = c[1];
                c[2] = -0.471053385409758;
                c[5] = c[2];
                c[3] = 0.687531682525198e-1;
                c[4] = c[3];

                d[0] = 0.784513610477560;
                d[6] = d[0];
                d[1] = 0.235573213359357;
                d[5] = d[1];
                d[2] = -0.117767998417887e1;
                d[4] = d[2];
                d[3] = 0.131518632068391e1;
                d[7] = 0.0;

                m_c = c;
                m_d = d;

            }

            /**
             * @brief ~yoshida6 deconstructor
             */
            ~yoshida6(){}

        };

    }
}

#endif // SMARTMATH_YOSHIDA6_H
