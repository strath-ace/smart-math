/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2017 University of Strathclyde and Authors-------
-------- e-mail: romain.serra@strath.ac.uk ---------------------------
--------- Author: Romain Serra ---------------------------------------
*/

#ifndef SMARTMATH_YOSHIDA6_MIXEDVAR_H
#define SMARTMATH_YOSHIDA6_MIXEDVAR_H

#include "symplectic_mixedvar.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The %yoshida6_mixedvar class is an adaptation for mixed variables of the Yoshida 6 integrator.
         *
         * The %yoshida6_mixedvar class is an adaptation for mixed variables of a 6th order integrator by Yoshida.
         */
        template < class T >
        class yoshida6_mixedvar: public symplectic_mixedvar<T>
        {

        protected:
            using symplectic_mixedvar<T>::m_ham;
            using symplectic_mixedvar<T>::m_stages;
            using symplectic_mixedvar<T>::m_c; 
            using symplectic_mixedvar<T>::m_d;   

        public:

            /**
             * @brief yoshida6_mixedvar constructor
             *
             * The constructor initializes a pointer to the dynamics to integrate and precomputes a parameter necessary for integration
             * @param dyn Hamiltonian system to integrate
             */
            yoshida6_mixedvar(const dynamics::hamiltonian_mixedvar<T> *dyn) : symplectic_mixedvar<T>("Yoshida 6 integrator with mixed variables", dyn, 8){
                
                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("YOSHIDA6_MIXEDVAR: symplectic integrator cannot operate on non-separable Hamiltonian");

                std::vector<double> c(m_stages), d(m_stages);
                c[0] = 0.392256805238780;
                d[0] = 0.784513610477560;
                c[1] = 0.510043411918458;
                d[1] = 0.235573213359357;
                c[2] = -0.471053385409758;
                d[2] = -0.117767998417887e1;
                c[3] = 0.687531682525198e-1;
                d[3] = 0.131518632068391e1;
                c[4] = 0.687531682525198e-1;
                d[4] = -0.117767998417887e1;
                c[5] = -0.471053385409758;
                d[5] = 0.235573213359357;
                c[6] = 0.510043411918458;
                d[6] = 0.784513610477560;
                c[7] = 0.392256805238780;
                d[7] = 0.0;

                m_c = c;
                m_d = d; 

            }

            /**
             * @brief ~yoshida6_mixedvar deconstructor
             */
            ~yoshida6_mixedvar(){}

        };

    }
}

#endif // SMARTMATH_YOSHIDA6_MIXEDVAR_H
