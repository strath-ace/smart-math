/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2017 University of Strathclyde and Authors-------
-------- e-mail: romain.serra@strath.ac.uk ---------------------------
--------- Author: Romain Serra ---------------------------------------
*/

#ifndef SMARTMATH_FOREST_MIXEDVAR_H
#define SMARTMATH_FOREST_MIXEDVAR_H

#include "symplectic_mixedvar.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The forest_mixedvar class is an adaptation for mixed variables of the Forest integrator.
         *
         * The forest_mixedvar class is an adaptation for mixed variables of a 4th order integrator by Forest (1987).
         */
        template < class T >
        class forest_mixedvar: public symplectic_mixedvar<T>
        {

        protected:
            using symplectic_mixedvar<T>::m_ham;
            using symplectic_mixedvar<T>::m_stages;
            using symplectic_mixedvar<T>::m_c; 
            using symplectic_mixedvar<T>::m_d;              

        public:

            /**
             * @brief forest_mixedvar constructor
             *
             * The constructor initializes a pointer to the dynamics to integrate and precomputes a parameter necessary for integration
             * @param dyn Hamiltonian system to integrate
             */
            forest_mixedvar(const dynamics::hamiltonian_mixedvar<T> *dyn) : symplectic_mixedvar<T>("Forest integrator with mixed variables", dyn, 4){
                
                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("FOREST_MIXEDVAR: symplectic integrator cannot operate on non-separable Hamiltonian");

                double beta = pow(2.0, 1.0 / 3.0);

                std::vector<double> c(m_stages), d(m_stages);
                c[0] = 0.5 / (2.0 - beta);
                d[0] = 1.0 / (2.0 - beta); 
                c[1] = 0.5 * (1.0 - beta) / (2.0 - beta);
                d[1] = - beta / (2.0 - beta); 
                c[2] = 0.5 * (1.0 - beta) / (2.0 - beta);
                d[2] = 1.0 / (2.0 - beta); 
                c[3] = 0.5 / (2.0 - beta);
                d[3] = 0.0;

                m_c = c;
                m_d = d;                

            }

            /**
             * @brief ~forest_mixedvar deconstructor
             */
            ~forest_mixedvar(){}

        };

    }
}

#endif // SMARTMATH_FOREST_MIXEDVAR_H
