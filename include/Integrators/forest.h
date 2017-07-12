/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_FOREST_H
#define SMARTMATH_FOREST_H

#include "base_symplectic.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The forest class is an instantiation of symplectic integrators with order 4.
         *
         * The forest class is a 4nd order instantiation of symplectic integrators from Forest (1987). 
         */
        template < class T >
        class forest: public base_symplectic<T>
        {

        protected:
            using base_symplectic<T>::m_dyn;
            using base_symplectic<T>::m_stages;
            using base_symplectic<T>::m_c;
            using base_symplectic<T>::m_d;

        public:

            /**
             * @brief forest constructor
             *
             * The constructor initializes a pointer to the dynamics to integrate and precomputes the integration coefficients
             * @param dyn Hamiltonian system to integrate
             */
            forest(const dynamics::base_hamiltonian<T> *dyn) : base_symplectic<T>("Forest scheme", dyn, 4){

                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("FOREST: symplectic integrator cannot operate on non-separable Hamiltonian");

                double beta = pow(2.0, 1.0 / 3.0);
                std::vector<double> c(m_stages), d(m_stages);

                c[0] = 0.5 / (2.0 - beta);
                c[3] = c[0];
                c[1] = 0.5 * (1.0 - beta) / (2.0 - beta);
                c[2] = c[1];

                d[0] = 1.0 / (2.0 - beta);
                d[2] = d[0];
                d[1] = -beta / (2.0 - beta);
                d[3] = 0.0;

                m_c = c;
                m_d = d;

            }

            /**
             * @brief ~forest deconstructor
             */
            ~forest(){}

        };

    }
}

#endif // SMARTMATH_FOREST_H
