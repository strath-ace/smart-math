/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_LEAPFROG_H
#define SMARTMATH_LEAPFROG_H

#include "base_symplectic.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The base_rungekutta class is a template abstract class. Any fixed-step Runge-Kutta algorithm added to the toolbox needs to inherit from it and implement the method integration_step()
         *
         * The base_rungekutta class is a template abstract class. Any fixed-step Runge-Kutta algorithm added to the toolbox needs to inherit from it and implement the method that performs on integration step between to given times given the initial state
         */
        template < class T >
        class leapfrog: public base_symplectic<T>
        {

        protected:
            using base_symplectic<T>::m_dyn;
            using base_symplectic<T>::m_order;
            using base_symplectic<T>::m_c;
            using base_symplectic<T>::m_d;
            bool m_flag;

        public:

            /**
             * @brief base_rungekutta constructor
             *
             * The constructor initialize the name of the integrator and a pointer to the dynamical system to be integrated
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             */
            leapfrog(const dynamics::base_hamiltonian<T> *dyn, const bool &flag) : base_symplectic<T>("leapfrog integrator", dyn, 2), m_flag(flag){
                
                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("LEAPFROG: symplectic integrator cannot operate on non-separable Hamiltonian");

                std::vector<double> c(m_order), d(m_order);

                if(flag){
                    c[0] = 1.0 / 2.0;
                    c[1] = 1.0 / 2.0;
                    d[0] = 1.0;
                    d[1] = 0.0;
                }
                else{
                    c[0] = 0.0;
                    c[1] = 1.0;
                    d[0] = 1.0 / 2.0;
                    d[1] = 1.0 / 2.0; 
                }

                m_c = c;
                m_d = d;
            }

            /**
             * @brief ~base_rungekutta deconstructor
             */
            ~leapfrog(){}

        };

    }
}

#endif // SMARTMATH_LEAPFROG_H
