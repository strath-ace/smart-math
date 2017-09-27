/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: romain.serra@strath.ac.uk -----------------------
----------------- Author: Romain Serra -------------------------------
*/

#ifndef SMARTMATH_RK4_H
#define SMARTMATH_RK4_H

#include "base_rungekutta.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The Runge Kutta 4 integrator scheme
         *
         * The class model the Runge Kutta fourth order integration scheme
         */
        template < class T >
        class rk4: public base_rungekutta<T>
        {

        private:
            using base_rungekutta<T>::m_name;
            using base_rungekutta<T>::m_dyn;
            using base_rungekutta<T>::m_stages;
            using base_rungekutta<T>::m_coeT;
            using base_rungekutta<T>::m_coeK;
            using base_rungekutta<T>::m_coeX;

        public:

            using base_rungekutta<T>::integrate;

            /**
             * @brief rk4 constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn pointer to dynamical system to be integrated
             */
            rk4(const dynamics::base_dynamics<T> *dyn) : base_rungekutta<T>("Runge Kutta 4 fixed time-step", dyn){

                m_stages = 4;

                std::vector<double> coe(m_stages);
                coe[0] = 0.0;
                coe[1] = 1.0 / 2.0;
                coe[2] = 1.0 / 2.0;
                coe[3] = 1.0;
                m_coeT = coe;

                coe[0] = 1.0 / 6.0;
                coe[1] = 1.0 / 3.0;
                coe[2] = 1.0 / 3.0;
                coe[3] = 1.0 / 6.0;
                m_coeX = coe;   

                std::vector<double> coe2(m_stages - 1);
                coe2[0] = 1.0 / 2.0;
                coe2[1] = 1.0 / 2.0;
                coe2[2] = 1.0;
                m_coeK = coe2;    
            }

            /**
              * @brief ~rk4 deconstructor
              */
            ~rk4(){}

        };

    }
}

#endif // SMARTMATH_RK4_H
