/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2017 University of Strathclyde and Authors-------
-------- e-mail: romain.serra@strath.ac.uk ---------------------------
--------- Author: Romain Serra ---------------------------------------
*/

#ifndef SMARTMATH_EULER_H
#define SMARTMATH_EULER_H

#include "base_rungekutta.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The Euler explicit integration scheme
         *
         * The class model the Euler explicit integration scheme
         */
        template < class T >
        class euler: public base_rungekutta<T>
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
            using base_rungekutta<T>::integrate_eigen;

            /**
             * @brief euler constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             */
            euler(const dynamics::base_dynamics<T> *dyn): base_rungekutta<T>("Explicit Euler integration scheme", dyn){

                m_stages = 1;

                std::vector<double> coe(m_stages);
                coe[0] = 0.0;
                m_coeT = coe;
                coe[0] = 1.0;
                m_coeX = coe;   

            }

            /**
              * @brief ~euler deconstructor
              */
            ~euler(){}

        };

    }
}

#endif // SMARTMATH_EULER_H
