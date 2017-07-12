/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: romain.serra@strath.ac.uk -----------------------
----------------- Author: Romain Serra -------------------------------
*/

#ifndef SMARTMATH_MIDPOINT_H
#define SMARTMATH_MIDPOINT_H

#include "base_rungekutta.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The midpoint explicit integration scheme
         *
         * The class models the midpoint explicit integration scheme
         */
        template < class T >
        class midpoint: public base_rungekutta<T>
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
             * @brief midpoint constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn pointer to the dynamical system to be integrated
             */
            midpoint(const dynamics::base_dynamics<T> *dyn): base_rungekutta<T>("Explicit midpoint integration scheme", dyn){

                m_stages = 2;

                std::vector<double> coe(m_stages);
                coe[0] = 0.0;
                coe[1] = 1.0 / 2.0;
                m_coeT = coe;

                coe[0] = 0.0;
                coe[1] = 1.0;
                m_coeX = coe;   

                std::vector<double> coe2(m_stages - 1);
                coe2[0] = 1.0 / 2.0;
                m_coeK = coe2;  
            }

            /**
              * @brief ~midpoint deconstructor
              */
            ~midpoint(){}

        };

    }
}

#endif // SMARTMATH_MIDPOINT_H
