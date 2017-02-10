/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_BASE_MULTISTEP_H
#define SMARTMATH_BASE_MULTISTEP_H

#include "base_multistep.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The base_multistep class is a template abstract class. Any fixed-step Runge-Kutta algorithm added to the toolbox needs to inherit from it and implement the method integration_step()
         *
         * The base_multistep class is a template abstract class. Any fixed-step Runge-Kutta algorithm added to the toolbox needs to inherit from it and implement the method that performs on integration step between to given times given the initial state 
         */
        template < class T >
        class base_multistep: public base_integrator<T>
        {

        protected:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            int m_order;

        public:

            using base_integrator<T>::integrate;

            /**
             * @brief base_multistep constructor
             *
             * The constructor initialize the name of the integrator and a pointer to the dynamical system to be integrated
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             */
            base_multistep(const std::string &name, const dynamics::base_dynamics<T> *dyn, const int &order) : base_integrator<T>(name, dyn), m_order(order){}

            /**
             * @brief ~base_multistep deconstructor
             */
            virtual ~base_multistep(){}

            /**
             * @brief initialize method to initialize integrator at initial time
             *
             * The method initializes the multistep scheme for an integration with step-size h starting at given initial time and condition 
             * @param[in] m number of saved steps
             * @param[in] ti initial time instant
             * @param[in] h step size
             * @param[in] x0 vector of initial states
             * @param[out] f vector of saved state vectors for multistep scheme
             * @return
             */     
            virtual  int initialize(const int &m, const double &ti, const double &h, const std::vector<T> &x0, std::vector<std::vector<T> > &f) const = 0;



        };

    }
}

#endif // SMARTMATH_BASE_RUNGEKUTTA_H
