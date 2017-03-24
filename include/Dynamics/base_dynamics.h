/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef SMARTMATH_BASE_DYNAMICS_H
#define SMARTMATH_BASE_DYNAMICS_H

#include <vector>
#include <cmath>
#include "../exception.h"

#include "../LinearAlgebra/Eigen/Core"

namespace smartmath
{
    namespace dynamics {

        /**
         * @brief The base_dynamics class is a template abstract class. Any dynamics added to the toolbox needs to inherit from it and implement the method evaluate()
         *
         * The base_dynamics class is an abstract class. Any dynamical system added to the toolbox need to extend this class and implement the method evaluate.
         */
        template < class T >
        class base_dynamics
        {

        public:

            /**
             * @brief base_dynamics constructor.
             *
             * The base dynamics classed is constructed with just the name of the dynamical system. Being an abstract class it is not possible to instantiate an object of this class.
             * The constructor will be called by the sublcasses.
             * The dynamics are implemented as template classes because they can operate in the space pf real number or in the algebra of polynomials.
             * @param name dynamical system name
             */
            base_dynamics(const std::string &name): m_name(name){}

            virtual ~base_dynamics(){}

            /**
             * @brief evaluate evaluate the dinamics at a given instant of time and a given state.
             *
             * Function to evaluate the dinamics at a given instant of time and a given state. It is a virtual function so any class that inherites from base_dynamics need to implement it.
             * @param[in] t time
             * @param[in] state state values at time t
             * @param[out] dstate derivative of the states at time t
             * @return
             */
            virtual int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const = 0;

            virtual int evaluate(const double &t, const Eigen::VectorXd &state, Eigen::Ref<Eigen::VectorXd> dstate) const
            { smartmath_throw("evaluate_function using Eigen not implemented "); return 1; }

            /**
             * @brief get_name return dynamical system name
             *
             * Function to get the name of the dynamical system
             * @return
             */
            std::string get_name() const{return m_name;}

        protected:
            /**
             * @brief m_name dynamical system name (used in error messages)
             */
            std::string m_name;

        };

    }
}

#endif // SMARTMATH_BASE_DYNAMICS_H
