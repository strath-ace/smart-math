/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
-------------------- Author: Francesco Torre -------------------------
-------------- e-mail: francesco.torre@strath.ac.uk ------------------
*/


#ifndef SMARTMATH_SPRING_H
#define SMARTMATH_SPRING_H

#include "base_dynamics.h"
#include "../exception.h"

namespace smartmath
{
    namespace dynamics {

        /**
         * @brief The spring problem
         *
         * The spring problem models the dynamics of an object of mass \f$m\f$ linked to the equilibrium point by a spring.
         *
         */
        template < class T >
        class spring: public base_dynamics<T>
        {

        private:
            double m_k;
            double m_t_scale;
            double m_r_scale;

        public:

            /**
             * @brief spring constructor
             *
             * The constructor initializes the problem parameters and scaling factors to the value supplied by the user.
             * Default values for scaling factors are 1 and parameters is a vector of zero values or zero polynomials.
             * @param k elastic coefficient of the spring
             * @param t_scale time scaling factor
             * @param r_scale position scaling factor
             */
            spring(const double &k, const double &t_scale=1, const double &r_scale=1) : base_dynamics<double>("Spring Problem"),
                m_k(k), m_t_scale(t_scale), m_r_scale(r_scale)
            {
                if(m_k < 0.0)
                    smartmath_throw(this->m_name+": the elastic coefficient must be >= 0");

            }

            /**
              * @brief ~spring deconstructor
              */
            ~spring(){}

            /**
             * @brief evaluate evaluate the dinamics of the Two-body problem at a given instant of time and a given state.
             *
             * Function to evaluate the dinamics of the Two-body problem at a given instant of time and a given state. It is a virtual function so any class that inherites from base_dynamics need to implement it.
             * @param[in] t time
             * @param[in] state state values at time t
             * @param[out] dstate derivative of the states at time t
             * @return
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const{
                //sanity checks
                if(t<0)
                    smartmath_throw(this->m_name+": negative time supplied in evaluation of the dynamical system");
                if(state.size() > 6)
                    smartmath_throw(this->m_name+": it's difficult to imagine a spring acting ot time or higher dimensions...");
                if(state.size()%2 != 0)
                    smartmath_throw(this->m_name+": half dimensions are not contemplated in this version. Wait for update #42.");

                unsigned int dimension = state.size()/2;

                dstate.clear();

                // velocities
                for(unsigned int index = 0; index < dimension; ++index)
                {
                    dstate.push_back(state[dimension+index]);
                }

                // acceerations
                for(unsigned int index = 0; index < dimension; ++index)
                {
                    dstate.push_back(-m_k*state[index]);
                }

                if(dstate.size() != state.size())
                    smartmath_throw(this->m_name+": state and dstate dimensions mismatch.");

                return 0;
}
        };

    }
}

#endif // SMARTMATH_SPRING_H
