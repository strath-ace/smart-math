/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef SMARTMATH_VANDERPOL_H
#define SMARTMATH_VANDERPOL_H

#include "base_dynamics.h"
#include "../exception.h"

namespace smartmath
{
    namespace dynamics {

        /**
         * @brief The Van der Pol oscillator
         *
         * The dynamics of the Van der Pol oscillator is defined as
         * \f{eqnarray*}{
         * \ddot{x} - \mu (1-x^2)\dot{x} + x = 0
         * \f}
         * where \f$\mu\f$ is user supplied parameter that can be either  aconstant or a polynomial value
         */
        template < class T >
        class vanderpol: public base_dynamics<T>
        {

        private:
            using base_dynamics<T>::m_name;

        public:
            /**
             * @brief vanderpol constructor
             * @param mu problem parameter (constant or polynomial value)
             */
            vanderpol(const T &mu=.5): base_dynamics<T>("Van der Pol dynamical system"), m_mu(mu){}

            /**
              * @brief ~vanderpol
              */
            ~vanderpol(){}

            /**
             * @brief evaluate evaluates the dynamics of the Van der Pol oscillator at a given instant of time and a given state.
             *
             * Function to evaluate the dynamics of the Van der Pol Oscillator at a given instant of time and a given state. It is a virtual function so any class that inherits from base_dynamics need to implement it.
             * @param[in] t time
             * @param[in] state state values at time t
             * @param[out] dstate derivative of the states at time t
             * @return
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const
		{
		    //sanity checks
		    if(t<0)
			smartmath_throw(m_name+": negative time supplied in evaluation of the dynamical system");
		    if(state.size()!=2)
			smartmath_throw(m_name+": the state dimension needs to be 2");

		    dstate.clear();

		    dstate.push_back(state[1]);
		    dstate.push_back(m_mu*(1.0-state[0]*state[0])*state[1]-state[0]);

		    return 0;
		}


        private:
            T m_mu;
        };

    }
}


#endif // SMARTMATH_VANDERPOL_H
