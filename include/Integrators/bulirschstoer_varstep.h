/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_BULIRSCHSTOER_VARSTEP_H
#define SMARTMATH_BULIRSCHSTOER_VARSTEP_H

#include "base_integrationwevent.h"
#include "bulirschstoer.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The bulirschstoer class is an implementation of the Burlish-Stoer method with polynomial extrapolation
         *
         * The bulirschstoer class is a fixed-stepsize implementation of the Burlish-Stoer method with polynomial extrapolation and the original Burlish sequence
         */
        template < class T >
        class bulirschstoer_varstep: public bulirschstoer<T>, public base_integrationwevent<T>
        {

        protected:
            using base_integrationwevent<T>::m_name;
            using base_integrationwevent<T>::m_dyn;
            using bulirschstoer<T>::m_sequence;
            using bulirschstoer<T>::m_extrapol;

        public:

            using base_integrationwevent<T>::integrate;

            /**
             * @brief bulirschstoer constructor
             *
             * The constructor initializes the name of the integrator and a pointer to the dynamical system to be integrated
             * @param dyn pointer to the dynamical system to be integrated
             * @param extrapol size of the extrapolation table (the order equals twice that number)
             */
            bulirschstoer_varstep(const dynamics::base_dynamics<T> *dyn, const unsigned int &extrapol = 7) : base_integrationwevent<T>("Bulirsch-Stoer method with variable stepsize", dyn), m_extrapol(extrapol){

                /* sanity checks */
                if(m_extrapol < 1)
                    smartmath_throw("bulirschstoer: number of extrapolations needs to be non negative"); 

                /* defining Bulirsch sequence */
                std::vector<unsigned int> orders(m_extrapol);
                orders[0] = 2;
                if(m_extrapol > 1)
                    orders[1] = 4;
                if(m_extrapol > 2)
                    orders[2] = 6;
                if(m_extrapol > 3)
                {
                    for(int k = 3; k < m_extrapol; k++)
                        orders[k] = 2 * orders[k - 2];
                }
                m_sequence = orders;
            }

            /**
             * @brief ~bulirschstoer deconstructor
             */
            ~bulirschstoer_varstep(){}

            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method implements the multistep scheme to integrate with given initial time,
             * final time, initial state condition and number of steps (constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] x_history vector of intermediate state vector (including final one)
             * @param[out] t_history vector of intermediate times (including final one)
             * @return
             */
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history) const{
    
                t_history.clear();
                x_history.clear();

                std::vector<T> x(x0), xp(x0);
                double t = ti, H = (tend - ti) / double(nsteps);

                for(int k = 0; k < nsteps; k++){

                    integration_step(t, H, x, xp);

                    /* Saving states */
                    t += H;
                    x = xp;
                    t_history.push_back(t);
                    x_history.push_back(x);

                }

                return 0;
            }

        };

    }
}

#endif // SMARTMATH_BULIRSCHSTOER_VARSTEP_H
