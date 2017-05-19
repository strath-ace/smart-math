/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_LEAPFROG_MIXEDVAR_H
#define SMARTMATH_LEAPFROG_MIXEDVAR_H

#include "symplectic_mixedvar.h"
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
        class leapfrog_mixedvar: public symplectic_mixedvar<T>
        {

        protected:
            using symplectic_mixedvar<T>::m_dyn;
            using symplectic_mixedvar<T>::m_order;
            using symplectic_mixedvar<T>::m_c;
            using symplectic_mixedvar<T>::m_d;
            bool m_flag;

        public:

            /**
             * @brief base_rungekutta constructor
             *
             * The constructor initialize the name of the integrator and a pointer to the dynamical system to be integrated
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             */
            leapfrog_mixedvar(const dynamics::hamiltonian_mixedvar<T> *dyn, const bool &flag) : symplectic_mixedvar<T>("leapfrog integrator with mixed variables", dyn, 2), m_flag(flag){
                
                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("LEAPFROG_MIXEDVAR: symplectic integrator cannot operate on non-separable Hamiltonian");

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
            ~leapfrog_mixedvar(){}


            /**
             * @brief integration_step performs one integration step from the Runge-Kutta scheme
             *
             * The method implements one step of a Runge-Kutta scheme to integrate with given initial time,
             * final time, initial state condition(constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] h time step
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integration_step(const double &ti, const double &tau, const std::vector<T> &x0, std::vector<T> &xfinal) const{

                std::vector<T> q0, p0;
                int n = m_dyn->get_dim();
                for(int i = 0; i < n; i++){
                    q0.push_back(x0[i]);
                    p0.push_back(x0[i + n]);
                }
                std::vector<T> q = q0, p = p0, dq = q0, dp = p0;

                m_dyn->DHp(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += m_c[0] * tau * dp[i];
                m_dyn->DHq(ti, q, p0, dq);
                for(int i = 0; i < n; i++)
                    p[i] -= m_d[0] * tau * dq[i];
                m_dyn->conversion(q, p, q0, p0);
                
                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += m_c[1] * tau * dp[i];
                m_dyn->DHq2(ti, q, p0, dq);
                for(int i = 0; i < n; i++)
                    p[i] -= m_d[1] * tau * dq[i];
                m_dyn->conversion2(q, p, q0, p0);

                xfinal.clear();
                for(int i = 0; i < n; i++)
                    xfinal.push_back(q[i]);
                for(int i = 0; i < n; i++)
                    xfinal.push_back(p[i]);
                               
                return 0;
            }

        };

    }
}

#endif // SMARTMATH_LEAPFROG_MIXEDVAR_H
