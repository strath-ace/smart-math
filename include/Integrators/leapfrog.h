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
            using base_symplectic<T>::m_name;
            using base_symplectic<T>::m_dyn;
            int m_order;
            std::vector<double> m_c, m_d;
            bool m_flag;

        public:

            using base_symplectic<T>::integrate;

            /**
             * @brief base_rungekutta constructor
             *
             * The constructor initialize the name of the integrator and a pointer to the dynamical system to be integrated
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             */
            leapfrog(const dynamics::base_hamiltonian<T> *dyn, const bool &flag) : base_symplectic<T>("leapfrog integrator", dyn), m_order(2), m_flag(flag){
                
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

                for(int j = 0; j < m_order; j++){
                    m_dyn->DHp(ti, q0, p0, dp);
                    for(int i = 0; i < n; i++)
                        q[i] += m_c[j] * tau * dp[i];
                    m_dyn->DHq(ti, q, p0, dq);
                    for(int i = 0; i < n; i++)
                        p[i] -= m_d[j] * tau * dq[i];
                    q0 = q;
                    p0 = p;
                }

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

#endif // SMARTMATH_LEAPFROG_H
