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
         * @brief The leapfrog_mixedvar class is an adaptation of the leapfrog algorithm for mixed variables.
         *
         * The leapfrog class is an adaptation of the leapfrog algorithm for mixed variables. It has two versions: kick-drift-kick and drift-kick-drift.
         */
        template < class T >
        class leapfrog_mixedvar: public symplectic_mixedvar<T>
        {

        protected:
            using symplectic_mixedvar<T>::m_ham;
            /**
             * @brief m_flag boolean to chose between drift-kick and kick-drift versions
             */             
            bool m_flag;

        public:

            /**
             * @brief leapfrog_mixedvar constructor
             *
             * The constructor initialize a pointer to the dynamics to integrate and a flag to decide on the integration scheme to use
             * @param dyn Hamiltonian system to integrate
             * @param flag boolean to know what algorithm to use (kick-drift-kick or drift-kick-drift)
             */
            leapfrog_mixedvar(const dynamics::hamiltonian_mixedvar<T> *dyn, const bool &flag) : symplectic_mixedvar<T>("leapfrog integrator with mixed variables", dyn), m_flag(flag){
                
                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("LEAPFROG_MIXEDVAR: symplectic integrator cannot operate on non-separable Hamiltonian");

            }

            /**
             * @brief ~leapfrog_mixedvar deconstructor
             */
            ~leapfrog_mixedvar(){}


            /**
             * @brief integration_step performs one integration step from the leapfrog for mixed variables
             *
             * The method implements one step of the leapfrog with mixed variables to integrate with given initial time,
             * final time, initial state condition(constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] tau time step
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integration_step(const double &ti, const double &tau, const std::vector<T> &x0, std::vector<T> &xfinal) const{

                /* sanity checks */
                if(x0.size() != 2 * m_ham->get_dim())
                    smartmath_throw("BASE_SYMPLECTIC: state vector must have consistent dimension with Hamiltonian system");
                if(xfinal.size() != x0.size())
                    smartmath_throw("BASE_SYMPLECTIC: initial and final states must have same dimension");      
                
                std::vector<T> q0, p0;
                int n = m_ham->get_dim();
                for(int i = 0; i < n; i++)
                {
                    q0.push_back(x0[i]);
                    p0.push_back(x0[i + n]);
                }
                std::vector<T> q = q0, p = p0, dq = q0, dp = p0;

                if(m_flag)
                { // drift-kick-drift
                    m_ham->conversion(q, p, q0, p0);
                    q = q0;
                    p = p0;

                    m_ham->DHp2(ti, q0, p0, dp);
                    for(int i = 0; i < n; i++)
                        q[i] += tau * dp[i] / 2.0;

                    m_ham->conversion2(q, p, q0, p0);
                    q = q0;
                    p = p0;
                    
                    m_ham->DHq(ti, q, p0, dq);
                    for(int i = 0; i < n; i++)
                        p[i] -= tau * dq[i];

                    m_ham->conversion(q, p, q0, p0);
                    q = q0;
                    p = p0; 

                    m_ham->DHp2(ti, q0, p0, dp);
                    for(int i = 0; i < n; i++)
                        q[i] += tau * dp[i] / 2.0;

                    m_ham->conversion2(q, p, q0, p0);
                    q = q0;
                    p = p0; 
                }
                else
                { // kick-drift-kick
                    m_ham->DHq(ti, q, p0, dq);
                    for(int i = 0; i < n; i++)
                        p[i] -= tau * dq[i] / 2.0;

                    m_ham->conversion(q, p, q0, p0);
                    q = q0;
                    p = p0;
                    
                    m_ham->DHp2(ti, q0, p0, dp);
                    for(int i = 0; i < n; i++)
                        q[i] += tau * dp[i];

                    m_ham->conversion2(q, p, q0, p0);
                    q = q0;
                    p = p0; 

                    m_ham->DHq(ti, q, p0, dq);
                    for(int i = 0; i < n; i++)
                        p[i] -= tau * dq[i] / 2.0;  
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

#endif // SMARTMATH_LEAPFROG_MIXEDVAR_H
