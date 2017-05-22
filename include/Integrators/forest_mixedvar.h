/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_FOREST_MIXEDVAR_H
#define SMARTMATH_FOREST_MIXEDVAR_H

#include "symplectic_mixedvar.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The forest_mixedvar class is an adaptation for mixed variables of the Forest integrator.
         *
         * The forest_mixedvar class is an adaptation for mixed variables by Yoshida of the Forest integrator (4th order).
         */
        template < class T >
        class forest_mixedvar: public symplectic_mixedvar<T>
        {

        protected:
            using symplectic_mixedvar<T>::m_dyn;
            double m_beta;

        public:

            /**
             * @brief forest_mixedvar constructor
             *
             * The constructor initialize a pointer to the dynamics to integrate
             * @param dyn Hamiltonian system to integrate
             */
            forest_mixedvar(const dynamics::hamiltonian_mixedvar<T> *dyn) : symplectic_mixedvar<T>("Forest integrator with mixed variables", dyn){
                
                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("FOREST_MIXEDVAR: symplectic integrator cannot operate on non-separable Hamiltonian");

                double beta = pow(2.0, 1.0 / 3.0);
                m_beta = beta;

            }

            /**
             * @brief ~forest_mixedvar deconstructor
             */
            ~forest_mixedvar(){}


            /**
             * @brief integration_step performs one integration step from the Forest scheme with mixed variables
             *
             * The method implements one step of the Forest algorithm with mixed variables to integrate with given initial time,
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

                m_dyn->conversion(q, p, q0, p0);
                q = q0;
                p = p0;                
                
                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += tau * dp[i] * 0.5 / (2.0 - m_beta);

                m_dyn->conversion2(q, p, q0, p0);
                q = q0;
                p = p0; 

                m_dyn->DHq(ti, q, p0, dq);
                for(int i = 0; i < n; i++)
                    p[i] -= tau * dq[i] * 1.0 / (2.0 - m_beta); 

                m_dyn->conversion(q, p, q0, p0);
                q = q0;
                p = p0; 

                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += tau * dp[i] * 0.5 * (1.0 - m_beta) / (2.0 - m_beta);

                m_dyn->conversion2(q, p, q0, p0);
                q = q0;
                p = p0; 

                m_dyn->DHq(ti, q, p0, dq);
                for(int i = 0; i < n; i++)
                    p[i] -= tau * dq[i] * (- m_beta / (2.0 - m_beta)); 

                m_dyn->conversion(q, p, q0, p0);
                q = q0;
                p = p0; 

                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += tau * dp[i] * 0.5 * (1.0 - m_beta) / (2.0 - m_beta);

                m_dyn->conversion2(q, p, q0, p0);
                q = q0;
                p = p0; 

                m_dyn->DHq(ti, q, p0, dq);
                for(int i = 0; i < n; i++)
                    p[i] -= tau * dq[i] * 1.0 / (2.0 - m_beta); 

                m_dyn->conversion(q, p, q0, p0);
                q = q0;
                p = p0; 

                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += tau * dp[i] * 0.5 / (2.0 - m_beta);

                m_dyn->conversion2(q, p, q0, p0);
                q = q0;
                p = p0;                 

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

#endif // SMARTMATH_FOREST_MIXEDVAR_H
