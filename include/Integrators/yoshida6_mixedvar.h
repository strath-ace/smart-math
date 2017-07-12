/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_YOSHIDA6_MIXEDVAR_H
#define SMARTMATH_YOSHIDA6_MIXEDVAR_H

#include "symplectic_mixedvar.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The yoshida6_mixedvar class is an adaptation for mixed variables of the Yoshida 6 integrator.
         *
         * The yoshida6_mixedvar class is an adaptation for mixed variables of a 6th order integrator by Yoshida.
         */
        template < class T >
        class yoshida6_mixedvar: public symplectic_mixedvar<T>
        {

        protected:
            using symplectic_mixedvar<T>::m_dyn;
            /**
             * @brief m_c1 first coefficient for drift
             */            
            double m_c1;
            /**
             * @brief m_c2 second coefficient for drift
             */            
            double m_c2;            
            /**
             * @brief m_c3 third coefficient for drift
             */            
            double m_c3;
            /**
             * @brief m_c4 fourth coefficient for drift
             */            
            double m_c4;            
            /**
             * @brief m_d1 first coefficient for kicks
             */                  
            double m_d1;
            /**
             * @brief m_d2 second coefficient for kicks
             */                  
            double m_d2;
            /**
             * @brief m_d3 third coefficient for kicks
             */                  
            double m_d3;
            /**
             * @brief m_d4 fourth coefficient for kicks
             */                  
            double m_d4;

        public:

            /**
             * @brief yoshida6_mixedvar constructor
             *
             * The constructor initializes a pointer to the dynamics to integrate and precomputes a parameter necessary for integration
             * @param dyn Hamiltonian system to integrate
             */
            yoshida6_mixedvar(const dynamics::hamiltonian_mixedvar<T> *dyn) : symplectic_mixedvar<T>("Yoshida 6 integrator with mixed variables", dyn){
                
                /* sanity checks */
                if(dyn->is_separable() == false)
                    smartmath_throw("YOSHIDA6_MIXEDVAR: symplectic integrator cannot operate on non-separable Hamiltonian");

                m_c1 = 0.392256805238780;
                m_c2 = 0.510043411918458;
                m_c3 = -0.471053385409758;
                m_c4 = 0.687531682525198e-1;

                m_d1 = 0.784513610477560;
                m_d2 = 0.235573213359357;
                m_d3 = -0.117767998417887e1;
                m_d4 = 0.131518632068391e1;

            }

            /**
             * @brief ~yoshida6_mixedvar deconstructor
             */
            ~yoshida6_mixedvar(){}


            /**
             * @brief integration_step performs one integration step from the Forest scheme with mixed variables
             *
             * The method implements one step of the Forest algorithm with mixed variables to integrate with given initial time,
             * final time, initial state condition(constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] tau time step
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integration_step(const double &ti, const double &tau, const std::vector<T> &x0, std::vector<T> &xfinal) const{

                /* sanity checks */
                if(x0.size() != 2 * m_dyn->get_dim())
                    smartmath_throw("INTEGRATION_STEP: state vector must have consistent dimension with Hamiltonian system");
                if(xfinal.size() != x0.size())
                    smartmath_throw("INTEGRATION_STEP: initial and final states must have same dimension");      

                std::vector<T> q0, p0;
                int n = m_dyn->get_dim();
                for(int i = 0; i < n; i++)
                {
                    q0.push_back(x0[i]);
                    p0.push_back(x0[i + n]);
                }
                std::vector<T> q = q0, p = p0, dq = q0, dp = p0;

                m_dyn->conversion(q, p, q0, p0);
                q = q0;
                p = p0;                
                
                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += tau * dp[i] * m_c1;

                m_dyn->conversion2(q, p, q0, p0);
                q = q0;
                p = p0; 

                m_dyn->DHq(ti, q, p0, dq);
                for(int i = 0; i < n; i++)
                    p[i] -= tau * dq[i] * m_d1; 

                m_dyn->conversion(q, p, q0, p0);
                q = q0;
                p = p0; 

                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += tau * dp[i] * m_c2;

                m_dyn->conversion2(q, p, q0, p0);
                q = q0;
                p = p0; 

                m_dyn->DHq(ti, q, p0, dq);
                for(int i = 0; i < n; i++)
                    p[i] -= tau * dq[i] * m_d2; 

                m_dyn->conversion(q, p, q0, p0);
                q = q0;
                p = p0; 

                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += tau * dp[i] * m_c3;

                m_dyn->conversion2(q, p, q0, p0);
                q = q0;
                p = p0; 

                m_dyn->DHq(ti, q, p0, dq);
                for(int i = 0; i < n; i++)
                    p[i] -= tau * dq[i] * m_d3; 

                m_dyn->conversion(q, p, q0, p0);
                q = q0;
                p = p0; 

                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += tau * dp[i] * m_c4;

                m_dyn->conversion2(q, p, q0, p0);
                q = q0;
                p = p0;   

                m_dyn->DHq(ti, q, p0, dq);
                for(int i = 0; i < n; i++)
                    p[i] -= tau * dq[i] * m_d4; 

                m_dyn->conversion(q, p, q0, p0);
                q = q0;
                p = p0;     

                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += tau * dp[i] * m_c4; // c5 = c4

                m_dyn->conversion2(q, p, q0, p0);
                q = q0;
                p = p0;   

                m_dyn->DHq(ti, q, p0, dq);
                for(int i = 0; i < n; i++)
                    p[i] -= tau * dq[i] * m_d3; // d5 = d3

                m_dyn->conversion(q, p, q0, p0);
                q = q0;
                p = p0;                                                

                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += tau * dp[i] * m_c3; // c6 = c3

                m_dyn->conversion2(q, p, q0, p0);
                q = q0;
                p = p0;   

                m_dyn->DHq(ti, q, p0, dq);
                for(int i = 0; i < n; i++)
                    p[i] -= tau * dq[i] * m_d2; // d6 = d2

                m_dyn->conversion(q, p, q0, p0);
                q = q0;
                p = p0;   

                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += tau * dp[i] * m_c2; // c7 = c2

                m_dyn->conversion2(q, p, q0, p0);
                q = q0;
                p = p0;   

                m_dyn->DHq(ti, q, p0, dq);
                for(int i = 0; i < n; i++)
                    p[i] -= tau * dq[i] * m_d1; // d7 = d1

                m_dyn->conversion(q, p, q0, p0);
                q = q0;
                p = p0;   

                m_dyn->DHp2(ti, q0, p0, dp);
                for(int i = 0; i < n; i++)
                    q[i] += tau * dp[i] * m_c1; // c8 = c1

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

#endif // SMARTMATH_YOSHIDA6_MIXEDVAR_H
