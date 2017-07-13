/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: romain.serra@strath.ac.uk -----------------------
----------------- Author: Romain Serra -------------------------------
*/

#ifndef SMARTMATH_BASE_SYMPLECTIC_H
#define SMARTMATH_BASE_SYMPLECTIC_H

#include "../Dynamics/base_hamiltonian.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The base_symplectic class is a template abstract class. Any sympletic integrator added to the toolbox needs to inherit from it
         *
         * The base_symplectic class is a template abstract class. Any symplectic integrator added to the toolbox needs to inherit from it
         */
        template < class T >
        class base_symplectic: public base_integrator<T>
        {

        public:

            using base_integrator<T>::integrate;

            /**
             * @brief base_integrator constructor
             *
             * The constructor initialize the name of the integrator, a pointer to the dynamical system to be integrated and the number of stages of the integration
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             * @param stages integer stating the number of integration stages in one step
             */
            base_symplectic(const std::string &name, const dynamics::base_hamiltonian<T> *dyn, const int &stages): base_integrator<T>(name, NULL), m_ham(dyn), m_stages(stages){}

            /**
             * @brief ~base_symplectic deconstructor
             */
            virtual ~base_symplectic(){}

            /**
             * @brief integration_step performs one integration step from the symplectic integrator
             *
             * The method implements one step of a symplectic scheme to integrate with given initial time,
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
                    smartmath_throw("INTEGRATION_STEP: state vector must have consistent dimension with Hamiltonian system");               

                /* reconstituting the initial canonical variables from the state vector */
                std::vector<T> q0, p0;
                int n = m_ham->get_dim();
                for(int i = 0; i < n; i++)
                {
                    q0.push_back(x0[i]);
                    p0.push_back(x0[i + n]);
                }
                std::vector<T> q = q0, p = p0, dq = q0, dp = p0;

                /* performing the integration step per se using the precomputed coefficients */
                for(int j = 0; j < m_stages; j++)
                {
                    if(m_c[j] != 0.0)
                    {
                        m_ham->DHp(ti, q0, p0, dp);
                        for(int i = 0; i < n; i++)
                            q[i] += m_c[j] * tau * dp[i];
                    }

                    if(m_d[j] != 0.0)
                    {
                        m_ham->DHq(ti, q, p0, dq);    
                        for(int i = 0; i < n; i++)
                            p[i] -= m_d[j] * tau * dq[i];
                    }

                    q0 = q;
                    p0 = p;
                }

                /* reconstituting the final state vector from the propagated canonical variables */
                xfinal.clear();
                for(int i = 0; i < n; i++)
                    xfinal.push_back(q[i]);
                for(int i = 0; i < n; i++)
                    xfinal.push_back(p[i]);
                               
                return 0;
            }
            
            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method implements a fixed-step symplectic scheme to integrate with given initial time,
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

                /* sanity checks */
                if(x0.size() != 2 * m_ham->get_dim())
                    smartmath_throw("INTEGRATE: state vector must have consistent dimension with Hamiltonian system"); 

                t_history.clear();
                x_history.clear();

                std::vector<T> dx = x0, x = x0, x_temp = x0;

                double t = ti, h = (tend-ti) / nsteps;

                for(int i = 0; i < nsteps; i++)
                {
                    integration_step(t, h, x, x_temp);
                    t += h;
                    x = x_temp;
                    t_history.push_back(t);
                    x_history.push_back(x);
                }

                return 0;
            }

        protected:
            using base_integrator<T>::m_name;
            /**
             * @brief m_ham pointer to Hamiltonian dynamics
             */
            const dynamics::base_hamiltonian<T> *m_ham;
            /**
             * @brief m_stages number of stages in integration step
             */
            int m_stages;
            /**
             * @brief m_c coefficients for drifts
             */
            std::vector<double> m_c;
            /**
             * @brief m_s coefficients for kicks
             */
            std::vector<double> m_d;

        };
    }
}

#endif // SMARTMATH_BASE_SYMPLECTIC_H
