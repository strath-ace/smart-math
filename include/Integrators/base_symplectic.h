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
             * The constructor initializes the name of the integrator, a pointer to the dynamical system to be integrated and the number of stages of the integration
             * @param name integrator name
             * @param dyn pointer to a Hamiltonian dynamics
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
             * @param[in] q0 vector of initial coordinates
             * @param[in] p0 vector of initial momenta
             * @param[out] qf vector of final coordinates
             * @param[out] pf vector of final momenta
             * @return
             */
            int integration_step(const double &ti, const double &tau, const std::vector<T> &q0, const std::vector<T> &p0, std::vector<T> &qf, std::vector<T> &pf) const{

                /* sanity checks */
                if(q0.size() != m_ham->get_dim())
                    smartmath_throw("INTEGRATION_STEP: position vector must have consistent dimension with Hamiltonian system");               
                if(p0.size() != m_ham->get_dim())
                    smartmath_throw("INTEGRATION_STEP: momentum vector must have consistent dimension with Hamiltonian system");     

                int n = m_ham->get_dim();
                std::vector<T> q = q0, p = p0, dq = q0, dp = p0;
                qf = q0;
                pf = p0;

                for(int j = 0; j < m_stages; j++)
                {
                    if(m_c[j] != 0.0)
                    {
                        m_ham->DHp(ti, q, p, dp);
                        for(int i = 0; i < n; i++)
                            qf[i] += m_c[j] * tau * dp[i];
                    }

                    if(m_d[j] != 0.0)
                    {
                        m_ham->DHq(ti, qf, p, dq);    
                        for(int i = 0; i < n; i++)
                            pf[i] -= m_d[j] * tau * dq[i];
                    }

                    q = qf;
                    p = pf;
                }
                               
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
                    smartmath_throw("INTEGRATION: state vector must have consistent dimension with Hamiltonian system"); 

                t_history.clear();
                x_history.clear();

                double t = ti, h = (tend-ti) / double(nsteps);

                std::vector<T> x = x0;

                /* splitting the initial state vector */
                T zero = 0.0 * x0[0];
                std::vector<T> q0(3, zero), p0(3, zero);
                int n = m_ham->get_dim();
                for(int j = 0; j < n; j++)
                {
                    q0[j] = x0[j];
                    p0[j] = x0[j + n];
                }
                std::vector<T> q = q0, p = p0;

                for(int i = 0; i < nsteps; i++)
                {
                    integration_step(t, h, q0, p0, q, p);
                    t += h;
                    q0 = q;
                    p0 = p;
                    for(int j = 0; j < n; j++)
                    {
                        x[j] = q0[j];
                        x[j + n] = p0[j];
                    }                    
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
