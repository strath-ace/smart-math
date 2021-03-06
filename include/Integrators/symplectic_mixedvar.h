/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2017 University of Strathclyde and Authors-------
-------- e-mail: romain.serra@strath.ac.uk ---------------------------
--------- Author: Romain Serra ---------------------------------------
*/

#ifndef SMARTMATH_SYMPLECTIC_MIXEDVAR_H
#define SMARTMATH_SYMPLECTIC_MIXEDVAR_H

#include "../Dynamics/hamiltonian_mixedvar.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The %symplectic_mixedvar class is a template abstract class. Any symplectic integrator using mixed variables needs to inherit from it and implement the method integrate()
         *
         * The %symplectic_mixedvar class is a template abstract class. Any symplectic integrator using mixed variables needs to inherit from it and implement the method that integrates between to given times, initial state and stepsize.
         * The use of mixed variables requires the dynamics to be formed in two parts: one accounting for the perturbed Hamiltonian and computed with the primary set (position-momentum) while the rest is evaluated with the secondary variables (action-angles). 
         * Usual symplectic schemes are then modified to take advantage of this formulation by performing part of the integration with the primary variables and the rest with the secondaries, requiring back and forth conversions between the two sets.
         */
        template < class T >
        class symplectic_mixedvar: public base_symplectic<T>
        {

        public:

            using base_symplectic<T>::integrate;

            /**
             * @brief symplectic_mixedvar constructor
             *
             * The constructor initializes the name of the integrator and a pointer to the dynamical system to be integrated
             * @param name integrator name
             * @param dyn pointer to a Hamiltonian dynamics with mixed variables
             * @param stages integer stating the number of integration stages in one step
             */
            symplectic_mixedvar(const std::string &name, const dynamics::hamiltonian_mixedvar<T> *dyn, const unsigned int &stages): base_symplectic<T>(name, NULL, stages), m_mix(dyn){}

            /**
             * @brief ~symplectic_mixedvar deconstructor
             */
            virtual ~symplectic_mixedvar(){}

            /**
             * @brief integration_step performs one integration step from the symplectic scheme with mixed variables
             *
             * The method implements one step of a symplectic scheme with mixed variables to integrate with given initial time,
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
                if(q0.size() != m_mix->get_dim())
                    smartmath_throw("INTEGRATION_STEP: position vector must have consistent dimension with Hamiltonian system");               
                if(p0.size() != m_mix->get_dim())
                    smartmath_throw("INTEGRATION_STEP: momentum vector must have consistent dimension with Hamiltonian system");     

                unsigned int n = m_mix->get_dim();
                std::vector<T> q = q0, p = p0, dq = q0, dp = p0;
                qf = q0;
                pf = p0;                

                /* performing the integration step per se using the precomputed coefficients */
                for(unsigned int j = 0; j < m_stages; j++)
                {

                    if(m_c[j] != 0.0)
                    { // drift with second set of coordinates
                        
                        m_mix->conversion(qf, pf, q, p);
                        qf = q;
                        pf = p;                        

                        m_mix->DHp2(ti, q, p, dp);
                        for(unsigned int i = 0; i < n; i++)
                            qf[i] += m_c[j] * tau * dp[i];

                        m_mix->conversion2(qf, pf, q, p);
                        qf = q;
                        pf = p;
                    }

                    if(m_d[j] != 0.0)
                    { // kick with first set of coordinates
                        
                        m_mix->DHq(ti, q, p, dq);
                        for(unsigned int i = 0; i < n; i++)
                            pf[i] -= m_d[j] * tau * dq[i];

                        p = pf;
                    }       

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
                if(x0.size() != 2 * m_mix->get_dim())
                    smartmath_throw("INTEGRATION: state vector must have consistent dimension with Hamiltonian system"); 

                t_history.clear();
                x_history.clear();

                double t = ti, h = (tend-ti) / double(nsteps);

                std::vector<T> x = x0;

                /* splitting the initial state vector */
                T zero = 0.0 * x0[0];
                unsigned int n = m_mix->get_dim();
                std::vector<T> q0(n, zero), p0(n, zero);
                for(unsigned int j = 0; j < n; j++)
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
                    for(unsigned int j = 0; j < n; j++)
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
            using base_symplectic<T>::m_name;
            using base_symplectic<T>::m_c;
            using base_symplectic<T>::m_d;
            using base_symplectic<T>::m_stages;
            /**
             * @brief m_mix pointer to Hamiltonian dynamics with mixed variables
             */            
            const dynamics::hamiltonian_mixedvar<T> *m_mix;

        };
    }
}

#endif // SMARTMATH_SYMPLECTIC_MIXEDVAR_H
