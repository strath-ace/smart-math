/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2017 University of Strathclyde and Authors-------
-------- e-mail: romain.serra@strath.ac.uk ---------------------------
--------- Author: Romain Serra ---------------------------------------
*/

#ifndef SMARTMATH_BULIRSCHSTOER_H
#define SMARTMATH_BULIRSCHSTOER_H

#include "base_integrator.h"
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
        class bulirschstoer: public base_integrator<T>
        {

        protected:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            /**
             * @brief m_sequence Bulirsch sequence
             */              
            std::vector<unsigned int> m_sequence;
            /**
             * @brief m_extrapol size of the extrapolation table
             */               
            unsigned int m_extrapol;

        public:

            using base_integrator<T>::integrate;

            /**
             * @brief bulirschstoer constructor
             *
             * The constructor initializes the name of the integrator and a pointer to the dynamical system to be integrated
             * @param dyn pointer to the dynamical system to be integrated
             * @param extrapol size of the extrapolation table (the order equals twice that number)
             */
            bulirschstoer(const dynamics::base_dynamics<T> *dyn, const unsigned int &extrapol = 7) : base_integrator<T>("Bulirsch-Stoer algorithm", dyn), m_extrapol(extrapol){

                /* sanity checks */
                if(m_extrapol < 1)
                    smartmath_throw("BULIRSCHSTOER: number of extrapolations needs to be non negative"); 

                /* defining Bulirsch sequence */
                std::vector<unsigned int> sequence(m_extrapol);
                sequence[0] = 2;
                if(m_extrapol > 1)
                    sequence[1] = 4;
                if(m_extrapol > 2)
                    sequence[2] = 6;
                if(m_extrapol > 3)
                {
                    for(unsigned int k = 3; k < m_extrapol; k++)
                        sequence[k] = 2 * sequence[k - 2];
                }
                m_sequence = sequence;
            }

            /**
             * @brief ~bulirschstoer deconstructor
             */
            ~bulirschstoer(){}

            /**
             * @brief integration_step performs one integration step from the Bulirsch-Stoer method
             *
             * The method implements one step of the Adam Bashforth 6 scheme to integrate with given initial time,
             * final time, initial state condition(constant stepsize)
             * @param[in] t initial time for integration step 
             * @param[in] H stepsize
             * @param[in] x0 vector of initial states          
             * @param[out] xfinal vector of final states
             * @return
             */
            int integration_step(const double &t, const double &H, const std::vector<T> &x0, std::vector<T> &xfinal) const{

                std::vector< std::vector<T> > M;
                extrapolation(m_extrapol, H, x0, t, M);
                xfinal = M[m_extrapol-1];
                
                return 0;
            }

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

            /**
             * @brief midpoint computes the mid-point rule for all the micro-steps 
             *
             * The method computes the mid-point rule for all the micro-steps 
             * @param[in] n number of micro-steps
             * @param[in] H step-size                          
             * @param[in] y state vector at time t
             * @param[in] t time
             * @param[out] eta estimation of the state at t + H
             * @return
             */ 
            int midpoint(const unsigned int &n, const double &H, const std::vector<T> &y, const double &t, std::vector<T> &eta) const{

                /* sanity checks */
                if(n < 2)
                    smartmath_throw("MIDPOINT: number of micro-steps needs to be higher or equal to 2"); 
                if(floor(double(n) / 2.0) != double(n) / 2.0)
                    smartmath_throw("MIDPOINT: number of micro-steps has to be even");                 

                unsigned int s = y.size();
                std::vector<T> f = y, u0 = y, u1 = y;
                double h = H / double(n), h2 = 2.0 * h;

                m_dyn->evaluate(t, y, f);
                for(unsigned int j = 0; j < s; j++)
                    u1[j] += h * f[j];

                std::vector<T> u2 = u0;
                m_dyn->evaluate(t + h, u1, f);
                for(unsigned int j = 0; j < s; j++)
                    u2[j] += h2 * f[j];

                std::vector<T> v = y, w = y;
                for(unsigned int i = 2; i <= n; i++)
                {
                    v = u2;
                    u2 = u1;
                    m_dyn->evaluate(t + i * h, v, f);
                    for(unsigned int j = 0; j < s; j++)
                        u2[j] += h2 * f[j];
                    w = u1;
                    u1 = v;
                    u0 = w;
                }

                eta = y;
                for(unsigned int j = 0; j < s; j++)
                    eta[j] = u0[j] / 4.0 + u1[j] / 2.0 + u2[j] / 4.0;

                return 0;
            } 

            /**
             * @brief extrapolation computes backward differences  
             *
             * The method computes the backward differences for the Adam scheme
             * @param[in] i size of the desired extrapolation table
             * @param[in] H step-size                          
             * @param[in] y state vector at time t
             * @param[in] t time
             * @param[out] M last line of extrapolatio table (vector of state vectors)
             * @return
             */ 
            int extrapolation(const unsigned int &i, const double &H, const std::vector<T> &y, const double &t, std::vector<std::vector<T> > &M) const{

                /* sanity checks */
                if(i < 1)
                    smartmath_throw("EXTRAPOLATION: the extrapolation scheme needs to have at least one step"); 
                
                double aux1, aux2;
                std::vector<T> eta = y;
                unsigned int s = y.size();

                midpoint(m_sequence[i-1], H, y, t, eta);
                M.clear();
                M.push_back(eta);

                if(i > 1)
                {
                    std::vector<std::vector<T> > Mp;
                    extrapolation(i - 1, H, y, t, Mp); // recursive call

                    for(unsigned int j = 1; j < i; j++)
                    {
                        eta = M[j - 1];
                        M.push_back(eta);
                        aux1 = double(m_sequence[i-1]) / double(m_sequence[i-1-j]);
                        aux2 = aux1 * aux1 - 1.0;
                        for(unsigned int k = 0; k < s; k++)
                            M[j][k] += (eta[k] - Mp[j-1][k]) / aux2;
                    }
                }

                return 0;
            } 


        };

    }
}

#endif // SMARTMATH_BULIRSCHSTOER_H
