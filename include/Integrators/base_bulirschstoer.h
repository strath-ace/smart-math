/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_BASE_BULIRSCHSTOER_H
#define SMARTMATH_BASE_BULIRSCHSTOER_H

#include "base_integrator.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The base_multistep class is a template abstract class. Any fixed-size, fixed order multistep integrator added to the toolbox needs to inherit from it and implement the method integration_step() as well as initialize()
         *
         * The base_multistep class is a template abstract class. Any fixed-size, fixed order multistep integrator added to the toolbox needs to inherit from it and implement the methods integration_step() and initialize() 
         */
        template < class T >
        class base_bulirschstoer: public base_integrator<T>
        {

        protected:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            std::vector<int> m_orders;
            int m_stage;

        public:

            using base_integrator<T>::integrate;

            /**
             * @brief base_multistep constructor
             *
             * The constructor initialize the name of the integrator and a pointer to the dynamical system to be integrated
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             * @param order number of saved steps
             */
            base_bulirschstoer(const dynamics::base_dynamics<T> *dyn, const int &stage = 7) : base_integrator<T>("Bulirsch-Stoer algorithm", dyn), m_stage(stage){

                /* sanity checks */
                if(m_stage < 3)
                    smartmath_throw("base_bulirschstoer: "); 

                std::vector<int> orders(m_stage);
                orders[0] = 2;
                orders[1] = 4;
                orders[2] = 6;
                for(int k = 3; k < m_stage; k++)
                    orders[k] = 2 * orders[k - 2];
                m_orders = orders;

            }

            /**
             * @brief ~base_multistep deconstructor
             */
            ~base_bulirschstoer(){}

            /**
             * @brief integration_step performs one integration step from the Adam Bashforth (order 6) method
             *
             * The method implements one step of the Adam Bashforth 6 scheme to integrate with given initial time,
             * final time, initial state condition(constant stepsize)
             * @param[in] t initial time for integration step 
             * @param[in] m order
             * @param[in] h step-size
             * @param[in] x0 vector of initial states
             * @param[in] f vector of saved state vectors for multistep scheme             
             * @param[out] xfinal vector of final states
             * @return
             */
            int integration_step(const double &t, const double &H, const std::vector<T> &x0, std::vector<T> &xfinal) const{

                std::vector< std::vector<T> > M;
                extrapolation(m_stage, H, x0, t, M);
                xfinal = M[m_stage-1];
                
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

                std::vector<T> x(x0),xp(x0);
                std::vector<std::vector<T> > f;
                double t=ti, H = (tend-ti)/nsteps;

                for(int k=0; k<nsteps; k++){

                    integration_step(t,H,x,xp);

                    /* Saving states */
                    t+=H;
                    x=xp;
                    t_history.push_back(t);
                    x_history.push_back(x);

                }

                return 0;
            }

            /**
             * @brief backward_differences computes backward differences  
             *
             * The method computes the backward differences for the Adam scheme
             * @param[in] f vector of saved state vectors for multistep scheme
             * @param[in] m number of saved steps
             * @param[out] Df vector of backward differences
             * @return
             */ 
            int midpoint(const std::vector<T> &y, const int &n, const double &H, const double &t, std::vector<T> &eta) const{

                double h = H / double(n);
                eta = y;

                std::vector<T> f = y, u0 = y, u1 = y;

                m_dyn->evaluate(t, y, f);
                for(unsigned int j = 0; j < y.size(); j++)
                    u1[j] += h * f[j];

                std::vector<T> u2 = u0;
                m_dyn->evaluate(t + h, u1, f);
                for(unsigned int j = 0; j < y.size(); j++)
                    u2[j] += 2.0 * h * f[j];

                std::vector<T> v = y, w = y;
                for(int i = 2; i < n - 1; i++)
                {
                    v = u2;
                    m_dyn->evaluate(t + i * h, u2, f);
                    for(unsigned int j = 0; j < y.size(); j++)
                        u2[j] = u1[j] + 2.0 * h * f[j];
                    w = u1;
                    u1 = v;
                    u0 = w;
                }

                for(unsigned int j = 0; j < y.size(); j++)
                    eta[j] = u0[j] / 4.0 + u1[j] / 2.0 + u2[j] / 4.0;

                return 0;
            } 

            /**
             * @brief backward_differences computes backward differences  
             *
             * The method computes the backward differences for the Adam scheme
             * @param[in] f vector of saved state vectors for multistep scheme
             * @param[in] m number of saved steps
             * @param[out] Df vector of backward differences
             * @return
             */ 
            int extrapolation(const int &i, const double &H, const std::vector<T> &y, const double &t, std::vector<std::vector<T> > &M) const{

                M.clear();
                double aux;
                std::vector<T> eta = y;

                midpoint(y, m_orders[i-1], H, t, eta);
                M.push_back(eta);

                if(i > 1)
                {
                    std::vector<std::vector<T> > Mp;
                    extrapolation(i - 1, H, y, t, Mp); // recursive call

                    for(int j = 1; j < i; j++){
                        M.push_back(M[j - 1]);
                        aux = double(m_orders[i-1]) / double(m_orders[i-1-j]);
                        for(unsigned int k = 0; k < y.size(); k++)
                            M[j][k] += (M[j-1][k] - Mp[j-1][k]) / (aux * aux - 1.0);
                    }
                }

                return 0;
            } 


        };

    }
}

#endif // SMARTMATH_BASE_BULIRSCHSTOER_H
