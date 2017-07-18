/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_BASE_MULTISTEP_H
#define SMARTMATH_BASE_MULTISTEP_H

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
        class base_multistep: public base_integrator<T>
        {

        protected:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            /**
             * @brief m_order order of the multistep integrator
             */
            int m_order;

        public:

            using base_integrator<T>::integrate;

            /**
             * @brief base_multistep constructor
             *
             * The constructor initializes the name of the integrator and a pointer to the dynamical system to be integrated
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             * @param order number of saved steps
             */
            base_multistep(const std::string &name, const dynamics::base_dynamics<T> *dyn, const int &order) : base_integrator<T>(name, dyn), m_order(order){}

            /**
             * @brief ~base_multistep deconstructor
             */
            virtual ~base_multistep(){}

            /**
             * @brief integration_step performs one integration step 
             *
             * The method implements one step of a multi-step algorithm to integrate with given initial time,
             * final time, initial state condition(constant stepsize) returning the full history of propagation
             * @param[in] t initial time for integration step 
             * @param[in] m order
             * @param[in] h step-size
             * @param[in] x0 vector of initial states
             * @param[in] f vector of saved state vectors for multistep scheme             
             * @param[out] xfinal vector of final states
             * @return
             */
            virtual int integration_step(const double &t, const int &m, const double &h, const std::vector<T> &x0, std::vector<std::vector<T> > &f, std::vector<T> &xfinal) const=0;

            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method implements a multistep scheme to integrate with given initial time,
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
                double t=ti, h = (tend-ti)/nsteps;

                initialize(m_order,ti,h,x0,f);

                for(int k=0; k<nsteps; k++)
                {
                    integration_step(t,m_order,h,x,f,xp);

                    /* Saving states */
                    t+=h;
                    x=xp;
                    t_history.push_back(t);
                    x_history.push_back(x);
                }

                return 0;
            }

            /**
             * @brief initialize method to initialize integrator at initial time
             *
             * The method initializes the multistep scheme for an integration with step-size h starting at given initial time and condition 
             * @param[in] m number of saved steps
             * @param[in] ti initial time instant
             * @param[in] h step size
             * @param[in] x0 vector of initial states
             * @param[out] f vector of saved state vectors for multistep scheme
             * @return
             */     
            virtual  int initialize(const int &m, const double &ti, const double &h, const std::vector<T> &x0, std::vector<std::vector<T> > &f) const = 0;

        };

    }
}

#endif // SMARTMATH_BASE_MULTISTEP_H
