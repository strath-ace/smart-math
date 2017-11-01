/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2017 University of Strathclyde and Authors-------
-------- e-mail: romain.serra@strath.ac.uk ---------------------------
--------- Author: Romain Serra ---------------------------------------
*/

#ifndef SMARTMATH_BASE_EMBEDDEDRK_H
#define SMARTMATH_BASE_EMBEDDEDRK_H

#include "base_integrator.h"
#include "../exception.h"
#include <type_traits>

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The %base_embeddedRK class is a template abstract class. Any variable step-size Runge-Kutta algorithm added to the toolbox needs to inherit from it and implement the method integration_step()
         *
         * The %base_embeddedRK class is a template abstract class. Any variable step-size Runge-Kutta algorithm added to the toolbox needs to inherit from it and implement the method that performs on integration step between to given times given the initial state 
         */
        template < class T >
        class base_embeddedRK: public base_integrationwevent<T>
        {

        protected:
            using base_integrationwevent<T>::m_name;
            using base_integrationwevent<T>::m_dyn;
            using base_integrationwevent<T>::m_minstep_events;
            using base_integrationwevent<T>::m_maxstep_events;               
            /**
             * @brief m_tol tolerance for step-size control
             */
            double m_tol;
            /**
             * @brief m_multiplier maximum multiplying factor used when increasing step-size
             */            
            double m_multiplier;
            /**
             * @brief m_control order of scheme used for state estimation in step-size control
             */            
            unsigned int m_control;
            

        public:

            using base_integrationwevent<T>::integrate;

            /**
             * @brief base_embeddedRK constructor
             *
             * The constructor specifically initializes a tolerance for integration error and a maximum multiplier to increase the stepsize
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             * @param tol threshold used for acceptable estimated error
             * @param multiplier factor used to increase step-sized when judged necessary
             * @param minstep_events minimum step-size to detect an event
             * @param maxstep_events maximum step-size
             */
            base_embeddedRK(const std::string &name, const dynamics::base_dynamics<T> *dyn, const double &tol, const double &multiplier, const double &minstep_events, const double &maxstep_events) : base_integrationwevent<T>(name, dyn, minstep_events, maxstep_events), m_tol(tol), m_multiplier(multiplier){
                
                /** sanity checks **/
                if(tol <= 0.0)
                   smartmath_throw("BASE_embeddedRK: tolerance for estimated error must be non negative");
                if((multiplier > 5.0) || (multiplier < 2.0))
                   smartmath_throw("BASE_embeddedRK: maximum step-multiplier must be between 2 and 5");

            }

            /**
             * @brief ~base_embeddedRK deconstructor
             */
            virtual ~base_embeddedRK(){}

            /**
             * @brief integration_step performs one integration step from the integration scheme
             *
             * The method implements one step of a variable step-size algorithm to integrate with given initial time,
             * final time, initial state condition(constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] m method order
             * @param[in] h time step
             * @param[in] x0 vector of initial states
             * @param[in] f vector of saved state vectors (for multistep scheme only) 
             * @param[out] xfinal vector of final states
             * @param[out] er estimated error
             * @return
             */
            virtual int integration_step(const double &ti, const unsigned int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal, T &er) const = 0;

            /**
             * @brief integrate method to integrate bewteen two given time steps, with initial condition and initial guess for step-size while handling events
             *
             * The method implements a variable step-size scheme to integrate with given initial time,
             * final time, initial state condition and initial guess for step-size
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps initial guess for number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] x_history vector of intermediate states
             * @param[out] t_history vector of intermediate times
             * @param[in] g event function             
             * @return
             */
            int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history, std::vector<int> (*g)(std::vector<T> x, double d)) const{

                x_history.clear();
                t_history.clear();

                unsigned int k;
                int check = 0;
                std::vector<T> x(x0), xtemp(x0);
                std::vector<std::vector<T> > f;

                double factor = 1.0, value = 0.0, t = ti, h = (tend - ti) / double(nsteps);
                T er = 0.0 * x0[0];

                std::vector<int> events, events2;
                events = g(x0, ti);
                events2 = events;
                unsigned int m = events.size();           

                unsigned int i = 0;
                while(sqrt(pow(t - ti, 2)) < sqrt(pow(tend - ti, 2)))
                {

                    if((h * h > m_maxstep_events*m_maxstep_events) && (m_maxstep_events > 0.0))
                    {
                        if(h > 0.0)
                            h = m_maxstep_events;
                        else
                            h = -m_maxstep_events;
                    }   

                    if(sqrt(pow(tend - t, 2)) < sqrt(h * h))
                        h = tend - t;

                    integration_step(t, m_control, h, x, f, xtemp, er);
                    
                    /* Step-size control */
                    value = evaluate_squarerootintegrationerror(er);
                    factor=pow(m_tol / value, 1.0 / (double(m_control) + 1.0));
                    if(er > m_tol) // unsucessful step
                        h *= 0.9 * factor;  
                    else
                    { // sucessful step
                        /* Checking for the events */
                        events2 = g(xtemp, t + h);

                        k = 0;
                        check = 0;        
                        while(k < m)
                        {
                            if(events2[k] - events[k] != 0)
                            {
                                check = 1;
                                k = m;
                            }
                            else
                                k++;
                        }

                        if(check == 1)
                        {
                            if(sqrt(h * h) > m_minstep_events)
                                h *= 0.5;
                            else
                            {
                                tend = t + h; // saving the termination time    
                                t = tend; // trick to get out of the while loop   
                                x_history.push_back(xtemp);
                                t_history.push_back(t);

                                if(this->m_comments)
                                    std::cout << "Propagation interrupted by terminal event at time " << tend << " after " << i << " steps" << std::endl;
                            }
                        }
                        else
                        {
                            x = xtemp; // updating state
                            t += h; // updating current time  
                            events = events2; 
                            x_history.push_back(x);
                            t_history.push_back(t); 
                            /* Step-size control */
                            if(factor > m_multiplier)
                                factor = m_multiplier;  
                            h *= factor; // updating step-size
                            i++; // counting the number of steps
                        }           
        
                    }
                }

                return 0;
            }

        };

    }
}

#endif // SMARTMATH_BASE_EMBEDDEDRK_H
