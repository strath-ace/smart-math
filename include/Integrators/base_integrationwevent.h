/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_BASE_INTEGRATIONWEVENT_H
#define SMARTMATH_BASE_INTEGRATIONWEVENT_H

#include "base_integrator.h"
#include "../exception.h"
#include "../Events/base_event.h"
#include <type_traits>

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The base_integrationwevent class is a template abstract class. Any integrator handling events added to the toolbox needs to inherit from it and implement the method integrate with events
         *
         * The base_integrationwevent class is a template abstract class. Any integrator handling events added to the toolbox needs to inherit from it and implement the method that integrates while dealing with events 
         */
        template < class T >
        class base_integrationwevent: public base_integrator<T>
        {

        protected:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            double m_minstep_events;
            double m_maxstep_events;
            std::vector<events::base_event<T>*> m_event_list;

        public:

            using base_integrator<T>::integrate;

            /**
             * @brief base_integrationwevent constructor
             *
             * The constructor specifically initializes a mimimum and maximum (0.0 corresponds to +infinity) stepzise for event detection
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             * @param m_minstep_events minimum step-size to detect an event
             * @param m_maxstep_events maximum step-size
             */
            base_integrationwevent(const std::string &name, const dynamics::base_dynamics<T> *dyn, const double &minstep_events, const double &maxstep_events) : base_integrator<T>(name, dyn), m_minstep_events(minstep_events), m_maxstep_events(maxstep_events){
                
                if(minstep_events<=0.0)
                   smartmath_throw("BASE_INTEGRATIONWEVENT: minimum step-size for events must be non negative");
                if(maxstep_events<0.0)
                   smartmath_throw("BASE_INTEGRATIONWEVENT: maximum step-size for events must be positive");

            }

            /**
             * @brief ~base_stepsizecontrol deconstructor
             */
            virtual ~base_integrationwevent(){}


            /**
             * @brief integrate method to integrate bewteen two given time steps, with initial condition and initial guess for step-size while handling events
             *
             * The method implements a variable step-size scheme to integrate with given initial time,
             * final time, initial state condition and initial guess for step-size
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps initial guess for number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of intermediate states
             * @param[out] t_history vector of intermediate times
             * @param[in] event function             
             * @return
             */
            virtual int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history, std::vector<int> (*g)(std::vector<T> x, double d)) const=0;

            /**
             * @brief integrate method to integrate bewteen two given time steps, with initial condition and initial guess for step-size
             *
             * The method implements a variable step-size scheme to integrate with given initial time,
             * final time, initial state condition and initial guess for step-size
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps initial guess for number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history) const{

                double t0=ti, tf=tend, n=nsteps;
                std::vector<T> x(x0);

                integrate(t0,tf,n,x,x_history,t_history,dummy_event);
    
                return 0;
            }


            /**
             * @brief integrate method to integrate bewteen two given time steps, with initial condition and initial guess for step-size while handling events
             *
             * The method implements a variable step-size scheme to integrate with given initial time,
             * final time, initial state condition and initial guess for step-size
             * @param[in] ti initial time instant
             * @param[out] tend final time instant
             * @param[in] nsteps initial guess for number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @param[in] event function
             * @return
             */
            int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal, std::vector<int> (*g)(std::vector<T> x, double d)) const{

                std::vector<std::vector<T> > x_history;
                std::vector<double> t_history;

                integrate(ti, tend, nsteps, x0, x_history, t_history, *g);

                xfinal=x_history.back();

                return 0;
            }

            /**
             * @brief returns a vector with one component equal to integer 0
             *
             * @param[in] x state vector
             * @param[in] d time
             * @return vector with one 0
             */
            static std::vector<int> dummy_event(std::vector<T> x, double d){

                std::vector<int> output(1,0);

                return output;
            }

            /**
             * @brief assings a list of events to the integrators if the events are handled via this approach (other way is using a function as dummy_event)
             *
             * @param[in] list of events
             */
            void set_event_list(std::vector<events::base_event<T>*> &event_list){

                m_event_list=event_list;

            }

            /**
             * @brief returns a double equal to the input for real numbers and something meaningful for polynomials
             *
             * @param[in] x estimated error
             * @return double equal to x for real numbers and something else for polynomials in smartuq
             */
            double evaluate_squarerootintegrationerror(const float &x) const{
                return x;
            }
            double evaluate_squarerootintegrationerror(const double &x) const{
                return x;
            }
            double evaluate_squarerootintegrationerror(const long double &x) const{
                return x;
            }    
	
        };

    }
}

#endif // SMARTMATH_BASE_INTEGRATIONWEVENT_H
