/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_BASE_STEPSIZECONTROL_H
#define SMARTMATH_BASE_STEPSIZECONTROL_H

#include "base_integrator.h"
#include "../exception.h"
#include <type_traits>

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The base_stepsizecontrol class is a template abstract class. Any variable step-size algorithm added to the toolbox needs to inherit from it and implement the method integration_step()
         *
         * The base_stepsizecontrol class is a template abstract class. Any variable step-size algorithm added to the toolbox needs to inherit from it and implement the method that performs on integration step between to given times given the initial state 
         */
        template < class T >
        class base_stepsizecontrol: public base_integrator<T>
        {

        protected:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            double m_tol;
            double m_multiplier;
            double m_control;
            double m_minstep_events;
            double m_maxstep_events;

        public:

            using base_integrator<T>::integrate;

            /**
             * @brief base_stepsizecontrol constructor
             *
             * The constructor initialize the name of the integrator and a pointer to the dynamical system to be integrated
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             * @param tol threshold used for acceptable estimated error
             * @param multiplier factor used to increase step-sized when judged necessary
             * @param m_minstep_events minimum step-size to detect an event
             * @param m_maxstep_events maximum step-size
             */
            base_stepsizecontrol(const std::string &name, const dynamics::base_dynamics<T> *dyn, const double &tol, const double &multiplier, const double &minstep_events, const double &maxstep_events) : base_integrator<T>(name, dyn), m_tol(tol), m_multiplier(multiplier), m_minstep_events(minstep_events), m_maxstep_events(maxstep_events){
                
                /** sanity checks **/
                if(tol<=0.0)
                   smartmath_throw("tolerance for estimated error must be non negative");
                if((multiplier>5.0)||(multiplier<2.0))
                   smartmath_throw("maximum step-multiplier must be between 2 and 5");
                if(minstep_events<=0.0)
                   smartmath_throw("minimum step-size for events must be non negative");

            }

            /**
             * @brief ~base_stepsizecontrol deconstructor
             */
            virtual ~base_stepsizecontrol(){}

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
             * @param[out] estimated error
             * @return
             */
            virtual int integration_step(const double &ti, const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal, T &er) const = 0;

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
            int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> t_history, std::vector<int> (*g)(std::vector<T> x, double d)) const{

                x_history.clear();
                t_history.clear();

                int k;
                int check=0;
                std::vector<T> x(x0), xtemp(x0);
                std::vector<std::vector<T> > f;

                double factor=1.0, value=0.0, t=ti, h = (tend-ti)/nsteps;
                T er=0.0*x0[0];

                std::vector<int> events=g(x0,ti), events2=events;
                int m=events.size();             

                int i=0;
                while(sqrt(pow(t-ti,2))<sqrt(pow(tend-ti,2))){

                    if((h*h>m_maxstep_events*m_maxstep_events)&&(m_maxstep_events>0.0)){
                        if(h>0.0){
                            h=m_maxstep_events;
                        }
                        else{
                            h=-m_maxstep_events;
                        }   
                    }   

                    if(sqrt(pow(tend-t,2))<sqrt(h*h))
                        h=tend-t;

                    integration_step(t,m_control,h,x,f,xtemp,er);
                    
                    /* Step-size control */
                    error(er,value);
                    factor=pow(m_tol/value,1.0/(m_control+1.0));
                    if(value>m_tol){ // unsucessful step
                        h*=0.9*factor;  
                    }
                    else{ // sucessful step
                        /* Checking for the events */
                        events2=g(xtemp,t+h);
                        k=0;
                        check=0;        
                        while(k<m){
                            if(events2[k]-events[k]!=0){
                                check=1;
                                k=m;
                            }
                            else{
                                k++;
                            }
                        }

                        if(check==1){
                            if(sqrt(h*h)>m_minstep_events){
                                h*=0.5;
                            }
                            else{
                                tend=t+h; // saving the termination time    
                                t=tend; // trick to get out of the while loop   
                                x_history.push_back(xtemp);
                                t_history.push_back(t);
                                std::cout << "Propagation interrupted by terminal event at time " << tend << " after " << i << " steps" << std::endl;
                            }
                        }
                        else{
                            x=xtemp; // updating state
                            t+=h; // updating current time          
                            events=events2;             
                            x_history.push_back(x);
                            t_history.push_back(t); 
                            /* Step-size control */
                            if(factor>m_multiplier){
                                factor=m_multiplier;
                            }   
                            h*=factor; // updating step-size
                            i++; // counting the number of steps
                        }           
        
                    }
                }

                return 0;
            }

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
             * @param[out] vector with one 0
             * @return
             */
            static std::vector<int> dummy_event(std::vector<T> x, double d){

                std::vector<int> output(1,0);

                return output;

            }

            /**
             * @brief returns a double equal to the input for real numbers and something meaningful for polynomials
             *
             * @param[in] x estimated error
             * @param[out] val double equal to x for real numbers and something else for polynomials
             * @return
             */
            int error(const T &x, double &val) const{
                val=x;  
                return 0;
            }

            #ifdef ENABLE_SMARTUQ
                int error(const smartuq::polynomial::chebyshev_polynomial<double> &x, double &val) const{
                    val=x.get_range()[1];
                    return 0;
                }
                int error(const smartuq::polynomial::chebyshev_polynomial<float> &x, double &val) const{
                    val=x.get_range()[1];
                    return 0;
                }
                int error(const smartuq::polynomial::chebyshev_polynomial<long double> &x, double &val) const{
                    val=x.get_range()[1];
                    return 0;
                }                        
                int error(const smartuq::polynomial::taylor_polynomial<double> &x, double &val) const{
                    val=x.get_coeffs()[0];
                    return 0;
                }
                int error(const smartuq::polynomial::taylor_polynomial<float> &x, double &val) const{
                    val=x.get_coeffs()[0];
                    return 0;
                }
                int error(const smartuq::polynomial::taylor_polynomial<long double> &x, double &val) const{
                    val=x.get_coeffs()[0];
                    return 0;
                }            
            #endif            
	
        };

    }
}

#endif // SMARTMATH_BASE_STEPSIZECONTROL_H
