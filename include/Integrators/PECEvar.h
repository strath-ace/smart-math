/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_PECEVAR_H
#define SMARTMATH_PECEVAR_H

#include "base_stepsizecontrol.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        template < class T >
        class PECEvar: public base_stepsizecontrol<T>
        { 

        private:
            using base_stepsizecontrol<T>::m_name;
            using base_stepsizecontrol<T>::m_dyn;
            using base_stepsizecontrol<T>:: m_tol;
            using base_stepsizecontrol<T>:: m_multiplier;
            using base_stepsizecontrol<T>:: m_control;
            using base_stepsizecontrol<T>::m_minstep_events;
            using base_stepsizecontrol<T>::m_maxstep_events;
            int m_order_max, m_order_min;
            std::vector<double> m_gamma;

        public:

            using base_stepsizecontrol<T>::integrate;
        	
            /**
             * @brief Adam Bashforth Moulton constructor
             *
             * The integrator is initialized with the super class constructor. 
             * @param order_max maximum order for predictor-corrector
             * @param order_min minimum order for predictor-corrector         
             */
            PECEvar(const dynamics::base_dynamics<T> *dyn, const int order_max=8, const int order_min=4, const double tol=1.0e-7, const double multiplier=5.0, const double minstep_events=1.0e-4, const double maxstep_events=0.0): base_stepsizecontrol <T>("Adam Bashforth Moulton integration scheme with variable order", dyn, tol, multiplier, minstep_events, maxstep_events), m_order_max(order_max), m_order_min(order_min)
            {

	            if(order_min>order_max)
                	smartmath_throw("minimum order must be smaller or equal to maximum order");
	            if((order_min<1)||(order_min>8))
                	smartmath_throw("minimum order must be between 1 and 8");    
	            if((order_max<1)||(order_max>8))
                	smartmath_throw("maximum order must be between 1 and 8");

	            double gamma[9]={1.0,-1.0/2.0,-1.0/12.0,-1.0/24.0,-19.0/720.0,-3.0/160.0,-863.0/60480.0,-275.0/24192.0,-33953.0/3628800.0};
	            for(int i=0; i<=m_order_max; i++)
		            m_gamma.push_back(gamma[i]);

            }

            /**
              * @brief ~ABM deconstructor
              */
            ~PECEvar(){}


            /**
             * @brief integrate method to integrate bewteen two given time steps, with initial condition and initial guess for step-size while handling events
             *
             * The method implements the ABM scheme to integrate with given initial time,
             * final time, initial state condition and initial guess for step-size
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps initial guess for number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of intemrediate states
             * @param[out] t_history vector of intermediate times
             * @param[in] event function             
             * @return
             */
            int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> t_history, std::vector<int> (*g)(std::vector<T> x, double d)) const{

	            t_history.clear();
	            x_history.clear();

                int i;
                int check=0;
                double factor;	

                integrator::AB<T> predictor(m_dyn,m_order_max); 
	            std::vector<T> x(x0),xp(x0),dx(x0);
	            std::vector<std::vector<T> > f_max, f, Df;
	            T er;
	            int m=m_order_min, mold=m_order_min;
	            double t=ti, h = (tend-ti)/nsteps;
	            double value;

	            std::vector<int> events=g(x0,ti), events2=events;
	            int q=events.size();

	            predictor.initialize(m_order_max,ti,h,x0,f_max);

	            int k=0;
                while(sqrt(pow(t-ti,2))<sqrt(pow(tend-ti,2))){

                    if((h*h>m_maxstep_events*m_maxstep_events)&&(m_maxstep_events>0.0)){
                        if(h>=0.0){
                            h=m_maxstep_events;
                        }
                        else{
                            h=-m_maxstep_events;  
                        }   
                        predictor.initialize(m_order_max,t,h,x,f_max);
                    }  

		            if(sqrt(pow(tend-t,2))<sqrt(h*h)){
			            h=tend-t;
			            predictor.initialize(m_order_max,t,h,x,f_max);	
		            }
			
		            f.clear();
		            for(int j=0; j<m; j++)
			            f.push_back(f_max[m_order_max-m+j]);  	

                	integration_step(t,m,h,x,f,xp,er);
		
		            /* Step-size and order control */
		            error(er,value);
		            factor=pow(m_tol/value,1.0/(double(m)+1.0));
		            if(value>m_tol){ // unsucessful step
			            if(m<m_order_max){
				            mold=m;
				            m++;
				            //std::cout << "increasing order" << std::endl;
			            }
			            else{				
				            h*=0.9*factor;			
				            predictor.initialize(m_order_max,t,h,x,f_max);
				            //std::cout << "decreasing time-step" << std::endl;	
			            }
		            }
		            else{ // sucessful step
			            /* Checking for the events */
			            events2=g(xp,t+h);
			            i=0;
			            check=0;		
			            while(i<q){
				            if(events2[i]-events[i]!=0){
					            check=1;
					            i=q;
				            }
				            else{
					            i++;
				            }
			            }
			            if(check==1){
				            if(sqrt(h*h)>m_minstep_events){
					            h*=0.5;
					            predictor.initialize(m_order_max,t,h,x,f_max);
					            //std::cout << "decreasing time-step for event detection" << std::endl;
				            }
				            else{
					            tend=t+h; // saving the termination time	
					            t=tend; // trick to get out of the while loop	
					            x_history.push_back(xp);
					            t_history.push_back(t);
					            std::cout << "Propagation interrupted by terminal event at time " << tend << " after " << k << " steps" << std::endl;
				            }
			            }
			            else{
				            k++; // counting the number of steps
				            x=xp; // updating state
				            t+=h; // updating current time	
                            predictor.update_saved_steps(m_order_max,t,x,f_max);	
				            events=events2;				
				            x_history.push_back(x);
				            t_history.push_back(t);	
				            /* Step-size and order control */
				            if((m>m_order_min)&&(mold!=m-1)){
					            mold=m;
					            m--;
					            //std::cout << "decreasing order" << std::endl;
				            }
				            else if((m==m_order_min)&&(factor>=2.0)){
					            if(factor>m_multiplier){
						            factor=m_multiplier;
					            }	
					            h*=factor; // updating step-size
					            predictor.initialize(m_order_max,t,h,x,f_max);
					            //std::cout << "increasing time-step" << std::endl;
				            }
				            else{
					            //std::cout << "keeping both time-step and order constant" << std::endl;
				            }
				            //std::cout << k << std::endl;
			            }	
		            }	
	            }

	            return 0;
            }

            
            /**
             * @brief integration_step method to perform one step of integration
             *
             * The method implements one step of the Adam Bashforth Moulton scheme 
             * @param[in] t initial time for integration step 
             * @param[in] m method order
             * @param[in] h step-size
             * @param[in] x0 vector of initial states
             * @param[in] f vector of saved state vectors (for multistep scheme only)            
             * @param[out] xfinal vector of final states
             * @param[out] er estimated error
             * @return
             */
            int integration_step(const double &t, const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal, T &er) const{
                	
	            if(f.size()!=m)
                	smartmath_throw("wrong number of saved states in multistep integration"); 

	            integrator::AB<T> predictor(m_dyn, m);	
                integrator::ABM<T> corrector(m_dyn, m);      	
	            std::vector<T> x(x0);

                predictor.integration_step(t,m,h,x0,f,x);

               	std::vector<std::vector<T> > fp=f;
                predictor.update_saved_steps(m,t+h,x,fp);

	            corrector.correction(m,h,x0,fp,xfinal);

	            unsigned int l=x.size();
	            er=0.0;
	            for(unsigned int i=0; i<l; i++)
			            er+=pow(xfinal[i]-x[i],2);
	            er=sqrt(er);	

	            return 0;
            }                                             
      	
        };

    }
}

#endif // SMARTMATH_PECEVAR_H
