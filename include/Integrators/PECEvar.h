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

#include "base_integrationwevent.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        template < class T >
        class PECEvar: public base_integrationwevent<T>
        { 

        private:
            using base_integrationwevent<T>::m_name;
            using base_integrationwevent<T>::m_dyn;
            using base_integrationwevent<T>::m_minstep_events;
            using base_integrationwevent<T>::m_maxstep_events;  
            using base_integrationwevent<T>::m_event_list;          
            double m_tol;
            double m_multiplier;
            int m_order_max, m_order_min;
            std::vector<double> m_gamma_Bashforth, m_gamma_Moulton;
            integrator::AB<T> *m_initializer;

        public:

            using base_integrationwevent<T>::integrate;
            using base_integrationwevent<T>::set_event_list;
        	
            /**
             * @brief Adam Bashforth Moulton with variable order and stepsize constructor
             *
             * The constructor specifically initializes a tolerance for integration error, a maximum multiplier to increase the stepsize as well as a minimum and maximum order  
             * @param dyn pointer to a base_dynamics object
             * @param order_max maximum order for predictor-corrector
             * @param order_min minimum order for predictor-corrector               
             * @param tol threshold used for acceptable estimated error
             * @param multiplier factor used to increase step-sized when judged necessary            
             * @param m_minstep_events minimum step-size to detect an event
             * @param m_maxstep_events maximum step-size     
             */
            PECEvar(const dynamics::base_dynamics<T> *dyn, const int order_min=4, const int order_max=8, const double tol=1.0e-7, const double multiplier=5.0, const double minstep_events=1.0e-4, const double maxstep_events=0.0): base_integrationwevent <T>("Adam Bashforth Moulton integration scheme with variable order", dyn, minstep_events, maxstep_events), m_tol(tol), m_multiplier(multiplier), m_order_max(order_max), m_order_min(order_min)
            {

	            if((order_min<1)||(order_min>8))
                	smartmath_throw("PECEVAR: minimum order must be between 1 and 8");    
	            if((order_max<1)||(order_max>8))
                	smartmath_throw("PECEVAR: maximum order must be between 1 and 8");
                if(order_min>order_max)
                    smartmath_throw("PECEVAR: minimum order must be smaller or equal to maximum order");                

	            double gamma_Bashforth[9]={1.0,-1.0/2.0,-1.0/12.0,-1.0/24.0,-19.0/720.0,-3.0/160.0,-863.0/60480.0,-275.0/24192.0,-33953.0/3628800.0};
	            for(int i=0; i<=m_order_max; i++)
		            m_gamma_Bashforth.push_back(gamma_Bashforth[i]);

                double gamma_Moulton[9]={1.0,-1.0/2.0,-1.0/12.0,-1.0/24.0,-19.0/720.0,-3.0/160.0,-863.0/60480.0,-275.0/24192.0,-33953.0/3628800.0};
                for(int i=0; i<=m_order_max; i++)
                    m_gamma_Moulton.push_back(gamma_Moulton[i]);

                m_initializer = new integrator::AB<T>(m_dyn,m_order_max);

            }

            /**
              * @brief ~PECEvar deconstructor
              */
            ~PECEvar(){}


            /**
             * @brief integrate method to integrate between two given time steps, with initial condition and initial guess for step-size while handling events
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
            int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history, std::vector<int> (*g)(std::vector<T> x, double d)) const{

	            t_history.clear();
	            x_history.clear();

                int i;
                int check=0;
                double factor;	
                
	            std::vector<T> x(x0),xp(x0);
	            std::vector<std::vector<T> > f_max, f;
	            T er;
	            int m=m_order_min, mold=m_order_min;
	            double t=ti, h = (tend-ti)/nsteps;
	            double value;

                std::vector<int> events, events2;
                if(m_event_list.size()==0)
                {
                    events=g(x0,ti);
                }
                else{
                    events=std::vector<int>(m_event_list.size(),0);
                    for(unsigned int index = 0; index < m_event_list.size(); ++index)
                    {   
                        events[index] = m_event_list[index]->evaluate(ti, x0);
                    }
                }
                events2=events;
                int q=events.size();  
                
	            m_initializer->initialize(m_order_max,ti,h,x0,f_max);

	            int k=0;
                while(sqrt(pow(t-ti,2))<sqrt(pow(tend-ti,2))){

                    if((h*h>m_maxstep_events*m_maxstep_events)&&(m_maxstep_events>0.0)){
                        if(h>=0.0){
                            h=m_maxstep_events;
                        }
                        else{
                            h=-m_maxstep_events;  
                        }   
                        m_initializer->initialize(m_order_max,t,h,x,f_max);
                    }  

		            if(sqrt(pow(tend-t,2))<sqrt(h*h)){
			            h=tend-t;
			            m_initializer->initialize(m_order_max,t,h,x,f_max);	
		            }
			
		            f.clear();
		            for(int j=0; j<m; j++)
			            f.push_back(f_max[m_order_max-m+j]);  	

                	integration_step(t,m,h,x,f,xp,er);
		
		            /* Step-size and order control */
		            value=evaluate_squarerootintegrationerror(er);
		            factor=pow(m_tol/value,1.0/(double(m)+1.0));
		            if(value>m_tol){ // unsucessful step
			            if(m<m_order_max){
				            mold=m;
				            m++;
				            //std::cout << "increasing order" << std::endl;
			            }
			            else{				
				            h*=0.9*factor;			
				            m_initializer->initialize(m_order_max,t,h,x,f_max);
				            //std::cout << "decreasing time-step" << std::endl;	
			            }
		            }
		            else{ // sucessful step
			            /* Checking for the events */
                        if(m_event_list.size()==0)
                        {
                            events2=g(xp,t+h);
                        }
                        else{
                            for(unsigned int index = 0; index < q; ++index)
                            {
                                events2[index] = m_event_list[index]->evaluate(t+h, xp);
                            }
                        }			            
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
					            m_initializer->initialize(m_order_max,t,h,x,f_max);
					            //std::cout << "decreasing time-step for event detection" << std::endl;
				            }
				            else{
					            tend=t+h; // saving the termination time	
					            t=tend; // trick to get out of the while loop	
					            x_history.push_back(xp);
					            t_history.push_back(t);
                                if(this->m_comments)
					                std::cout << "Propagation interrupted by terminal event at time " << tend << " after " << k << " steps" << std::endl;
				            }
			            }
			            else{
				            k++; // counting the number of steps
				            x=xp; // updating state
				            t+=h; // updating current time	
                            m_initializer->update_saved_steps(m_order_max,t,x,f_max);	
				            events=events2;				
				            x_history.push_back(x);
				            t_history.push_back(t);	
                            if(m_event_list.size()!=0)
                            {
                                for(unsigned int index = 0; index < q; ++index)
                                {
                                    if(events2[index]-events[index]!=0)
                                        m_event_list[index]->switch_trigger_on(t_history.back());
                                }
                            }                             
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
					            m_initializer->initialize(m_order_max,t,h,x,f_max);
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
                	smartmath_throw("INTEGRATION_STEP: wrong number of saved states for multistep integration"); 
     	
	            std::vector<T> x(x0), dx(x0);

                half_step(m,h,x0,f,m_gamma_Bashforth,x); // prediction

                std::vector<std::vector<T> > fp=f;
                for(int j=0; j<m-1; j++){
                    fp[j]=f[j+1];
                }
                m_dyn->evaluate(t+h, x, dx);
                fp[m-1]=dx;

	            half_step(m,h,x0,fp,m_gamma_Moulton,xfinal); // correction

	            unsigned int l=x.size();
	            er=0.0;
	            for(unsigned int i=0; i<l; i++)
			            er+=pow(xfinal[i]-x[i],2);
	            er=sqrt(er);	

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
            int backward_differences(const std::vector<std::vector<T> > &f, const int &m, std::vector<std::vector<T> > &Df) const{

                if(f.size()!=m)
                    smartmath_throw("BACKWARD_DIFFERENCES: wrong number of saved states for multistep integration"); 

                Df.clear();
                Df.push_back(f[m-1]);

                if(m>1){
                    std::vector<std::vector<T> > fp,Dfp;
                    for(int j=0; j<m-1; j++)
                        fp.push_back(f[j]);
                    backward_differences(fp,m-1,Dfp); // recursive call
                    for(int j=1; j<m; j++){
                        Df.push_back(Df[j-1]);
                        for(int i=0; i<f[0].size(); i++){
                            Df[j][i]-=Dfp[j-1][i];
                        }
                    }
                }

                return 0;
            }                                                     


             /**
             * @brief correction method to integrate between two given time steps, initial condition and number of steps
             *
             * The method implements either the prediction or the correction step in the Adam Bashforth Moulton scheme 
             * @param[in] m order
             * @param[in] h step-size
             * @param[in] x0 vector of initial states
             * @param[in] Df vector of finite differences vectors            
             * @param[out] xfinal vector of final states
             * @return
             */
            int half_step(const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, const std::vector<double> &gamma, std::vector<T> &xfinal) const{

                if(f.size()!=m)
                    smartmath_throw("HALF_STEP: wrong number of previously saved states for multistep integration");  

                std::vector<std::vector<T> > Df;
                backward_differences(f,m,Df);

                xfinal=x0;
                for(int i=0; i<x0.size(); i++){
                    for(int j=0; j<m; j++){     
                        xfinal[i]+=h*gamma[j]*Df[j][i];
                    }
                }

                return 0;
            }
      	
        };

    }
}

#endif // SMARTMATH_PECEVAR_H
