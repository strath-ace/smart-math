/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_ABM_H
#define SMARTMATH_ABM_H

#include "base_integrator.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        template < class T >
        class ABM: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            int m_order;
            std::vector<double> m_gamma;

        public:

        using base_integrator<T>::integrate;

            /**
             * @brief Adam Bashforth Moulton constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             * @param order order for Adam Bashforth Moulton scheme        
             */
            ABM(const dynamics::base_dynamics<T> *dyn, const int order=6): base_integrator<T>("Adam Bashforth Moulton integration scheme with user-defined order", dyn), m_order(order)
            {

	            if((order<1)||(order>8))
                	smartmath_throw("order must be between 1 and 8");    

	            double gamma[9]={1.0,-1.0/2.0,-1.0/12.0,-1.0/24.0,-19.0/720.0,-3.0/160.0,-863.0/60480.0,-275.0/24192.0,-33953.0/3628800.0};
	            for(int i=0; i<=m_order; i++)
		            m_gamma.push_back(gamma[i]);

            }

            /**
              * @brief ~ABM deconstructor
              */
            ~ABM(){}

            /**
             * @brief computes backward differences 
             *
             * The method computes the backward differences for the Adam Moulton scheme
             * @param[in] f vector of saved state vectors for multistep scheme
             * @param[in] m order
             * @param[out] Df vector of backward differences
             * @return
             */ 
            int backward_differences(const std::vector<std::vector<T> > &f, const int &m, std::vector<std::vector<T> > &Df) const{

	            if(f.size()!=m)
                	smartmath_throw("wrong number of saved states in multistep integration"); 

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
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method implements the Adam Bashforth Moulton 6 scheme to integrate with given initial time,
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

	            std::vector<T> x(x0),xp(x0),dx(x0);
	            std::vector<std::vector<T> > f, fp, Df;
	            integrator::AB<T> predictor(m_dyn, m_order);	
	            double t=ti, h = (tend-ti)/nsteps;

	            predictor.initialize(m_order,ti,h,x0,f); // initializing predictor-corrector by initializing predictor

                for(int k=0; k<nsteps; k++){

                	integration_step(t,m_order,h,x,f,xp);
		            x=xp;

		            /* Updating the saved steps */
		            fp=f;
		            for(int j=0; j<m_order-1; j++)
			            f[j]=fp[j+1];
               		m_dyn->evaluate(t+h, xp, dx);
               		f[m_order-1]=dx;

               		/* Saving states */
		            t+=h;
		            t_history.push_back(t);
               		x_history.push_back(x);
	            }

	            return 0;
            }


             /**
             * @brief correction method to integrate between two given time steps, initial condition and number of steps
             *
             * The method implements the correction step in the Adam Bashforth Moulton scheme 
             * @param[in] m order
             * @param[in] h step-size
             * @param[in] x0 vector of initial states
             * @param[in] Df vector of finite differences vectors            
             * @param[out] xfinal vector of final states
             * @return
             */
            int correction(const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &Df, std::vector<T> &xfinal) const{

	            xfinal=x0;
	            for(int i=0; i<x0.size(); i++){
		            for(int j=0; j<m; j++){		
			            xfinal[i]+=h*m_gamma[j]*Df[j][i];
	                }
	            }

	            return 0;
            }
            
            /**
             * @brief integration_step method to perform one step of integration
             *
             * The method implements one step of the Adam Bashforth Moulton scheme 
             * @param[in] ti initial time instant
             * @param[in] m order
             * @param[in] h step-size
             * @param[in] x0 vector of initial states
             * @param[in] f vector of saved state vectors for multistep scheme             
             * @param[out] xfinal vector of final states 
             * @return
             */
            int integration_step(const double &t, const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal) const{
                	
	            if(f.size()!=m)
                	smartmath_throw("wrong number of saved states in multistep integration"); 

	            integrator::AB<T> predictor(m_dyn, m);	    	
	            std::vector<T> x(x0), dx(x0);
	            std::vector<std::vector<T> > Df;

	            /* Prediction */
                predictor.backward_differences(f,m,Df);
                predictor.integration_step(m,h,x0,Df,x);

	            /* Updating the saved steps */
               	std::vector<std::vector<T> > fp=f;
	            for(int j=0; j<m-1; j++)
		            fp[j]=f[j+1];
               	m_dyn->evaluate(t+h, x, dx);
	            fp[m-1]=dx;

	            /* Correction */
	            backward_differences(fp,m,Df);
	            correction(m,h,x0,Df,xfinal);	

	            return 0;
            }
             	
        };

    }
}

#endif // SMARTMATH_ABM_H
