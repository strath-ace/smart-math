/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_AB_H
#define SMARTMATH_AB_H

#include "base_multistep.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        template < class T >
        class AB: public base_multistep<T>
        {

        private:
            using base_multistep<T>::m_name;
            using base_multistep<T>::m_dyn;
            using base_multistep<T>::m_order;
            std::vector<double> m_gamma;

        public:

            using base_multistep<T>::integrate;

            /**
             * @brief Adam Bashforth constructor
             *
             * The integrator is initialized with the super class constructor. 
             * @param dyn
             * @param order chosen order for Adam Bashforth method
             */
            AB(const dynamics::base_dynamics<T> *dyn, const int order=6): base_multistep<T>("Adam Bashford integration scheme with user-defined order", dyn, order)
            {

	            if((order<1)||(order>8))
                	smartmath_throw("order must be between 1 and 8");    

	            double gamma[9]={1.0,1.0/2.0,5.0/12.0,3.0/8.0,251.0/720.0,95.0/288.0,19087.0/60480.0,5257.0/17280.0,1070017.0/3628800.0};
	            for(int i=0; i<=m_order; i++)
		            m_gamma.push_back(gamma[i]);
	
            }

            /**
              * @brief ~AB deconstructor
              */
            ~AB(){}

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
             * The method implements the Adam Bashforth scheme to integrate with given initial time,
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
	            double t=ti, h = (tend-ti)/nsteps;

	            initialize(m_order,ti,h,x0,f);

                for(int k=0; k<nsteps; k++){

		            backward_differences(f,m_order,Df);
		            integration_step(m_order,h,x,Df,xp);
		
                	/* Updating saved steps */
		            fp=f;
		            for(int j=0; j<m_order-1; j++)
			            f[j]=fp[j+1];
               		m_dyn->evaluate(t+h, xp, dx);
               		f[m_order-1]=dx;

               		/* Saving states */
		            t+=h;
		            x=xp;
		            t_history.push_back(t);
               		x_history.push_back(x);
	            }

	            return 0;
            }
    
            /**
             * @brief integration-step performs one integration step from the Adam Bashforth method
             *
             * The method implements one step of the Adam Bashforth scheme 
             * @param[in] m number of saved steps
             * @param[in] h time step
             * @param[in] x0 vector of initial states
             * @param[in] Df vector of finite differences vectors 
             * @param[out] xfinal vector of final states
             * @return
             */	
            int integration_step(const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &Df, std::vector<T> &xfinal) const{

	            xfinal=x0;
	            for(int i=0; i<x0.size(); i++){
		            for(int j=0; j<m; j++){		
			            xfinal[i]+=h*m_gamma[j]*Df[j][i];
	                }
	            }

	            return 0;
            }
            
            /**
             * @brief initialize method to initialize integrator at initial time
             *
             * The method initializes via RK4 the Adam Bashforth scheme for an integration with step-size h starting at given initial time and condition 
             * @param[in] m number of saved steps
             * @param[in] ti initial time instant
             * @param[in] h step size
             * @param[in] x0 vector of initial states
             * @param[out] f vector of saved state vectors for multistep scheme
             * @return
             */     
            int initialize(const int &m, const double &ti, const double &h, const std::vector<T> &x0, std::vector<std::vector<T> > &f) const{

	            f.clear();

	            std::vector<T> dx(x0), x(x0), xp(x0);
	            std::vector< std::vector<T> > fp;

	            integrator::rk4<T> RK(m_dyn); // Runge kutta schemed used for initialization (here RK4)

	            /* Computing the initial saved steps */
	            m_dyn->evaluate(ti,x,dx);
	            fp.push_back(dx);
	            double t=ti;
	            for(int j=0; j<m-1; j++){
		            RK.integration_step(t,-h,x,xp);
		            t-=h;
		            x=xp;
		            m_dyn->evaluate(t,x,dx);
		            fp.push_back(dx);
	            }

	            /* Putting the saved steps in the right order */
	            for(int j=0; j<m; j++)
		            f.push_back(fp[m-1-j]);

	            return 0;
            }
    
        };

    }
}

#endif // SMARTMATH_AB_H
