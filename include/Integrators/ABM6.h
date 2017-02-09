/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_ABM6_H
#define SMARTMATH_ABM6_H

#include "base_integrator.h"
#include "AB6.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        template < class T >
        class ABM6: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            int m_order;
            std::vector<double> m_beta;

        public:

            using base_integrator<T>::integrate;

            /**
             * @brief Adam Bashforth Moulton 6 constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             */
            ABM6(const dynamics::base_dynamics<T> *dyn): base_integrator<T>("Adam Bashforth Moulton algorithm with order 6", dyn)
            {
	            m_order=6;
	            double prebeta[6]={27.0,-173.0,482.0,-798.0,1427.0,475.0};
	            for(int i=0; i<m_order; i++)
		            m_beta.push_back(prebeta[i]/1440.0);

            }

            /**
              * @brief ~ABM6 deconstructor
              */
            ~ABM6(){}
	        
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
	            std::vector<std::vector<T> > f, fp;

	            double t=ti, h = (tend-ti)/nsteps;

	            integrator::AB6<T> predictor(m_dyn);

	            predictor.initialize(ti,h,x0,f); // initializing predictor-corrector by initializing predictor

                for(int k=0; k<nsteps; k++){

                	predictor.integration_step(h,x,f,xp); // prediction 

                	/* Updating the saved steps */
                	fp=f;
		            for(int j=0; j<m_order-1; j++){
			            f[j]=fp[j+1];
	                }
	                m_dyn->evaluate(t+h, xp, dx);
	                f[m_order-1]=dx;

		            correction(h,x,f,xp); 

		            /* Updating the last saved step */
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
             * @brief correction performs the correction according to Adam Moulton 6 
             *
             * The method implements the correction step in the Adam Bashforth Moulton 6 algorithm
             * @param[in] h step-size
             * @param[in] x0 vector of initial states
             * @param[in] f vector of saved state vectors for multistep scheme
             * @param[out] xfinal vector of final states
             * @return
             */
            int correction(const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal) const{

	            xfinal=x0;
	            for(int i=0; i<x0.size(); i++){
		            for(int j=0; j<m_order; j++){
			            xfinal[i]+=h*m_beta[j]*f[j][i];
	                }
	            }

	            return 0;
            }

        };

    }
}

#endif // SMARTMATH_ABM6_H
