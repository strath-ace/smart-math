/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_AB6_H
#define SMARTMATH_AB6_H

#include "base_multistep.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        template < class T >
        class AB6: public base_multistep<T>
        {

        private:
            using base_multistep<T>::m_name;
            using base_multistep<T>::m_dyn;
            using base_multistep<T>::m_order;
            std::vector<double> m_beta;

        public:

            using base_multistep<T>::integrate;
            
            /**
             * @brief Adam Bashforth 6 constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             */
            AB6(const dynamics::base_dynamics<T> *dyn): base_multistep<T>("Adam Bashforth 6 integration scheme", dyn, 6)
            {
                
	            double prebeta[6]={-475.0,2877.0,-7298.0,9982.0,-7923.0,4277.0};
	            for(int i=0; i<m_order; i++){
		            m_beta.push_back(prebeta[i]/1440.0);
	            }

            }

            /**
              * @brief ~AB6 deconstructor
              */
            ~AB6(){}

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
            int integration_step(const double &t, const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal) const{

                if(f.size()!=m)
                    smartmath_throw("wrong number of saved states in multistep integration"); 

	            xfinal=x0;
	            for(int i=0; i<x0.size(); i++){
		            for(int j=0; j<m; j++){
			            xfinal[i]+=h*m_beta[j]*f[j][i];
	                }
	            }

	            return 0;
            }
  

            /**
             * @brief initialize method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method initializes via RK4 the Adam Bashforth 6 scheme for an integration with step-size h starting at given initial time and condition 
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
	            for(int j=0; j<m; j++){
		            f.push_back(fp[m-j-1]);
	            }

	            return 0;
            }
    
        };

    }
}

#endif // SMARTMATH_AB6_H
