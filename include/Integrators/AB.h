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
            std::vector<double> m_beta;
            integrator::rk4<T> *m_initializer;

        public:

            using base_multistep<T>::integrate;
            
            /**
             * @brief Adam Bashforth constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             */
            AB(const dynamics::base_dynamics<T> *dyn, const int order=6): base_multistep<T>("Adam Bashforth integration scheme", dyn, order)
            {

                if((order<1)||(order>8))
                    smartmath_throw("order must be between 1 and 8");  

	            double prebeta[8][8]={
                {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                {-1.0/2.0,3.0/2.0,0.0,0.0,0.0,0.0,0.0,0.0},
                {5.0/12.0,-16.0/12.0,23.0/12.0,0.0,0.0,0.0,0.0,0.0},
                {-9.0/24.0,37.0/24.0,-59.0/24.0,55.0/24.0,0.0,0.0,0.0,0.0},
                {251.0/720.0,-1274.0/720.0,2616.0/720.0,-2774.0/720.0,1901.0/720.0,0.0,0.0,0.0},
                {-475.0/1440.0,2877.0/1440.0,-7298.0/1440.0,9982.0/1440.0,-7923.0/1440.0,4277.0/1440.0,0.0,0.0},
                {19087.0/60480.0,-134472.0/60480.0,407139.0/60480.0,-688256.0/60480.0,705549.0/60480.0,-447288.0/60480.0,198721.0/60480.0,0.0},
                {-36799.0/120960.0,295767.0/120960.0,-1041723.0/120960.0,2102243.0/120960.0,-2664477.0/120960.0,2183877.0/120960.0,-1152169.0/120960.0,434241.0/120960.0}
                };        
	            for(int i=0; i<m_order; i++){
		            m_beta.push_back(prebeta[m_order-1][i]);
	            }

                m_initializer = new integrator::rk4<T>(m_dyn);

            }

            /**
              * @brief ~AB deconstructor
              */
            ~AB(){}

            /**
             * @brief integration_step performs one integration step from the Adam Bashforth (order 6) method
             *
             * The method implements one step of the Adam Bashforth scheme to integrate with given initial time,
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
             * The method initializes via RK4 the Adam Bashforth scheme for an integration with step-size h starting at given initial time and condition 
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
	            
                /* Computing the initial saved steps */
	            m_dyn->evaluate(ti,x,dx);
	            fp.push_back(dx);
	            double t=ti;
	            for(int j=0; j<m-1; j++){
		            m_initializer->integration_step(t,-h,x,xp);
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

#endif // SMARTMATH_AB_H
