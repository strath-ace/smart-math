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

#include "base_multistep.h"
#include "AB6.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        template < class T >
        class ABM6: public base_multistep<T>
        {

        private:
            using base_multistep<T>::m_name;
            using base_multistep<T>::m_dyn;
            using base_multistep<T>::m_order;
            std::vector<double> m_beta;

        public:

            using base_integrator<T>::integrate;

            /**
             * @brief Adam Bashforth Moulton 6 constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             */
            ABM6(const dynamics::base_dynamics<T> *dyn): base_multistep<T>("Adam Bashforth Moulton algorithm with order 6", dyn, 6)
            {

	            double prebeta[6]={27.0,-173.0,482.0,-798.0,1427.0,475.0};
	            for(int i=0; i<m_order; i++)
		            m_beta.push_back(prebeta[i]/1440.0);

            }

            /**
              * @brief ~ABM6 deconstructor
              */
            ~ABM6(){}
	        

            /**
             * @brief integration_step method to perform one step of integration
             *
             * The method implements one step of the Adam Bashforth Moulton 6 scheme 
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

                std::vector<T> dx=x0;    
                integrator::AB6<T> predictor(m_dyn);

                predictor.integration_step(t,m,h,x0,f,xfinal); // prediction 

                /* Updating the saved steps */
                std::vector<std::vector<T> > fp=f;
                predictor.update_saved_steps(m,t+h,xfinal,fp);

                correction(h,x0,fp,xfinal); 

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

            /**
             * @brief initialize method to initialize integrator at initial time
             *
             * The method initializes via Adam Bashforth the Adam Bashforth Moulton scheme for an integration with step-size h starting at given initial time and condition 
             * @param[in] m number of saved steps
             * @param[in] ti initial time instant
             * @param[in] h step size
             * @param[in] x0 vector of initial states
             * @param[out] f vector of saved state vectors for multistep scheme
             * @return
             */     
            int initialize(const int &m, const double &ti, const double &h, const std::vector<T> &x0, std::vector<std::vector<T> > &f) const{

                integrator::AB6<T> predictor(m_dyn);    

                predictor.initialize(6,ti,h,x0,f);

                return 0;
            }            

        };

    }
}

#endif // SMARTMATH_ABM6_H
