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

#include "base_multistep.h"
#include "AB.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        template < class T >
        class ABM: public base_multistep<T>
        {

        private:
            using base_multistep<T>::m_name;
            using base_multistep<T>::m_dyn;
            using base_multistep<T>::m_order;
            std::vector<double> m_beta_Moulton;
            integrator::AB<T> *m_predictor;

        public:

            using base_integrator<T>::integrate;

            /**
             * @brief Adam Bashforth Moulton constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             */
            ABM(const dynamics::base_dynamics<T> *dyn, const int order=8): base_multistep<T>("Adam Bashforth Moulton algorithm", dyn, order)
            {

                if((order<2)||(order>8))
                    smartmath_throw("ABM: order must be between 2 and 8"); 

                double prebeta[7][8]={
                {1.0/2.0,1.0/2.0,0.0,0.0,0.0,0.0,0.0,0.0},
                {-1.0/12.0,8.0/12.0,5.0/12.0,0.0,0.0,0.0,0.0,0.0},
                {1.0/24.0,-5.0/24.0,19.0/24.0,9.0/24.0,0.0,0.0,0.0,0.0},
                {-19.0/720.0,106.0/720.0,-264.0/720.0,646.0/720.0,251.0/720.0,0.0,0.0,0.0},
                {27.0/1440.0,-173.0/1440.0,482.0/1440.0,-798.0/1440.0,1427.0/1440.0,475.0/1440.0,0.0,0.0},
                {-863.0/60480.0,6312.0/60480.0,-20211.0/60480.0,37504.0/60480.0,-46461.0/60480.0,65112.0/60480.0,19087.0/60480.0,0.0},
                {1375.0/120960.0,-11351.0/120960.0,41499.0/120960.0,-88547.0/120960.0,123133.0/120960.0,-121797.0/120960.0,139849.0/120960.0,36799.0/120960.0}
                };        
                for(int i=0; i<m_order; i++){
                    m_beta_Moulton.push_back(prebeta[m_order-2][i]);
                }

                m_predictor = new integrator::AB<T>(m_dyn,m_order);

            }

            /**
              * @brief ~ABM deconstructor
              */
            ~ABM(){}
	        

            /**
             * @brief integration_step method to perform one step of integration
             *
             * The method implements one step of the Adam Bashforth Moulton scheme 
             * @param[in] t initial time for integration step 
             * @param[in] m order
             * @param[in] h step-size
             * @param[in] x0 vector of initial states
             * @param[in] f vector of saved state vectors for multistep scheme             
             * @param[out] xfinal vector of final states
             * @return
             */
            int integration_step(const double &t, const int &m, const double &h, const std::vector<T> &x0, std::vector<std::vector<T> > &f, std::vector<T> &xfinal) const{
                
                if(f.size()!=m)
                    smartmath_throw("INTEGRATION_STEP: wrong number of saved states in multistep integration"); 

                std::vector<T> dx=x0;    

                m_predictor->integration_step(t,m,h,x0,f,xfinal); // prediction 

                correction(h,x0,f,xfinal); 

                m_dyn->evaluate(t, xfinal, dx);
                f[m-1]=dx;

                return 0;
            }            
            
             /**
             * @brief correction performs the correction according to Adam Moulton 
             *
             * The method implements the correction step in the Adam Bashforth Moulton algorithm
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
			            xfinal[i]+=h*m_beta_Moulton[j]*f[j][i];
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

                m_predictor->initialize(m,ti,h,x0,f);

                return 0;
            }            

        };

    }
}

#endif // SMARTMATH_ABM_H
