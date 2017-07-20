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
#include "bulirschstoer.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The AB class is an implementation of the Adam-Bashforth algorithm 
         *
         * The AB class is an implementation of a multistep integrator namely the Adam-Bashforth algorithm with fixed step-size 
         */
        template < class T >
        class AB: public base_multistep<T>
        {

        private:
            using base_multistep<T>::m_name;
            using base_multistep<T>::m_dyn;
            using base_multistep<T>::m_order;
            /**
             * @brief m_beta coefficients used in integration step
             */             
            std::vector<double> m_beta;
            /**
             * @brief m_initializerRK Runge-Kutta scheme used to initialize Adam-Bashforth
             */             
            integrator::rk4<T> *m_initializerRK;
            /**
             * @brief m_initializerBS Bulirsch-Stoer scheme used to initialize Adam-Bashforth
             */             
            integrator::bulirschstoer<T> *m_initializerBS;
            /**
             * @brief m_init boolean defining the type of initializer (true is B-S, false is R-K)
             */             
            bool m_init;

        public:

            using base_multistep<T>::integrate;
            
            /**
             * @brief Adam Bashforth constructor
             *
             * The integrator is initialized with the super class constructor. The user can choose the order of the method, the default value being 8 
             * @param dyn pointer to the dynamical system to be integrated
             * @param order order of the method
             */
            AB(const dynamics::base_dynamics<T> *dyn, const int order=8, const bool init = false): base_multistep<T>("Adam Bashforth integration scheme", dyn, order)
            {
                if((order<1)||(order>8))
                    smartmath_throw("AB: order must be between 1 and 8");  

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
                for(int i=0; i<m_order; i++)
                    m_beta.push_back(prebeta[m_order-1][i]);

                m_initializerRK = new integrator::rk4<T>(m_dyn);
                m_initializerBS = new integrator::bulirschstoer<T>(m_dyn, (order + 1) / 2);

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
            int integration_step(const double &t, const int &m, const double &h, const std::vector<T> &x0, std::vector<std::vector<T> > &f, std::vector<T> &xfinal) const{

                if(f.size()!= static_cast<unsigned int>(m))
                    smartmath_throw("INTEGRATION_STEP: wrong number of saved states for multistep integration"); 

                xfinal=x0;
                int size_x0 = x0.size();
                for(int i=0; i<size_x0; i++)
                {
                    for(int j=0; j<m; j++)
                        xfinal[i]+=h*m_beta[j]*f[j][i];
                }

                update_saved_steps(m,t+h,xfinal,f);

                return 0;
            }
  

            /**
             * @brief initialize method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method initializes with another integrator the Adam Bashforth scheme for an integration with step-size h starting at given initial time and conditions 
             * @param[in] m order
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
                for(int j=0; j<m-1; j++)
                {
                    if(m_init)
                        m_initializerBS->integration_step(t,-h,x,xp);
                    else
                        m_initializerRK->integration_step(t,-h,x,xp);
                    t-=h;
                    x=xp;
                    m_dyn->evaluate(t,x,dx);
                    fp.push_back(dx);
                }
                
                /* Putting the saved steps in the right order */
                for(int j=0; j<m; j++)
                    f.push_back(fp[m-j-1]);

                return 0;
            }
    
            /**
             * @brief update_saved_steps method to update saved integration steps
             *
             * The method updates the saved steps
             * @param[in] m number of saved steps
             * @param[in] t time of last state to save
             * @param[in] x vector of states at time t
             * @param[out] f vector of saved state vectors for multistep scheme
             * @return
             */     
            int update_saved_steps(const int &m, const double &t, const std::vector<T> &x, std::vector<std::vector<T> > &f) const{

                if(f[0].size()!=x.size())
                    smartmath_throw("UPDATE_SAVED_STEPS: wrong number of previously saved states for multistep integration"); 

                std::vector<T> dx=x;
                std::vector<std::vector<T> > fp=f;
                for(int j=0; j<m-1; j++)
                    f[j]=fp[j+1];
                m_dyn->evaluate(t, x, dx);
                f[m-1]=dx;

                return 0;
            }
    
        };

    }
}

#endif // SMARTMATH_AB_H
