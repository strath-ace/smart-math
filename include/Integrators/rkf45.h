/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_RKF45_H
#define SMARTMATH_RKF45_H

#include "base_stepsizecontrol.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The Runge Kutta Felhberg (5) integrator scheme
         *
         * The class model the Runge Kutta Felhberg integration scheme
         */
        template < class T >
        class rkf45: public base_stepsizecontrol<T>
        {

        private:
            using base_stepsizecontrol<T>::m_name;
            using base_stepsizecontrol<T>::m_dyn;
            using base_stepsizecontrol<T>::m_tol;
            using base_stepsizecontrol<T>::m_multiplier;
            using base_stepsizecontrol<T>::m_control;
            using base_stepsizecontrol<T>::m_minstep_events;
            using base_stepsizecontrol<T>::m_maxstep_events;

        public:

            using base_stepsizecontrol<T>::integrate;
            using base_stepsizecontrol<T>::dummy_event;
            
            /**
             * @brief rkf45 constructor
             *
             * @param dyn
             * @param tolerance for error estimation in step-size control
             * @param max multiplier for step-size control
             * @param min time step for events detection
             */
            rkf45(const dynamics::base_dynamics<T> *dyn, const double tol=1.0e-7, const double multiplier=5.0, const double minstep_events=1.0e-4, const double maxstep_events=0.0): base_stepsizecontrol<T>("Runge Kutta 4-5 variable step time", dyn, tol, multiplier, minstep_events, maxstep_events)
            {

               m_control=4.0;

            }


            /**
              * @brief ~rkf45 deconstructor
              */
            ~rkf45(){}

            /**
             * @brief integrate method to integrate bewteen two given time steps, initial condition and initial guess for step-size
             *
             * The method implements the RKF4(5) scheme to perform one integration step
             * @param[in] ti initial time
             * @param[in] m method order
             * @param[in] h step size
             * @param[in] x0 vector of initial states
             * @param[in] f vector of saved state vectors (for multistep scheme only) 
             * @param[out] xfinal vector of final states
             * @param[out] er estimated error 
             * @return
             */
            int integration_step(const double &ti, const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal, T &er) const{
		
		        int n = x0.size();
		        std::vector<T> k1(x0), k2(x0), k3(x0), k4(x0), k5(x0), k6(x0), xbar(x0), xtemp(x0);
		        double t1, t2, t3, t4, t5, t6;
                xfinal=x0;

		        t1 = ti;
		        t2 = t1 + h/4.0;
		        t3 = t1 + h*3.0/8.0;
		        t4 = t1 + h*12.0/13.0;
		        t5 = t1 + h;
		        t6 = t1 + h/2.0;

		        //* Evaluate k1 
		        m_dyn->evaluate(t1, x0, k1);

		        //* Evaluate k2 
		        for(int j=0; j<n; j++)
		            xtemp[j] = x0[j]+k1[j]*h/4.0;
		        m_dyn->evaluate(t2, xtemp, k2);

		        //* Evaluate k3 
		        for(int j=0; j<n; j++)
		            xtemp[j] = x0[j]+k1[j]*h*3.0/32.0+k2[j]*h*9.0/32.0;
		        m_dyn->evaluate(t3, xtemp, k3);

		        //* Evaluate k4 
		        for(int j=0; j<n; j++)
		            xtemp[j] = x0[j]+k1[j]*h*1932.0/2197.0-k2[j]*h*7200.0/2197.0+k3[j]*h*7296.0/2197.0;
		        m_dyn->evaluate(t4, xtemp, k4);

		        //* Evaluate k5
		        for(int j=0; j<n; j++)
		            xtemp[j] = x0[j]+k1[j]*h*439.0/216.0-k2[j]*h*8.0+k3[j]*h*3680.0/513.0-k4[j]*h*845.0/4104.0;
		        m_dyn->evaluate(t5, xtemp, k5);

		        //* Evaluate k6
		        for(int j=0; j<n; j++)
		            xtemp[j] = x0[j]-k1[j]*h*8.0/27.0+k2[j]*h*2.0-k3[j]*h*3544.0/2565.0+k4[j]*h*1859.0/4104.0-k4[j]*h*11.0/40.0;
		        m_dyn->evaluate(t6, xtemp, k6);

		        //* Return x(t+h) computed from fourth-order Runge Kutta.
		        er=0.0*x0[0];
		        for(int j=0; j<n; j++){
		            xbar[j] += (k1[j]*16.0/135.0+k3[j]*6656.0/12825.0+k4[j]*28561.0/56430.0-k5[j]*9.0/50.0+k6[j]*2.0/55.0)*h;
		            xfinal[j] += (k1[j]*25.0/216.0+k3[j]*1408.0/2565.0+k4[j]*2197.0/4104.0-k5[j]/5.0)*h;
		            er+=pow(xbar[j]-xfinal[j],2);
		        }
		        er=sqrt(er);

	        return 0;
        }


        };

    }
}

#endif // SMARTMATH_RKF45_H
