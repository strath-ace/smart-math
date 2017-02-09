/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_RK87_H
#define SMARTMATH_RK87_H

#include "base_integrator.h"
#include "../exception.h"
#include "../Events/base_event.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The Runge Kutta 8(7)-13 integrator scheme by Prince and Dormand
         *
         * The class model the Runge Kutta Felhberg integration scheme
         */
        template < class T >
        class rk87: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            double m_control;
            double m_tol;
            double m_multiplier;
            double m_minstep_events;
            double m_maxstep_events;

        public:
        	
        	using base_integrator<T>::integrate;

            /**
             * @brief rk87 constructor
             *
             * @param dyn
             * @param tolerance for error estimarion in step-size control
             * @param max multiplier for step-size control
             * @param min time step for events detection
             */
            rk87(const dynamics::base_dynamics<T> *dyn, const double tol=1.0e-7, const double multiplier=5.0, const double minstep_events=1.0e-4, const double maxstep_events=0.0): base_integrator<T>("Runge Kutta 8-7 variable step time", dyn), m_tol(tol), m_multiplier(multiplier), m_minstep_events(minstep_events), m_maxstep_events(maxstep_events)
            {
	            /** sanity checks **/
	            if(tol<=0.0)
                   smartmath_throw("tolerance for estimated error must be non negative");
                if((multiplier>5.0)||(multiplier<2.0))
                   smartmath_throw("maximum step-multiplier must be between 2 and 5");
	            if(minstep_events<=0.0)
                   smartmath_throw("minimum step for events must be non negative"); 

               m_control=8.0;

            }
            /**
              * @brief ~rk87 deconstructor
              */
            ~rk87(){}

            /**
             * @brief integrate method to integrate bewteen two given time steps, initial condition and step lenght
             *
             * The method implements the RK8(7) scheme to perform one integration step
             * @param[in] ti initial time
             * @param[in] h step size
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @param[out] er estimated error 
             * @return
             */
            int integration_step(const double &t, const double &h, const std::vector<T> &x, std::vector<T> &xfinal, T &er) const{
		
		        int n = x.size();
		        std::vector<T> xtemp(x), xbar(x), k1(x), k2(x), k3(x), k4(x), k5(x), k6(x), k7(x), k8(x), k9(x), k10(x), k11(x), k12(x), k13(x);
		        double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13;

		        t1 = t;
		        t2 = t + h/18.0;
		        t3 = t + h/12.0;
		        t4 = t + h/8.0;
		        t5 = t + h*5.0/16.0;
		        t6 = t + h*3.0/8.0;
		        t7 = t + h*59.0/400.0;
		        t8 = t + h*93.0/200.0;
		        t9 = t + h*5490023248.0/9719169821.0;
		        t10 = t + h*13.0/20.0;
		        t11 = t + h*1201146811.0/1299019798.0;
		        t12 = t + h;
		        t13 = t + h;

		        //* Evaluate k1 
		        m_dyn->evaluate(t1, x, k1);

		        //* Evaluate k2 
		        for(int j=0; j<n; j++)
		            xtemp[j] = x[j]+k1[j]*h/18.0;
		        m_dyn->evaluate(t2, xtemp, k2);

		        //* Evaluate k3 
		        for(int j=0; j<n; j++)
		            xtemp[j] = x[j]+k1[j]*h/48.0+k2[j]*h/16.0;
		        m_dyn->evaluate(t3, xtemp, k3);

		        //* Evaluate k4 
		        for(int j=0; j<n; j++)
		            xtemp[j] = x[j]+k1[j]*h/32.0+k3[j]*h*3.0/32.0;
		        m_dyn->evaluate(t4, xtemp, k4);

		        //* Evaluate k5
		        for(int j=0; j<n; j++)
		            xtemp[j] = x[j]+k1[j]*h*5.0/16.0-k3[j]*h*75.0/64.0+k4[j]*h*75.0/64.0;
		        m_dyn->evaluate(t5, xtemp, k5);		

		        //* Evaluate k6
		        for(int j=0; j<n; j++)
		            xtemp[j] = x[j]+k1[j]*h*3.0/80.0+k4[j]*h*3.0/16.0+k4[j]*h*3.0/20.0;
		        m_dyn->evaluate(t6, xtemp, k6);

		        //* Evaluate k7
		        for(int j=0; j<n; j++)
		            xtemp[j] = x[j]+k1[j]*h*29443841.0/614563906.0 +k4[j]*h*77736538.0/692538347.0 -k5[j]*h*28693883.0/1125000000.0 +k6[j]*h*23124283.0/1800000000.0;
		        m_dyn->evaluate(t7, xtemp, k7);

		        //* Evaluate k8
		        for(int j=0; j<n; j++)
		            xtemp[j] = x[j]+k1[j]*h*16016141.0/946692911.0 +k4[j]*h*61564180.0/158732637.0 +k5[j]*h*22789713.0/633445777.0
		             +k6[j]*h*545815736.0/2771057229.0 -k7[j]*h*180193667.0/1043307555.0;
		        m_dyn->evaluate(t8, xtemp, k8);

		        //* Evaluate k9
		        for(int j=0; j<n; j++)
		            xtemp[j] = x[j]+k1[j]*h*39632708.0/573591083.0 -k4[j]*h*433636366.0/683701615.0 -k5[j]*h*421739975.0/2616292301.0
		             +k6[j]*h*100302831.0/723423059.0 +k7[j]*h*790204164.0/839813087.0  +k8[j]*h* 800635310.0/3783071287.0;
		        m_dyn->evaluate(t9, xtemp, k9);

		        //* Evaluate k10
		        for(int j=0; j<n; j++)
		            xtemp[j] = x[j]+k1[j]*h*246121993.0/1340847787.0 -k4[j]*h*37695042795.0/15268766246.0-k5[j]*h*309121744.0/1061227803.0
		             -k6[j]*h*12992083.0/490766935.0 +k7[j]*h*6005943493.0/2108947869.0  +k8[j]*h*393006217.0/1396673457.0 +k9[j]*h*123872331.0/1001029789.0;
		        m_dyn->evaluate(t10, xtemp, k10);

		        //* Evaluate k11
		        for(int j=0; j<n; j++)
		            xtemp[j] = x[j]-k1[j]*h*1028468189.0/846180014.0 +k4[j]*h*8478235783.0/508512852.0  +k5[j]*h*1311729495.0/1432422823.0
		             -k6[j]*h*10304129995.0/1701304382.0 -k7[j]*h*48777925059.0/3047939560.0  +k8[j]*h*15336726248.0/1032824649.0 -k9[j]*h*45442868181.0/3398467696.0 +k10[j]*h*3065993473.0/597172653.0;
		        m_dyn->evaluate(t11, xtemp, k11);

		        //* Evaluate k12
		        for(int j=0; j<n; j++)
		            xtemp[j] = x[j]+k1[j]*h*185892177.0/718116043.0 -k4[j]*h*3185094517.0/667107341.0-k5[j]*h*477755414.0/1098053517.0
		             -k6[j]*h*703635378.0/230739211.0 +k7[j]*h*5731566787.0/1027545527.0  +k8[j]*h*5232866602.0/850066563.0 -k9[j]*h*4093664535.0/808688257.0 +k10[j]*h* 3962137247.0/1805957418.0
		             +k11[j]*h*65686358.0/487910083.0;
		        m_dyn->evaluate(t12, xtemp, k12);

		        //* Evaluate k13
		        for(int j=0; j<n; j++)
		            xtemp[j] = x[j]+k1[j]*h*403863854.0/491063109.0 -k4[j]*h*5068492393.0/434740067.0-k5[j]*h*411421997.0/543043805.0
		             +k6[j]*h*652783627.0/914296604.0 +k7[j]*h*11173962825.0/925320556.0  -k8[j]*h*13158990841.0/6184727034.0 +k9[j]*h*3936647629.0/1978049680.0 -k10[j]*h*160528059.0/685178525.0 
		             +k11[j]*h*248638103.0/1413531060.0 ;
		        m_dyn->evaluate(t13, xtemp, k13);

		        //* Return x(t+h) computed from Runge Kutta.
		        er=0.0;
		        for(int j=0; j<n; j++){
		            xbar[j] = x[j]+ (k1[j]*14005451.0/335480064.0 -k6[j]*59238493.0/1068277825.0 +k7[j]*181606767.0/758867731.0 +k8[j]* 561292985.0/797845732.0
		             -k9[j]*1041891430.0/1371343529.0 +k10[j]*760417239.0/1151165299.0  +k11[j]*118820643.0/751138087.0 -k12[j]*528747749.0/2220607170.0 + k13[j]/4.0)*h;
		            xtemp[j] = x[j] + (k1[j]*13451932.0/455176623.0 -k6[j]*808719846.0/976000145.0 +k7[j]*1757004468.0/5645159321.0 +k8[j]*656045339.0/265891186.0
		             -k9[j]*3867574721.0/1518517206.0 +k10[j]*465885868.0/322736535.0  +k11[j]*53011238.0/667516719.0 +k12[j]*2.0/45.0)*h;
		              er+=pow(xbar[j]-xtemp[j],2);
		        }
		        er=sqrt(er);

		        return 0;

            }

            /*
             * @brief integrate method to integrate bewteen two given time steps, initial condition and step lenght
             *
             * The method implements the RK8(7)-13 scheme to integrate with given initial time,
             * final time, initial state condition and initial guess for step-size
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &xfinal, std::vector<double> &t_history) const{

	            xfinal.clear();
	            t_history.clear();

                std::vector<T> x(x0), xtemp(x0);
	            double factor=1.0, value=0.0, t=ti, h = (tend-ti)/nsteps;
	            T er=0.0*x0[0];

	            int i=0;
                while(sqrt(pow(t-ti,2))<sqrt(pow(tend-ti,2))){

		            if(sqrt(pow(tend-t,2))<sqrt(h*h))
			            h=tend-t;

		            integration_step(t,h,x,xtemp,er);

		            /* Step-size control */
		            error(er,value);
		            factor=pow(m_tol/value,1.0/8.0);
		            if(value>m_tol){ // unsuccessful step
			            h*=0.9*factor; // correcting step	
			            //std::cout << i << std::endl;	
		            }
		            else{ // successful step
				            x=xtemp; // updating state	
				            t+=h; // updating current time
				            if(factor>m_multiplier)
					            factor=m_multiplier;	
				            i++; // counting number of integration steps	
				            xfinal.push_back(x);		
				            t_history.push_back(t);
		            }
	            }


	            return 0;
            }


            /**
             * @brief integrate method to integrate bewteen two given time steps, initial condition and initial guess for step-size
             *
             * The method implements the RK8(7)-13 scheme to integrate with given initial time,
             * final time, initial state condition and initial guess for step-size
             * @param[in] ti initial time instant
             * @param[out] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @param[in] event function
             * @return
             */
            int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal, std::vector<int> (*g)(std::vector<T> x, double d)) const{

	            xfinal.clear();

                int n = x0.size();
                int k;
                int check=0;

                std::vector<T> x(x0), xtemp(x0);

	            double factor=1.0, value=0.0, t=ti, h = (tend-ti)/nsteps;
	            T er=0.0*x0[0];

	            std::vector<int> events=g(x0,ti), events2=events;
	            int m=events.size();

	            int i=0;

	            if((h>m_maxstep_events)&&(m_maxstep_events>0.0)){
			            h=m_maxstep_events;
	            }

                while(sqrt(pow(t-ti,2))<sqrt(pow(tend-ti,2))){

		            if(sqrt(pow(tend-t,2))<sqrt(h*h))
			            h=tend-t;

		            integration_step(t,h,x,xtemp,er);

		            /* Step-size control */
		            error(er,value);
		            factor=pow(m_tol/value,1.0/(m_control+1.0));
		            if(value>m_tol){ // unsucessful step
			            h*=0.9*factor; 		
		            }
		            else{ // sucessful step
			            /* Checking for the events */
			            events2=g(xtemp,t+h);
			            k=0;
			            check=0;		
			            while(k<m){
				            if(events2[k]-events[k]!=0){
					            check=1;
					            k=m;
				            }
				            else{
					            k++;
				            }
			            }

			            if(check==1){
				            if(sqrt(h*h)>m_minstep_events){
					            h*=0.5;
				            }
				            else{
					            tend=t+h; // saving the termination time			
					            t=tend;
					            std::cout << "Propagation interrupted by terminal event at time " << tend << " after " << i << " steps" << std::endl;
				            }
			            }
			            else{
				            x=xtemp; // updating state
				            t+=h; // updating current time			
				            events=events2;					
				            /* Step-size control */
				            if(factor>m_multiplier){
					            factor=m_multiplier;
				            }	
				            h*=factor; // updating step-size
				            i++; // counting the number of steps
				            if((m_maxstep_events>0.0)&&(h>m_maxstep_events)){
						            h=m_maxstep_events;
				            }											
			            }				
		
		            }	
		            //std::cout << "current step is " << i << " with relative size of " << h/(tend-ti) << " and estimated error of " << value << std::endl;	
	            }

	            for(int j=0; j<n; j++)
	                xfinal.push_back(x[j]);

	            return 0;
            }

            /**
             * @brief integrate method to integrate bewteen two given time steps, initial condition and initial guess for step-size
             *
             * The method implements the RK8(7)-13 scheme to integrate with given initial time,
             * final time, initial state condition and initial guess for step-size
             * @param[in] ti initial time instant
             * @param[out] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @param[in] event function
             * @return
             */
            int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal, std::vector<T> &t_history, std::vector<events::base_event<T>*> &event_list) const{

	            xfinal.clear();
	            t_history.clear();

                int k;
                int check=0;

                std::vector<T> x(x0), xtemp(x0);

	            double factor=1.0, value=0.0, t=ti, h = (tend-ti)/nsteps;
	            T er=0.0*x0[0];

	            std::vector<int> events(event_list.size(), 0);
	            for(unsigned int index = 0; index < events.size(); ++index)
	            {
		            events[index] = event_list[index]->evaluate(ti, x0);
	            }
	            std::vector<int> events2 = events;
	            int m=events.size();

	            int i=0;

	            if((h>m_maxstep_events)&&(m_maxstep_events>0.0)){
			            h=m_maxstep_events;
	            }

                while(sqrt(pow(t-ti,2))<sqrt(pow(tend-ti,2))){

		            if(sqrt(pow(tend-t,2))<sqrt(h*h))
			            h=tend-t;

		            integration_step(t,h,x,xtemp,er);

		            /* Step-size control */
		            error(er,value);
		            factor=pow(m_tol/value,1.0/(m_control+1.0));
		            if(value>m_tol){ // unsucessful step
			            h*=0.9*factor; 		
		            }
		            else{ // sucessful step
			            /* Checking for the events */
			            for(unsigned int index = 0; index < event_list.size(); ++index)
			            {
				            events2[index] = event_list[index]->evaluate(t+h, xtemp);
			            }
			            k=0;
			            check=0;		
			            while(k<m){
				            if(events2[k]-events[k]!=0){
					            check=1;
					            k=m;
				            }
				            else{
					            k++;
				            }
			            }

			            if(check==1){
				            if(sqrt(h*h)>m_minstep_events){
					            h*=0.5;
				            }
				            else{
					            t_history.push_back(t+h); // saving the termination time					
					            t=tend; // trick to get out of the while loop
					            xfinal.push_back(x);				
					            for(unsigned int index = 0; index < m; ++index)
					            {
						            if(events2[index]-events[index]!=0)
							            event_list[index]->switch_trigger_on(t_history.back());
					            }
					            // std::cout << "Propagation interrupted by terminal event at time " << t_history.back() << " after " << i << " steps" << std::endl;
				            }
			            }
			            else{
				            x=xtemp; // updating state
				            xfinal.push_back(x);
				            t+=h; // updating current time	
				            t_history.push_back(t);		
				            events=events2;					
				            /* Step-size control */
				            if(factor>m_multiplier){
					            factor=m_multiplier;
				            }	
				            h*=factor; // updating step-size
				            i++; // counting the number of steps
				            if((m_maxstep_events>0.0)&&(h>m_maxstep_events)){
						            h=m_maxstep_events;
				            }											
			            }			
		
		            }	
		             //std::cout << "current step is " << i << " with relative size of " << h/(tend-ti) << " and estimated error of " << value << std::endl;	
	            }

	            return 0;
            }


             int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal, std::vector<events::base_event<T>*> &event_list) const{

	            std::vector<std::vector<T> > xf;
	            std::vector<double> t;

	            integrate(ti, tend, nsteps, x0, xf, t, event_list);

	            xfinal=xf.back();
	            tend=t.back();

	            return 0;
            }


            /**
             * @brief returns a double equal to the input for real numbers and something meaningful for polynomials
             *
             * 
             * @param[in] x estimated error
             * @param[out] val double equal to x for real numbers and something else for polynomials
             * @return
             */
            int error(const T &x, T &val) const{
	            val=x;
                return 0;
            }
            #ifdef ENABLE_SMARTUQ
                int error(const smartuq::polynomial::chebyshev_polynomial<double> &x, double &val) const{
	                val=x.get_coeffs()[0];
                return 0;
                }
                int error(const smartuq::polynomial::chebyshev_polynomial<float> &x, double &val) const{
	                val=x.get_coeffs()[0];
                return 0;
                }
                int error(const smartuq::polynomial::chebyshev_polynomial<long double> &x, double &val) const{
	                val=x.get_coeffs()[0];
                return 0;
                }                      
                int error(const smartuq::polynomial::taylor_polynomial<double> &x, double &val) const{
	                val=x.get_coeffs()[0];
                return 0;
                }
                int error(const smartuq::polynomial::taylor_polynomial<float> &x, double &val) const{
	                val=x.get_coeffs()[0];
                return 0;
                }
                int error(const smartuq::polynomial::taylor_polynomial<long double> &x, double &val) const{
	                val=x.get_coeffs()[0];
                return 0;
                }            
            #endif

        };

    }
}

#endif // SMARTMATH_RK87_H
