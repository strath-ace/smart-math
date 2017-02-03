#ifndef SMARTMATH_RK4_H
#define SMARTMATH_RK4_H

#include "base_integrator.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The Runge Kutta 4 integrator scheme
         *
         * The class model the Runge Kutta fourth order integration scheme
         */
        template < class T >
        class rk4: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;

        public:
            /**
             * @brief rk4 constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             */
            rk4(const dynamics::base_dynamics<T> *dyn) : base_integrator<T>("Runge Kutta 4 fixed step time", dyn){}

            /**
              * @brief ~rk4 deconstructor
              */
            ~rk4(){}

            /**
             * @brief performs one integration step from the RK4 method
             *
             * The method implements one step of the RK4 scheme to integrate with given initial time,
             * final time, initial state condition(constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] h time step
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integration_step(const double &ti, const double &h, const std::vector<T> &x0, std::vector<T> &xfinal) const{
	
	            std::vector<T> dx=x0;
	            std::vector<T> x_temp=x0, k1=x0, k2=x0, k3=x0, k4=x0;
	            unsigned int l = x0.size();
	            double t1, t2, t3, t4;
	            t1 = t;
	            t2 = t + h/2.0;
	            t3 = t + h/2.0;
	            t4 = t + h;

	            //* Evaluate k1 = f(x).
	            m_dyn->evaluate(t1, x0, k1);

	            //* Evaluate k2 = f(x+h/2*k1),
	            for(unsigned int j=0; j<l; j++)
	                x_temp[j] = x0[j]+k1[j]*h/2.0;
	            m_dyn->evaluate(t2, x_temp, k2);

	            //* Evaluate k3 = f(x+h/2*k2),
	            for(unsigned int j=0; j<l; j++)
	                x_temp[j] = x0[j]+k2[j]*h/2.0;
	            m_dyn->evaluate(t3, x_temp, k3);

	            //* Evaluate k4 = f(x+h*k3),
	            for(unsigned int j=0; j<l; j++)
	                x_temp[j] = x0[j]+k3[j]*h;
	            m_dyn->evaluate(t4, x_temp, k4);

	            //* Return x(t+h) computed from third-order Runge Kutta.
	            xfinal=x0;
	            for(unsigned int j=0; j<l; j++)
	                xfinal[j] +=  (k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j])*h/6.0;

	            return 0;
            }


            /**
             * @brief integrate method to integrate bewteen two given time steps, initial condition and step lenght
             *
             * The method implements the RK4 scheme to integrate with given initial time,
             * final time, initial state condition and number of steps (constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

	            std::vector<std::vector<T> > x_history;
	            std::vector<double> t_history;

	            integrate(ti,tend,nsteps,x0,x_history,t_history);

	            xfinal=x_history.back();

	            return 0;
            }
            
            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method implements the RK4 scheme to integrate with given initial time,
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

	            std::vector<T> dx=x0, x=x0, x_temp=x0;

	            double t=ti, h = (tend-ti)/nsteps;

                for(int i=0; i<nsteps; i++){
		            integration_step(t,h,x,x_temp);
		            t+=h;
		            x=x_temp;
		            t_history.push_back(t);
		            x_history.push_back(x);
	            }

	            return 0;
            }
        };

    }
}

#endif // SMARTMATH_RK4_H
