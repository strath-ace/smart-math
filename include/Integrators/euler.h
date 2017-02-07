#ifndef SMARTMATH_EULER_H
#define SMARTMATH_EULER_H

#include "base_integrator.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The Euler explicit integration scheme
         *
         * The class model the Euler explicit integration scheme
         */
        template < class T >
        class euler: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;

        public:
            /**
             * @brief euler constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             */
            euler(const dynamics::base_dynamics<T> *dyn): base_integrator<T>("Explicit Euler integration scheme", dyn){}

            /**
              * @brief ~euler deconstructor
              */
            ~euler(){}

            /**
             * @brief integration_step performs one integration step from the Euler explicit method
             *
             * The method implements one step of the explicit Euler scheme to integrate with given initial time,
             * final time, initial state condition(constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] h time step
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integration_step(const double &ti, const double &h, const std::vector<T> &x0, std::vector<T> &xfinal) const{
	
	            std::vector<T> dx=x0;
	            unsigned int l = x0.size();

	            m_dyn->evaluate(ti, x0, dx);

	            xfinal=x0;	
	            for(unsigned int j=0; j<l; j++){
		            xfinal[j] += h*dx[j];
	            }

	            return 0;
            }

            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps
             *
             * The method implements the explicit Euler scheme to integrate with given initial time,
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
             * The method implements the explicit Euler scheme to integrate with given initial time,
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

#endif // SMARTMATH_EULER_H
