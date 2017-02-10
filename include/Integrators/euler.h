#ifndef SMARTMATH_EULER_H
#define SMARTMATH_EULER_H

#include "base_rungekutta.h"
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
        class euler: public base_rungekutta<T>
        {

        private:
            using base_rungekutta<T>::m_name;
            using base_rungekutta<T>::m_dyn;

        public:

            using base_rungekutta<T>::integrate;
            
            /**
             * @brief euler constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             */
            euler(const dynamics::base_dynamics<T> *dyn): base_rungekutta<T>("Explicit Euler integration scheme", dyn){}

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
	
        };

    }
}

#endif // SMARTMATH_EULER_H
