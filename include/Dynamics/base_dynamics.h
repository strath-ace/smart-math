#ifndef SMARTMATH_BASE_DYNAMICS_H
#define SMARTMATH_BASE_DYNAMICS_H

#include <vector>
#include "../config.h"
#include "../exception.h"

#ifdef ENABLE_SMARTUQ
#include "smartuq.h"
#endif

namespace smartmath
{
    namespace dynamics {

        /**
         * @brief The base_dynamics class is an abstract class. Any dynamics added to the toolbox needs to inherit from it and implement the method evaluate()
         *
         * The base_dynamics class is an abstract class. Any dynamical system added to the toolbox need to extend this class and implement the method evaluate.
         * The class has been designed to allow external forces in the dynamics. These are defined through a Thrust vector or Control vector functions
         * given as an input in the constructor. It is possible to define for each dynamics also a scaling factor for the position and the time.
         * The planetary constants is updated accordingly [L^3/T^2] and the position scaled.
         */
	template < class T >
        class base_dynamics
        {

        public:

            /**
             * @brief base_dynamics constructor.
             *
             * In the constructor the name of the dynamics, the thrust and control profiles are initialized as class members.
             * @param name dynamical system name
             * @param r_scale position scaling factor. Default value = 1 (no scaling).
             * @param t_scale time scaling factor. Default value = 1 (no scaling)
             */
            base_dynamics(const std::string &name,
                          const double &r_scale = 1.0,
                          const double &t_scale = 1.0);
            virtual ~base_dynamics();

            /**
             * @brief Function to evaluate the dinamics at a given instant of time and a given state. It is a virtual function so any class that inherits from base_dynamics need to implement it
             * @param[in] t time
             * @param[in] state state values at time t
             * @param[out] dstate derivative of the states at time t
             * @return
             */
            virtual int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const = 0;
            
            /**
             * @brief Returns the name of the dynamics
             * @param[out] name name of the dynamics
             */
            std::string get_name() const;

        protected:
            std::string m_name;
            double m_r_scale;
            double m_t_scale;
        };

    }
}

#endif // SMARTMATH_BASE_DYNAMICS_H
