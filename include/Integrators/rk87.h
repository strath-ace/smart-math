#ifndef SMARTMATH_RK87_H
#define SMARTMATH_RK87_H

#include <cmath>
#include "base_integrator.h"
#include "../exception.h"

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
            double m_tol;
            double m_multiplier;
            double m_minstep_events;

        public:
            /**
             * @brief rk87 constructor
             *
             * @param dyn
             * @param tolerance for error estimarion in step-size control
             * @param max multiplier for step-size control
             * @param min time step for events detection
             */
            rk87(const dynamics::base_dynamics<T> *dyn, const double tol=1.0e-7, const double multiplier=5.0, const double minstep_events=1.0e-4);

            /**
              * @brief ~rk87 deconstructor
              */
            ~rk87();

            /**
             * @brief integrate method to integrate bewteen two given time steps, initial condition and step lenght
             *
             * The method implements the RK8(7)-13 scheme to integrate with given initial time,
             * final time, initial state condition and number of steps (constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const;

            /**
             * @brief integrate method to integrate bewteen two given time steps, initial condition and step lenght
             *
             * The method implements the RK8(7)-13 scheme to integrate with given initial time,
             * final time, initial state condition and number of steps (constant stepsize)
             * @param[in] ti initial time instant
             * @param[out] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @param[in] event function
             * @return
             */
            int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal, std::vector<int> (*g)(std::vector<T> x, double d)) const;

            /**
             * @brief computes estimated error for step-size control
             *
             * 
             * @param[in] x current state
             * @param[out] val value of estimated error
             * @return
             */
            int error(const double &x, double &val) const;
#ifdef ENABLE_SMARTUQ
            int error(const smartuq::polynomial::chebyshev_polynomial<double> &x, double &val) const;
            int error(const smartuq::polynomial::taylor_polynomial<double> &x, double &val) const;
#endif

        };

    }
}

#endif // SMARTMATH_RK87_H
