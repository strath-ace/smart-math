#ifndef SMARTMATH_RKF45_H
#define SMARTMATH_RKF45_H

#include "base_integrator.h"
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
        class rkf45: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            double m_tol;
            double m_multiplier;
            double m_minstep_events;

        public:
            /**
             * @brief rkf45 constructor
             *
             * @param dyn
             * @param tolerance for error estimation in step-size control
             * @param max multiplier for step-size control
             * @param min time step for events detection
             */
            rkf45(const dynamics::base_dynamics<T> *dyn, const double tol=1.0e-7, const double multiplier=5.0, const double minstep_events=1.0e-4);

            /**
              * @brief ~rkf45 deconstructor
              */
            ~rkf45();

            /**
             * @brief integrate method to integrate bewteen two given time steps, initial condition and initial guess for step-size
             *
             * The method implements the RKF4(5) scheme to perform one integration step
             * @param[in] ti initial time
             * @param[in] h step size
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @param[out] er estimated error 
             * @return
             */
            int integration_step(const double &ti, const double &h, const std::vector<T> &x0, std::vector<T> &xfinal, T &er) const;

            /**
             * @brief integrate method to integrate bewteen two given time steps, with initial condition and initial guess for step-size
             *
             * The method implements the RKF4(5) scheme to integrate with given initial time,
             * final time, initial state condition and initial guess for step-size
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps initial guess for number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &xfinal, std::vector<double> &t_history) const;

            /**
             * @brief integrate method to integrate bewteen two given time steps, with initial condition and initial guess for step-size
             *
             * The method implements the RKF4(5) scheme to integrate with given initial time,
             * final time, initial state condition and initial guess for step-size
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps initial guess for number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of intemrediate states
             * @param[out] t_history vector of intermediate times
             * @return
             */
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const;


            /**
             * @brief integrate method to integrate bewteen two given time steps, with initial condition and initial guess for step-size while handling events
             *
             * The method implements the RKF4(5) scheme to integrate with given initial time,
             * final time, initial state condition and initial guess for step-size
             * @param[in] ti initial time instant
             * @param[out] tend final time instant
             * @param[in] nsteps initial guess for number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @param[in] event function
             * @return
             */
            int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal, std::vector<int> (*g)(std::vector<T> x, double d)) const;
            
            /**
             * @brief integrate method to integrate bewteen two given time steps, with initial condition and initial guess for step-size while handling events
             *
             * The method implements the RKF4(5) scheme to integrate with given initial time,
             * final time, initial state condition and initial guess for step-size
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps initial guess for number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of intemrediate states
             * @param[out] t_history vector of intermediate times
             * @param[in] event function             
             * @return
             */
            int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> t_history, std::vector<int> (*g)(std::vector<T> x, double d)) const;

            static std::vector<int> dummy_event(std::vector<T> x, double d);

            /**
             * @brief returns a double equal to the input for real numbers and something meaningful for polynomials
             *
             * 
             * @param[in] x estimated error
             * @param[out] val double equal to x for real numbers and something else for polynomials
             * @return
             */
            int error(const float &x, double &val) const;
            int error(const double &x, double &val) const;
            int error(const long double &x, double &val) const;
            #ifdef ENABLE_SMARTUQ
                int error(const smartuq::polynomial::chebyshev_polynomial<double> &x, double &val) const;
                int error(const smartuq::polynomial::chebyshev_polynomial<float> &x, double &val) const;
                int error(const smartuq::polynomial::chebyshev_polynomial<long double> &x, double &val) const;                        
                int error(const smartuq::polynomial::taylor_polynomial<double> &x, double &val) const;
                int error(const smartuq::polynomial::taylor_polynomial<float> &x, double &val) const;
                int error(const smartuq::polynomial::taylor_polynomial<long double> &x, double &val) const;            
            #endif

        };

    }
}

#endif // SMARTMATH_RKF45_H
