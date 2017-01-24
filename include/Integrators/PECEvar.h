#ifndef SMARTMATH_PECEVAR_H
#define SMARTMATH_PECEVAR_H

#include "base_integrator.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        template < class T >
        class PECEvar: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            int m_order_max, m_order_min;
            std::vector<double> m_gamma;
            double m_tol, m_multiplier, m_minstep_events;

        public:
            /**
             * @brief Adam Bashforth Moulton constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             * @param order_max maximum order for predictor-corrector
             * @param order_min minimum order for predictor-corrector
             * @param tol tolerance for estimated integration error             
             */
            PECEvar(const dynamics::base_dynamics<T> *dyn, const int order_max=8, const int order_min=4, const double tol=1.0e-7, const double multiplier=5.0, const double minstep_events=1.0e-4);

            /**
              * @brief ~ABM deconstructor
              */
            ~PECEvar();

            /**
             * @brief computes backward differences 
             *
             * The method computes the backward differences for the Adam Moulton scheme
             * @param[in] f vector of saved state vectors for multistep scheme
             * @param[in] m order
             * @param[out] Df vector of backward differences
             * @return
             */ 
            int backward_differences(const std::vector<std::vector<T> > &f, const int &m, std::vector<std::vector<T> > &Df) const;
            
            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps
             *
             * The method implements the Adam Bashforth Moulton scheme to integrate with given initial time,
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
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method implements the Adam Bashforth Moulton 6 scheme to integrate with given initial time,
             * final time, initial state condition and number of steps (constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] x_history vector of intermediate state vector (including final one)
             * @param[out] t_history vector of intermediate times (including final one)
             * @return
             */
            int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history) const;

            /**
             * @brief integrate method to integrate bewteen two given time steps, with initial condition and initial guess for step-size while handling events
             *
             * The method implements the ABM scheme to integrate with given initial time,
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
            int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &xfinal, std::vector<double> t_history, std::vector<int> (*g)(std::vector<T> x, double d)) const;

            /**
             * @brief integrate method to integrate bewteen two given time steps, with initial condition and initial guess for step-size while handling events
             *
             * The method implements the ABM scheme to integrate with given initial time,
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
            int integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal, std::vector<int> (*g)(std::vector<T> x, double d)) const;

             /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps
             *
             * The method implements the correction step in the Adam Bashforth Moulton scheme 
             * @param[in] m order
             * @param[in] h step-size
             * @param[in] x0 vector of initial states
             * @param[in] Df vector of finite differences vectors            
             * @param[out] xfinal vector of final states
             * @return
             */
            int correction(const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &Df, std::vector<T> &xfinal) const;
            

            int integration_step(const double &ti, const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal, T &er) const;

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

#endif // SMARTMATH_PECEVAR_H
