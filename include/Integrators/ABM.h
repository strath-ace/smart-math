#ifndef SMARTMATH_ABM_H
#define SMARTMATH_ABM_H

#include "base_integrator.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        template < class T >
        class ABM: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            int m_order;
            std::vector<double> m_gamma;

        public:
            /**
             * @brief Adam Bashforth Moulton constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             * @param order order for Adam Bashforth Moulton scheme        
             */
            ABM(const dynamics::base_dynamics<T> *dyn, const int order=6);

            /**
              * @brief ~ABM deconstructor
              */
            ~ABM();

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
             * @brief correction method to integrate between two given time steps, initial condition and number of steps
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
            
            /**
             * @brief integration_step method to perform one step of integration
             *
             * The method implements one step of the Adam Bashforth Moulton scheme 
             * @param[in] ti initial time instant
             * @param[in] m order
             * @param[in] h step-size
             * @param[in] x0 vector of initial states
             * @param[in] f vector of saved state vectors for multistep scheme             
             * @param[out] xfinal vector of final states 
             * @return
             */
            int integration_step(const double &ti, const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal) const;
             	
        };

    }
}

#endif // SMARTMATH_ABM_H
