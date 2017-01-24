#ifndef SMARTMATH_AB_H
#define SMARTMATH_AB_H

#include "base_integrator.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        template < class T >
        class AB: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            int m_order;
            std::vector<double> m_gamma;

        public:
            /**
             * @brief Adam Bashforth constructor
             *
             * The integrator is initialized with the super class constructor. 
             * @param dyn
             * @param order chosen order for Adam Bashforth method
             */
            AB(const dynamics::base_dynamics<T> *dyn, const int order=6);

            /**
              * @brief ~AB deconstructor
              */
            ~AB();

            /**
             * @brief backward_differences computes backward differences  
             *
             * The method computes the backward differences for the Adam scheme
             * @param[in] f vector of saved state vectors for multistep scheme
             * @param[in] m number of saved steps
             * @param[out] Df vector of backward differences
             * @return
             */ 
            int backward_differences(const std::vector<std::vector<T> > &f, const int &m, std::vector<std::vector<T> > &Df) const;
    
            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps
             *
             * The method implements the Adam Bashforth scheme to integrate with given initial time,
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
             * The method implements the Adam Bashforth scheme to integrate with given initial time,
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
             * @brief integration-step performs one integration step from the Adam Bashforth method
             *
             * The method implements one step of the Adam Bashforth scheme 
             * @param[in] m number of saved steps
             * @param[in] h time step
             * @param[in] x0 vector of initial states
             * @param[in] Df vector of finite differences vectors 
             * @param[out] xfinal vector of final states
             * @return
             */	
            int integration_step(const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &Df, std::vector<T> &xfinal) const;
            
            /**
             * @brief initialize method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method initializes via RK4 the Adam Bashforth scheme for an integration with step-size h starting at given initial time and condition 
             * @param[in] m number of saved steps
             * @param[in] ti initial time instant
             * @param[in] h step size
             * @param[in] x0 vector of initial states
             * @param[out] f vector of saved state vectors for multistep scheme
             * @return
             */     
            int initialize(const int &m, const double &ti, const double &h, const std::vector<T> &x0, std::vector<std::vector<T> > &f) const;
    
        };

    }
}

#endif // SMARTMATH_AB_H
