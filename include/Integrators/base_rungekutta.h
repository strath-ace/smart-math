/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/

#ifndef SMARTMATH_BASE_RUNGEKUTTA_H
#define SMARTMATH_BASE_RUNGEKUTTA_H

#include "base_integrator.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The base_rungekutta class is a template abstract class. Any fixed-step Runge-Kutta algorithm added to the toolbox needs to inherit from it and implement the method integration_step()
         *
         * The base_rungekutta class is a template abstract class. Any fixed-step Runge-Kutta algorithm added to the toolbox needs to inherit from it and implement the method that performs on integration step between to given times given the initial state
         */
        template < class T >
        class base_rungekutta: public base_integrator<T>
        {

        protected:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            /**
             * @brief m_stages number of stages
             */            
            int m_stages;
            /**
             * @brief m_coeT coefficients for evaluation in time inside the integration step
             */             
            std::vector<double> m_coeT;
            /**
             * @brief m_coeK coefficients for evaluation in state inside the integration step (only valid for methods whose Butcher tableau has zeros everywhere except on the sub-diagonal)
             */             
            std::vector<double> m_coeK;
            /**
             * @brief m_coeX coefficients for state update inside the integration step
             */             
            std::vector<double> m_coeX;

        public:

            using base_integrator<T>::integrate;
            using base_integrator<T>::integrate_eigen;

            /**
             * @brief base_rungekutta constructor
             *
             * The constructor initializes the name of the integrator and a pointer to the dynamical system to be integrated
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             */
            base_rungekutta(const std::string &name, const dynamics::base_dynamics<T> *dyn) : base_integrator<T>(name, dyn){}

            /**
             * @brief ~base_rungekutta deconstructor
             */
            virtual ~base_rungekutta(){}

            /**
             * @brief integration_step performs one integration step from the Runge-Kutta scheme
             *
             * The method implements one step of a Runge-Kutta scheme to integrate with given initial time,
             * final time, initial state condition (constant stepsize)
             * This implementation only works for methods whose Butcher tableau has only a sub-diagonal of non-zero coefficients
             * It needs to be overwritten in other cases e.g. Kutta's third order method
             * @param[in] ti initial time instant
             * @param[in] h time step
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integration_step(const double &ti, const double &h, const std::vector<T> &x0, std::vector<T> &xfinal) const{

                std::vector<T> x_temp = x0, k = x0;
                unsigned int l = x0.size();

                xfinal = x0;

                for(int i = 0; i < m_stages; i++)
                {
                    x_temp = x0;
                    if(i == 0)
                    {
                        m_dyn->evaluate(ti + h * m_coeT[i], x_temp, k);
                        for(unsigned int j = 0; j < l; j++)
                            xfinal[j] += m_coeX[i] * k[j] * h;
                    }
                    else
                    { 
                        for(unsigned int j = 0; j < l; j++)
                            x_temp[j] += m_coeK[i - 1] * k[j] * h;    
                        m_dyn->evaluate(ti + h * m_coeT[i], x_temp, k);
                        for(unsigned int j = 0; j < l; j++)
                            xfinal[j] += m_coeX[i] * k[j] * h;        
                    }
                }              

                return 0;
            }

            /**
             * @brief integration_step_eigen performs one integration step from the Runge-Kutta scheme handling Eigen vectors
             *
             * The method implements one step of a Runge-Kutta scheme to integrate with given initial time,
             * final time, initial state condition (constant stepsize)
             * This implementation only works for methods whose Butcher tableau has only a sub-diagonal of non-zero coefficients
             * It needs to be overwritten in other cases e.g. Kutta's third order method
             * @param[in] ti initial time instant
             * @param[in] h time step
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integration_step_eigen(const double &ti, const double &h, const Eigen::VectorXd &x0, Eigen::Ref<Eigen::VectorXd> xfinal) const{

                Eigen::VectorXd x_temp=x0, k=x0;

                xfinal = x0;

                for(int i = 0; i < m_stages; i++)
                {
                    x_temp = x0;
                    if(i == 0)
                    {
                        m_dyn->evaluate_eigen(ti + h * m_coeT[i], x_temp, k);
                        xfinal += m_coeX[i] * k * h;
                    }
                    else
                    { 
                        x_temp += m_coeK[i - 1] * k * h;    
                        m_dyn->evaluate_eigen(ti + h * m_coeT[i], x_temp, k);
                        xfinal += m_coeX[i] * k * h;        
                    }
                }              

                return 0;
            }


            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method implements a fixed-step Runge-Kutta scheme to integrate with given initial time,
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

            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps (saving intermediate states) handling Eigen vectors
             *
             * The method implements a fixed-step Runge-Kutta scheme to integrate with given initial time,
             * final time, initial state condition and number of steps (constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] x_history vector of intermediate state vector (including final one)
             * @param[out] t_history vector of intermediate times (including final one)
             * @return
             */
            int integrate_eigen(const double &ti, const double &tend, const int &nsteps, const Eigen::VectorXd &x0, Eigen::Ref<Eigen::MatrixXd> x_history, Eigen::Ref<Eigen::VectorXd> t_history) const{

                t_history.setZero();
                x_history.setZero();

                Eigen::VectorXd x=x0, x_temp=x0;

                double t=ti, h = (tend-ti)/nsteps;

                for(int i=0; i<nsteps; i++){
                    integration_step_eigen(t,h,x,x_temp);
                    t+=h;
                    x=x_temp;
                    t_history(i) = t;
                    x_history.col(i) = x;
                }

                return 0;
            }

        };

    }
}

#endif // SMARTMATH_BASE_RUNGEKUTTA_H
