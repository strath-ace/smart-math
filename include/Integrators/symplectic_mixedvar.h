/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: romain.serra@strath.ac.uk -----------------------
----------------- Author: Romain Serra -------------------------------
*/



#ifndef SMARTMATH_SYMPLECTIC_MIXEDVAR_H
#define SMARTMATH_SYMPLECTIC_MIXEDVAR_H

#include "../LinearAlgebra/Eigen/Core"
#include "../Dynamics/hamiltonian_mixedvar.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The base_integrator class is a template abstract class. Any integrator added to the toolbox needs to inherit from it and implement the method integrate()
         *
         * The base_integrator class is a template abstract class. Any integrator added to the toolbox needs to inherit from it and implement the method that integrates between to given times, initial state and stepsize
         */
        template < class T >
        class symplectic_mixedvar
        {

        public:

            /**
             * @brief base_integrator constructor
             *
             * The constructor initialize the name of the integrator and a pointer to the dynamical system to eb integrated
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             */
            symplectic_mixedvar(const std::string &name, const dynamics::hamiltonian_mixedvar<T> *dyn, const int &order): m_name(name), m_dyn(dyn), m_order(order){}

            /**
             * @brief ~base_integrator deconstructor
             */
            virtual ~symplectic_mixedvar(){}

            /**
             * @brief integration_step performs one integration step from the Runge-Kutta scheme
             *
             * The method implements one step of a Runge-Kutta scheme to integrate with given initial time,
             * final time, initial state condition(constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] h time step
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            virtual int integration_step(const double &ti, const double &tau, const std::vector<T> &x0, std::vector<T> &xfinal) const = 0;

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

                std::vector<T> dx = x0, x = x0, x_temp = x0;

                double t = ti, h = (tend-ti) / nsteps;

                for(int i = 0; i < nsteps; i++){
                    integration_step(t, h, x, x_temp);
                    t += h;
                    x = x_temp;
                    t_history.push_back(t);
                    x_history.push_back(x);
                }

                return 0;
            }

            /**
             * @brief integrate method to integrate bewteen two given time steps, initial condition and step lenght
             *
             * The virtual method is inherited by any subclass. It implements the corresponding integration scheme with given initial time,
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

        protected:

            std::string m_name;
            const dynamics::hamiltonian_mixedvar<T> *m_dyn;
            int m_order;
            std::vector<T> m_c, m_d;

        };
    }
}

#endif // SMARTMATH_SYMPLECTIC_MIXEDVAR_H
