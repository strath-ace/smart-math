/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: romain.serra@strath.ac.uk -----------------------
----------------- Author: Romain Serra -------------------------------
*/



#ifndef SMARTMATH_BASE_SYMPLECTIC_H
#define SMARTMATH_BASE_SYMPLECTIC_H

#include "../LinearAlgebra/Eigen/Core"
#include "../Dynamics/base_hamiltonian.h"
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
        class base_symplectic
        {

        public:

            /**
             * @brief base_integrator constructor
             *
             * The constructor initialize the name of the integrator and a pointer to the dynamical system to eb integrated
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             */
            base_symplectic(const std::string &name, const dynamics::base_hamiltonian<T> *dyn):m_name(name),m_dyn(dyn){}

            /**
             * @brief ~base_integrator deconstructor
             */
            virtual ~base_symplectic(){}

            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method implements the scheme to integrate with given initial time,
             * final time, initial state condition and number of steps (constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] x_history vector of intermediate state vector (including final one)
             * @param[out] t_history vector of intermediate times (including final one)
             * @return
             */
            virtual int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history)  const = 0;

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
            const dynamics::base_hamiltonian<T> *m_dyn;

        };
    }
}

#endif // SMARTMATH_BASE_SYMPLECTIC_H
