/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
-------------- e-mail: romain.serra@strath.ac.uk ---------------------
--------- Author: Annalisa Riccardi, Carlos Ortega and Romain Serra --
*/


#ifndef SMARTMATH_BASE_INTEGRATOR_H
#define SMARTMATH_BASE_INTEGRATOR_H

#include "../LinearAlgebra/Eigen/Core"
#include "../Dynamics/base_dynamics.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The base_integrator class is a template abstract class. Any integrator added to the toolbox needs to inherit from it and implement the method integrate()
         *
         * The base_integrator class is a template abstract class. Any integrator added to the toolbox needs to inherit from it and implement the method that integrates between two given times, with initial state and stepsize
         */
        template < class T >
        class base_integrator
        {

        public:

            /**
             * @brief base_integrator constructor
             *
             * The constructor initialize the name of the integrator and a pointer to the dynamical system to eb integrated
             * @param name integrator name
             * @param dyn pointer to a base_dynamics object
             */
            base_integrator(const std::string &name, const dynamics::base_dynamics<T> *dyn):m_name(name),m_dyn(dyn){}

            /**
             * @brief ~base_integrator deconstructor
             */
            virtual ~base_integrator(){}

            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method implements the scheme to integrate with given initial time,
             * final time, initial state condition and number of steps (constant stepsize) returning the full history of propagation
             * @param[in] ti initial time instant
             * @param[in] tend final time instant
             * @param[in] nsteps number of integration steps
             * @param[in] x0 vector of initial states
             * @param[out] x_history vector of intermediate state vector (including final one)
             * @param[out] t_history vector of intermediate times (including final one)
             * @return
             */
            virtual int integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history)  const = 0;

            virtual int integrate_eigen(const double &ti, const double &tend, const int &nsteps, const Eigen::VectorXd &x0, Eigen::Ref<Eigen::MatrixXd> x_history, Eigen::Ref<Eigen::VectorXd> t_history)  const
            { smartmath_throw("integrate_function using Eigen not implemented "); return 1; }

            /**
             * @brief integrate method to integrate from initial conditions to a final time with a given number of steps
             *
             * The method implements the corresponding integration scheme with given initial time,
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

            int integrate_eigen(const double &ti, const double &tend, const int &nsteps, const Eigen::VectorXd &x0, Eigen::Ref<Eigen::VectorXd> xfinal) const{

                Eigen::MatrixXd x_history = Eigen::MatrixXd::Zero(x0.size(),nsteps);
                Eigen::VectorXd t_history = Eigen::VectorXd::Zero(nsteps);

                integrate_eigen(ti,tend,nsteps,x0,x_history,t_history);

                xfinal=x_history.rightCols(1);

                return 0;
            }

            /**
             * @brief get_name return integrator name
             *
             * Function to get the name of the integration scheme
             * @return
             */
            std::string get_name() const {return m_name;}

            /**
             * @brief set_comments changes the flag for printing comments
             *
             * The method works as a toggle command for printing comments.
             * @param[in] status new status for comments: true = ON, false = OFF
             */
            void set_comments(const bool &status)
            {
                m_comments = status;
            }


        protected:
            /**
             * @brief m_name integrator name
             */
            std::string m_name;
            /**
             * @brief m_dyn pointer to a base_dynamics system to be integrated
             */
            const dynamics::base_dynamics<T> *m_dyn;
            /**
             * @brief m_comments status for comment printing
             */
            bool m_comments = true;
        };
    }
}

#endif // SMARTMATH_BASE_INTEGRATOR_H
