/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2017 University of Strathclyde and Authors-------
-------- e-mail: romain.serra@strath.ac.uk ---------------------------
--------- Author: Romain Serra ---------------------------------------
*/


#ifndef SMARTMATH_HAMILTONIAN_MIXEDVAR_H
#define SMARTMATH_HAMILTONIAN_MIXEDVAR_H

#include "base_hamiltonian.h"
#include "../exception.h"

namespace smartmath
{
    namespace dynamics {

        /**
         * @brief The hamiltonian_mixedvar class is a template abstract class. Any Hamiltonian system with mixed variables needs to inherit from it
         *
         * The hamiltonian_mixedvar class is a template abstract class. Any Hamiltonian system added to the toolbox needs to inherit from it
         * The system has two sets of canonical variables and the Hamiltonian writes H(q, p) =  H(q2, p2).
         */
        template < class T >
        class hamiltonian_mixedvar: public base_hamiltonian<T>
        {

        protected:
        	using smartmath::dynamics::base_hamiltonian<T>::m_dim;

        public:
            
        	 /**
             * @brief hamiltonian_mixedvar constructor
             *
             * The constructor initializes the name of the Hamiltonian dynamics, its half-dimension and a flag about its separability
             * @param name integrator name
             * @param dim half-order of the Hamiltonian system
             * @param separable boolean precising whether the system is separable or not
             */
            hamiltonian_mixedvar(const std::string &name, const unsigned int &dim, const bool &separable = false): base_hamiltonian<T>(name, dim, separable){}

            /**
             * @brief ~hamiltonian_mixedvar deconstructor
             */
            ~hamiltonian_mixedvar(){}

            /**
             * @brief evaluate differential equations of the implemented Hamiltonian system
             * @param[in] t time in scaled units
             * @param[in] state vector in scaled units
             * @param[out] dstate vector of time derivatives in scaled units
             * @return exit flag (0=success)
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const{

                smartmath_throw("EVALUATE: do not use EVALUATE for Hamiltonian systems with mixed variables because part of the dynamics is missing");

                return 0;
            }          

            /**
             * @brief DHq2 computes the partial derivative of the Hamiltonian with respect to the second 'position' q2
             *
             * The method computes the partial derivative of the Hamiltonian with respect to q2
             * @param[in] t time in scaled units
             * @param[in] q2 vector in scaled units
             * @param[in] p2 vector in scaled units
             * @param[out] dH2 vector of partial derivatives of H w.r.t. the vector q
             * @return exit flag (0=success)
             */
            virtual int DHq2(const double &t, const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &dH2) const = 0;

            /**
             * @brief DHp2 computes the partial derivative of the Hamiltonian with respect to the second 'momentum' p2
             *
             * The method computes the partial derivative of the Hamiltonian with respect to p2
             * @param[in] t time in scaled units
             * @param[in] q2 vector in scaled units
             * @param[in] p2 vector in scaled units
             * @param[out] dH2 vector of partial derivatives of H w.r.t. the vector p
             * @return exit flag (0=success)
             */
            virtual int DHp2(const double &t, const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &dH2) const = 0;            

            /**
             * @brief conversion performs the conversion from the first set of canonical variables to the second one
             *
             * The method converts the first set of canonical variables into the second one
             * @param[in] q vector in scaled units
             * @param[in] p vector in scaled units            
             * @param[out] q2 vector in scaled units
             * @param[out] p2 vector in scaled units
             * @return exit flag (0=success)
             */
            virtual int conversion(const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &q2, std::vector<T> &p2) const = 0; 

            /**
             * @brief conversion performs the conversion from the second set of canonical variables to the first one
             *
             * The method converts the second set of canonical variables into the first one
             * @param[in] q2 vector in scaled units
             * @param[in] p2 vector in scaled units            
             * @param[out] q vector in scaled units
             * @param[out] p vector in scaled units
             * @return exit flag (0=success)
             */
            virtual int conversion2(const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &q, std::vector<T> &p) const = 0; 

        };

    }
}

#endif // SMARTMATH_HAMILTONIAN_MIXEDVAR_H
