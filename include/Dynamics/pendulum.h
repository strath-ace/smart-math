/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde--------------
-------------------- Author: Romain Serra -------------------------
-------------- e-mail: romain.serra@strath.ac.uk ------------------
*/


#ifndef SMARTMATH_PENDULUM_H
#define SMARTMATH_PENDULUM_H

#include "hamiltonian_momentum.h"
#include "../exception.h"

namespace smartmath
{
    namespace dynamics {

        /**
         * @brief The mathematical pendulum problem
         *
         * The problem models the dynamics of a pendulum with unit scales.
         *
         */
        template < class T >
        class pendulum: public hamiltonian_momentum<T>
        {

        private:
            using hamiltonian_momentum<T>::m_dim;

        public:

            /**
             * @brief pendulum constructor
             *
             * The constructor initializes the problem
             */
            pendulum() : hamiltonian_momentum<double>("Mathematical pendulum problem", 1, true)
            {

            }

            /**
              * @brief ~pendulum deconstructor
              */
            ~pendulum(){}

            /**
             * @brief DHq computes the partial derivative of the Hamiltonian with respect to the position
             *
             * The method computes the partial derivative of the Hamiltonian with respect to q the vector of position
             * @param[in] t time in scaled units
             * @param[in] q position vector in scaled units
             * @param[in] p momenta vector in scaled units
             * @param[out] dH vector of partial derivatives of H w.r.t. the vector q
             * @return exit flag (0=success)
             */
            int DHq(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const{

                if(q.size() != m_dim)
                    smartmath_exception("DHQ: the position must have the correct dimension");
                if(dH.size() != m_dim)
                    smartmath_exception("DHQ: the derivative w.r.t. position must have the correct dimension");                

                dH[0] = sin(q[0]);

                return 0;
            };
        };

    }
}

#endif // SMARTMATH_PENDULUM_H
