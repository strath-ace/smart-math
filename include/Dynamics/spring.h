/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2017 University of Strathclyde and Authors-------
-------- e-mail: romain.serra@strath.ac.uk ---------------------------
--------- Author: Romain Serra ---------------------------------------
*/

#ifndef SMARTMATH_SPRING_H
#define SMARTMATH_SPRING_H

#include "hamiltonian_mixedvar.h"
#include "../exception.h"

namespace smartmath
{
    namespace dynamics {

        /**
         * @brief The harmonic oscillator as a Hamiltonian system with mixed variables
         *
         * The problem models the Hamiltonian dynamics of a 1-D spring with unit scales i.e.
         * \f{eqnarray*}{
            \ddot{x} &=& -x
          \f}  
         * whose analytical solution is \f$x (t)= \sqrt{x(t_0)^2+\dot{x}(t_0)^2}\sin(t-t_0+\arctan2(x(t_0),\dot{x}(t_0))) \f$.
         * In order to introduce mixed variables for modified symplectic integration, the Hamiltonian is artificially written as \f$ H = H_0(\theta,\Theta) + V(q,p)\f$ with \f$ H_0 = \Theta \f$ and \f$ V = 0 \f$ where:          
         * \f{eqnarray*}{
            q &=& x \\
            p &=& \dot{x} \\
            \Theta &=& \frac{p^2}{2} + \frac{q^2}{2} \\
            \theta &= &\arctan2(q,p) 
          \f} 
         *    
         */
        template < class T >
        class spring: public hamiltonian_mixedvar<T>
        {

        private:
            using hamiltonian_mixedvar<T>::m_dim;

        public:

            /**
             * @brief spring constructor
             *
             * The constructor initializes the problem
             */
            spring() : hamiltonian_mixedvar<double>("Mathematical spring problem", 1, true)
            {

            }

            /**
              * @brief ~spring deconstructor
              */
            ~spring(){}

            /**
             * @brief evaluate differential equations of the implemented Hamiltonian system
             * @param[in] t time in scaled units
             * @param[in] state vector in scaled units
             * @param[out] dstate derivative in scaled units
             * @return exit flag (0=success)
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const{

                /* sanity checks */
                if(state.size() != 2 * m_dim)
                    smartmath_throw("EVALUATE: the Hamiltonian state must have a consistent dimension");

                dstate[0] = state[1];
                dstate[1] = -state[0];

                return 0;
            }

            /**
             * @brief DHq evaluates partial derivative of Hamiltonian with respect to q
             * @param[in] t time in TU
             * @param[in] q vector of primary canonical position
             * @param[in] p vector of primary canonical momenta 
             * @param[out] dH partial derivative of H with respect to q in scaled units 
             * @return exit flag (0=success)
             */
            int DHq(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const{

                dH[0] = 0.0;

                return 0;
            };

            /**
             * @brief DHq evaluates partial derivative of Hamiltonian with respect to p
             * @param[in] t time in TU
             * @param[in] q vector of primary canonical position
             * @param[in] p vector of primary canonical momenta 
             * @param[out] dH partial derivative of H with respect to p in scaled units 
             * @return exit flag (0=success)
             */
            int DHp(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const{

                dH[0] = p[0];

                return 0;
            };            

            /**
             * @brief DHq2 computes the partial derivative of the Hamiltonian with respect to the second 'position' q2
             *
             * The method computes the partial derivative of the Hamiltonian with respect to q2
             * @param[in] t time in scaled units
             * @param[in] q2 vector in scaled units
             * @param[in] p2 vector in scaled units
             * @param[out] dH2 vector of partial derivatives of H w.r.t. the vector q2
             * @return exit flag (0=success)
             */
            int DHq2(const double &t, const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &dH2) const{

                dH2[0] = 0.0;

                return 0;
            };

            /**
             * @brief DHp2 computes the partial derivative of the Hamiltonian with respect to the second 'momentum' p2
             *
             * The method computes the partial derivative of the Hamiltonian with respect to p2
             * @param[in] t time in scaled units
             * @param[in] q2 vector in scaled units
             * @param[in] p2 vector in scaled units
             * @param[out] dH2 vector of partial derivatives of H w.r.t. the vector p2
             * @return exit flag (0=success)
             */
            int DHp2(const double &t, const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &dH2) const{

                dH2[0] = 1.0;

                return 0;
            };

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
            int conversion(const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &q2, std::vector<T> &p2) const{

                p2[0] = (q[0] * q[0] + p[0] * p[0]) / 2.0;
                q2[0] = atan2(q[0], p[0]);

                return 0;
            }; 

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
            int conversion2(const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &q, std::vector<T> &p) const{

                T inter = sqrt(2.0 * p2[0]);
                p[0] = inter * cos(q2[0]);
                q[0] = inter * sin(q2[0]);

                return 0;
            }; 

        };

    }
}

#endif // SMARTMATH_SPRING_H
