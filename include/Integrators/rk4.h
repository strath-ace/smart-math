#ifndef SMARTMATH_RK4_H
#define SMARTMATH_RK4_H

#include "base_rungekutta.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        /**
         * @brief The Runge Kutta 4 integrator scheme
         *
         * The class model the Runge Kutta fourth order integration scheme
         */
        template < class T >
        class rk4: public base_rungekutta<T>
        {

        private:
            using base_rungekutta<T>::m_name;
            using base_rungekutta<T>::m_dyn;

        public:

            using base_rungekutta<T>::integrate;

            /**
             * @brief rk4 constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             */
            rk4(const dynamics::base_dynamics<T> *dyn) : base_rungekutta<T>("Runge Kutta 4 fixed time-step", dyn){}

            /**
              * @brief ~rk4 deconstructor
              */
            ~rk4(){}

            /**
             * @brief performs one integration step from the RK4 method
             *
             * The method implements one step of the RK4 scheme to integrate with given initial time,
             * final time, initial state condition(constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] h time step
             * @param[in] x0 vector of initial states
             * @param[out] xfinal vector of final states
             * @return
             */
            int integration_step(const double &ti, const double &h, const std::vector<T> &x0, std::vector<T> &xfinal) const{

                std::vector<T> dx=x0;
                std::vector<T> x_temp=x0, k1=x0, k2=x0, k3=x0, k4=x0;
                unsigned int l = x0.size();
                double t1, t2, t3, t4;
                t1 = ti;
                t2 = t1 + h/2.0;
                t3 = t1 + h/2.0;
                t4 = t1 + h;

                //* Evaluate k1 = f(x).
                m_dyn->evaluate(t1, x0, k1);

                //* Evaluate k2 = f(x+h/2*k1),
                for(unsigned int j=0; j<l; j++)
                    x_temp[j] = x0[j]+k1[j]*h/2.0;
                m_dyn->evaluate(t2, x_temp, k2);

                //* Evaluate k3 = f(x+h/2*k2),
                for(unsigned int j=0; j<l; j++)
                    x_temp[j] = x0[j]+k2[j]*h/2.0;
                m_dyn->evaluate(t3, x_temp, k3);

                //* Evaluate k4 = f(x+h*k3),
                for(unsigned int j=0; j<l; j++)
                    x_temp[j] = x0[j]+k3[j]*h;
                m_dyn->evaluate(t4, x_temp, k4);

                //* Return x(t+h) computed from third-order Runge Kutta.
                xfinal=x0;
                for(unsigned int j=0; j<l; j++)
                    xfinal[j] +=  (k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j])*h/6.0;

                return 0;
            }

            /**
             * @brief integration_step, the same as above but for Eigen
             */
            int integration_step(const double &ti, const double &h, const Eigen::VectorXd &x0, Eigen::Ref<Eigen::VectorXd> xfinal) const{

                Eigen::VectorXd x_temp=x0, k1=x0, k2=x0, k3=x0, k4=x0;

                double t1, t2, t3, t4;
                t1 = ti;
                t2 = t1 + h/2.0;
                t3 = t1 + h/2.0;
                t4 = t1 + h;

                //* Evaluate k1 = f(x).
                m_dyn->evaluate(t1, x0, k1);

                //* Evaluate k2 = f(x+h/2*k1),
                x_temp = x0 + k1 *h/2.0;
                m_dyn->evaluate(t2, x_temp, k2);

                //* Evaluate k3 = f(x+h/2*k2),
                x_temp = x0 + k2 * h / 2.0;
                m_dyn->evaluate(t3, x_temp, k3);

                //* Evaluate k4 = f(x+h*k3),
                x_temp = x0 + k3 * h;
                m_dyn->evaluate(t4, x_temp, k4);

                //* Return x(t+h) computed from fourth-order Runge Kutta.
                xfinal = x0 + ( k1 + 2.0 * k2 + 2.0 * k3 + k4 ) * h / 6.0 ;

                return 0;
            }

        };

    }
}

#endif // SMARTMATH_RK4_H
