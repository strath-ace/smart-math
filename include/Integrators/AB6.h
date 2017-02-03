#ifndef SMARTMATH_AB6_H
#define SMARTMATH_AB6_H

#include "base_integrator.h"
#include "../exception.h"

namespace smartmath
{
    namespace integrator {

        template < class T >
        class AB6: public base_integrator<T>
        {

        private:
            using base_integrator<T>::m_name;
            using base_integrator<T>::m_dyn;
            int m_order;
            std::vector<double> m_beta;

        public:
            /**
             * @brief Adam Bashforth 6 constructor
             *
             * The integrator is initialized with the super class constructor. No additional parameters are set.
             * @param dyn
             */
            AB6(const dynamics::base_dynamics<T> *dyn): base_integrator<T>("Adam Bashforth 6 integration scheme", dyn)
            {
	            m_order=6;

	            double prebeta[6]={-475.0,2877.0,-7298.0,9982.0,-7923.0,4277.0};
	            for(int i=0; i<m_order; i++){
		            m_beta.push_back(prebeta[i]/1440.0);
	            }

            }

            /**
              * @brief ~AB6 deconstructor
              */
            ~AB6(){}

            /**
             * @brief integration_step performs one integration step from the Adam Bashforth (order 6) method
             *
             * The method implements one step of the Adam Bashforth 6 scheme to integrate with given initial time,
             * final time, initial state condition(constant stepsize)
             * @param[in] ti initial time instant
             * @param[in] h time step
             * @param[in] x0 vector of initial states
             * @param[in] f vector of saved state vectors for multistep scheme             
             * @param[out] xfinal vector of final states
             * @return
             */
            int integration_step(const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal) const{

	            xfinal=x0;
	            for(int i=0; i<x0.size(); i++){
		            for(int j=0; j<m_order; j++){
			            xfinal[i]+=h*m_beta[j]*f[j][i];
	                }
	            }

	            return 0;
            }
            
            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps
             *
             * The method implements the Adam Bashforth 6 scheme to integrate with given initial time,
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
	       
            /**
             * @brief integrate method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method implements the Adam Bashforth 6 scheme to integrate with given initial time,
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

	            std::vector<T> x(x0), xp(x0), dx(x0);
	            std::vector<std::vector<T> > f, fp;

	            double t=ti, h = (tend-ti)/nsteps;

	            initialize(ti,h,x0,f);

                for(int k=0; k<nsteps; k++){
                	integration_step(h,x,f,xp);
                	x=xp;
                	t+=h;

                	/* Updating saved steps */
                	fp=f;
		            for(int j=0; j<m_order-1; j++){
			            f[j]=fp[j+1];
	                }
	                m_dyn->evaluate(t, x, dx);
	                f[m_order-1]=dx;
	                
	                /* Saving states */
	                t_history.push_back(t);
	                x_history.push_back(x);
	            }

	            return 0;
            }
    
            /**
             * @brief initialize method to integrate between two given time steps, initial condition and number of steps (saving intermediate states)
             *
             * The method initializes via RK4 the Adam Bashforth 6 scheme for an integration with step-size h starting at given initial time and condition 
             * @param[in] ti initial time instant
             * @param[in] h step size
             * @param[in] x0 vector of initial states
             * @param[out] f vector of saved state vectors for multistep scheme
             * @return
             */    
            int initialize(const double &ti, const double &h, const std::vector<T> &x0, std::vector<std::vector<T> > &f) const{

	            f.clear();

	            std::vector<T> dx(x0), x(x0), xp(x0);
	            std::vector< std::vector<T> > fp;

	            integrator::rk4<T> RK(m_dyn); // Runge kutta schemed used for initialization (here RK4)

	            /* Computing the initial saved steps */
	            m_dyn->evaluate(ti,x,dx);
	            fp.push_back(dx);
	            double t=ti;
	            for(int j=0; j<m_order-1; j++){
		            RK.integration_step(t,-h,x,xp);
		            t-=h;
		            x=xp;
		            m_dyn->evaluate(t,x,dx);
		            fp.push_back(dx);
	            }

	            /* Putting the saved steps in the right order */
	            for(int j=0; j<m_order; j++){
		            f.push_back(fp[m_order-j-1]);
	            }

	            return 0;
            }
    
        };

    }
}

#endif // SMARTMATH_AB6_H
