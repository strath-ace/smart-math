#include "../../include/Integrators/rk4.h"


using namespace smartmath;
using namespace smartmath::integrator;

template < class T >
rk4<T>::rk4(const dynamics::base_dynamics<T> *dyn) : base_integrator<T>("Runge Kutta 4 fixed step time", dyn)
{
}

template < class T >
rk4<T>::~rk4(){

}

template < class T >
int rk4<T>::integration_step(const double &t, const double &h, const std::vector<T> &x0, std::vector<T> &xfinal) const{
	
	std::vector<T> dx=x0;
	std::vector<T> x_temp=x0, k1=x0, k2=x0, k3=x0, k4=x0;
	unsigned int l = x0.size();
	double t1, t2, t3, t4;
	t1 = t;
	t2 = t + h/2.0;
	t3 = t + h/2.0;
	t4 = t + h;

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


template < class T >
int rk4<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history) const{
	
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

template < class T >
int rk4<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

	std::vector<std::vector<T> > x_history;
	std::vector<double> t_history;

	integrate(ti,tend,nsteps,x0,x_history,t_history);

	xfinal=x_history.back();

	return 0;
}


template class rk4<double>;
template class rk4<float>;
template class rk4<long double>;
//template class rk4<astrocore::base_state*>;
#ifdef ENABLE_SMARTUQ
template class rk4<smartuq::polynomial::chebyshev_polynomial<double> >;
template class rk4<smartuq::polynomial::chebyshev_polynomial<float> >;
template class rk4<smartuq::polynomial::chebyshev_polynomial<long double> >;
template class rk4<smartuq::polynomial::taylor_polynomial<double> >;
template class rk4<smartuq::polynomial::taylor_polynomial<float> >;
template class rk4<smartuq::polynomial::taylor_polynomial<long double> >;
#endif
