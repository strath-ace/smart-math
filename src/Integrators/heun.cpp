#include "../../include/Integrators/heun.h"


using namespace smartmath;
using namespace smartmath::integrator;

template < class T >
heun<T>::heun(const dynamics::base_dynamics<T> *dyn) : base_integrator<T>("Heun's method of order 2", dyn)
{
}

template < class T >
heun<T>::~heun(){

}

template < class T >
int heun<T>::integration_step(const double &t, const double &h, const std::vector<T> &x0, std::vector<T> &xfinal) const{
	
	std::vector<T> dx=x0, dx_temp=x0, x_temp=x0;
	unsigned int l = x0.size();

	m_dyn->evaluate(t, x0, dx);

	for(unsigned int j=0; j<l; j++){
	    x_temp[j] = x0[j] + h*dx[j];
	}
	m_dyn->evaluate(t+h, x_temp, dx_temp);
	
	xfinal=x0;
	for(unsigned int j=0; j<l; j++){
	    xfinal[j] += h/2.0*(dx[j] + dx_temp[j]);
	}

	return 0;
}

template < class T >
int heun<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history) const{
	
	x_history.clear();
	t_history.clear();

	std::vector<T> x(x0), x_temp(x0);

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
int heun<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

	std::vector<std::vector<T> > x_history;
	std::vector<double> t_history;

	integrate(ti,tend,nsteps,x0,x_history,t_history);

	xfinal=x_history.back();

	return 0;
}


template class heun<double>;
template class heun<float>;
template class heun<long double>;
#ifdef ENABLE_SMARTUQ
template class heun<smartuq::polynomial::chebyshev_polynomial<double> >;
template class heun<smartuq::polynomial::chebyshev_polynomial<float> >;
template class heun<smartuq::polynomial::chebyshev_polynomial<long double> >;
template class heun<smartuq::polynomial::taylor_polynomial<double> >;
template class heun<smartuq::polynomial::taylor_polynomial<float> >;
template class heun<smartuq::polynomial::taylor_polynomial<long double> >;
#endif


