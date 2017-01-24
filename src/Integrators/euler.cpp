#include "../../include/Integrators/euler.h"


using namespace smartmath;
using namespace smartmath::integrator;

template < class T >
euler<T>::euler(const dynamics::base_dynamics<T> *dyn) : base_integrator<T>("Explicit Euler integration scheme", dyn)
{
}

template < class T >
euler<T>::~euler(){

}

template < class T >
int euler<T>::integration_step(const double &t, const double &h, const std::vector<T> &x0, std::vector<T> &xfinal) const{
	
	std::vector<T> dx=x0;
	unsigned int l = x0.size();

	m_dyn->evaluate(t, x0, dx);

	xfinal=x0;	
	for(unsigned int j=0; j<l; j++){
		xfinal[j] += h*dx[j];
	}

	return 0;
}

template < class T >
int euler<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history) const{
	
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
int euler<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

	std::vector<std::vector<T> > x_history;
	std::vector<double> t_history;

	integrate(ti,tend,nsteps,x0,x_history,t_history);

	xfinal=x_history.back();

	return 0;
}


template class euler<double>;
template class euler<float>;
template class euler<long double>;
#ifdef ENABLE_SMARTUQ
template class euler<smartuq::polynomial::chebyshev_polynomial<double> >;
template class euler<smartuq::polynomial::chebyshev_polynomial<float> >;
template class euler<smartuq::polynomial::chebyshev_polynomial<long double> >;
template class euler<smartuq::polynomial::taylor_polynomial<double> >;
template class euler<smartuq::polynomial::taylor_polynomial<float> >;
template class euler<smartuq::polynomial::taylor_polynomial<long double> >;
#endif
