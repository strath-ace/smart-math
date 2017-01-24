#include "../../include/Integrators/AB6.h"
#include "../../include/Integrators/rk4.h"


using namespace smartmath;
using namespace smartmath::integrator;

template < class T >
AB6<T>::AB6(const dynamics::base_dynamics<T> *dyn) : base_integrator<T>("Adam Bashforth 6 integration scheme", dyn)
{
	m_order=6;

	double prebeta[6]={-475.0,2877.0,-7298.0,9982.0,-7923.0,4277.0};
	for(int i=0; i<m_order; i++){
		m_beta.push_back(prebeta[i]/1440.0);
	}

}

template < class T >
AB6<T>::~AB6(){

}


template < class T >
int AB6<T>::integration_step(const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal) const{

	xfinal=x0;
	for(int i=0; i<x0.size(); i++){
		for(int j=0; j<m_order; j++){
			xfinal[i]+=h*m_beta[j]*f[j][i];
	    }
	}

	return 0;
}

template < class T >
int AB6<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history) const{

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

template < class T >
int AB6<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

	std::vector<std::vector<T> > x_history;
	std::vector<double> t_history;

	integrate(ti,tend,nsteps,x0,x_history,t_history);

	xfinal=x_history.back();

	return 0;
}

template < class T >
int AB6<T>::initialize(const double &ti, const double &h, const std::vector<T> &x0, std::vector<std::vector<T> > &f) const{

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


template class AB6<double>;
template class AB6<float>;
template class AB6<long double>;
#ifdef ENABLE_SMARTUQ
template class AB6<smartuq::polynomial::chebyshev_polynomial<double> >;
template class AB6<smartuq::polynomial::chebyshev_polynomial<float> >;
template class AB6<smartuq::polynomial::chebyshev_polynomial<long double> >;
template class AB6<smartuq::polynomial::taylor_polynomial<double> >;
template class AB6<smartuq::polynomial::taylor_polynomial<float> >;
template class AB6<smartuq::polynomial::taylor_polynomial<long double> >;
#endif
