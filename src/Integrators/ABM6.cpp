#include "../../include/Integrators/ABM6.h"


using namespace smartmath;
using namespace smartmath::integrator;

template < class T >
ABM6<T>::ABM6(const dynamics::base_dynamics<T> *dyn) : base_integrator<T>("Adam Bashforth Moulton algorithm with order 6", dyn)
{
	m_order=6;
	double prebeta[6]={27.0,-173.0,482.0,-798.0,1427.0,475.0};
	for(int i=0; i<m_order; i++)
		m_beta.push_back(prebeta[i]/1440.0);

}

template < class T >
ABM6<T>::~ABM6(){

}

template < class T >
int ABM6<T>::correction(const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal) const{

	xfinal=x0;
	for(int i=0; i<x0.size(); i++){
		for(int j=0; j<m_order; j++){
			xfinal[i]+=h*m_beta[j]*f[j][i];
	    }
	}

	return 0;
}

template < class T >
int ABM6<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history) const{
	
	t_history.clear();
	x_history.clear();

	std::vector<T> x(x0),xp(x0),dx(x0);
	std::vector<std::vector<T> > f, fp;

	double t=ti, h = (tend-ti)/nsteps;

	integrator::AB6<T> predictor(m_dyn);

	predictor.initialize(ti,h,x0,f); // initializing predictor-corrector by initializing predictor

    for(int k=0; k<nsteps; k++){

    	predictor.integration_step(h,x,f,xp); // prediction 

    	/* Updating the saved steps */
    	fp=f;
		for(int j=0; j<m_order-1; j++){
			f[j]=fp[j+1];
	    }
	    m_dyn->evaluate(t+h, xp, dx);
	    f[m_order-1]=dx;

		correction(h,x,f,xp); 

		/* Updating the last saved step */
	    m_dyn->evaluate(t+h, xp, dx);
	    f[m_order-1]=dx;
	    
		/* Saving states */
		t+=h;
		x=xp;
		t_history.push_back(t);
		x_history.push_back(x);
	}

	return 0;
}

template < class T >
int ABM6<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

	std::vector<std::vector<T> > x_history;
	std::vector<double> t_history;

	integrate(ti,tend,nsteps,x0,x_history,t_history);

	xfinal=x_history.back();

	return 0;
}


template class ABM6<double>;
template class ABM6<float>;
template class ABM6<long double>;
#ifdef ENABLE_SMARTUQ
template class ABM6<smartuq::polynomial::chebyshev_polynomial<double> >;
template class ABM6<smartuq::polynomial::chebyshev_polynomial<float> >;
template class ABM6<smartuq::polynomial::chebyshev_polynomial<long double> >;
template class ABM6<smartuq::polynomial::taylor_polynomial<double> >;
template class ABM6<smartuq::polynomial::taylor_polynomial<float> >;
template class ABM6<smartuq::polynomial::taylor_polynomial<long double> >;
#endif
