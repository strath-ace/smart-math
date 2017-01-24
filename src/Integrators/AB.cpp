#include "../../include/Integrators/AB.h"
#include "../../include/Integrators/rk4.h"


using namespace smartmath;
using namespace smartmath::integrator;

template < class T >
AB<T>::AB(const dynamics::base_dynamics<T> *dyn, const int order) : base_integrator<T>("Adam Bashford integration scheme with user-defined order", dyn), m_order(order)
{

	if((order<1)||(order>8))
    	smart_throw("order must be between 1 and 8");    

	double gamma[9]={1.0,1.0/2.0,5.0/12.0,3.0/8.0,251.0/720.0,95.0/288.0,19087.0/60480.0,5257.0/17280.0,1070017.0/3628800.0};
	for(int i=0; i<=m_order; i++)
		m_gamma.push_back(gamma[i]);
	
}

template < class T >
AB<T>::~AB(){

}

template < class T >
int AB<T>::backward_differences(const std::vector<std::vector<T> > &f, const int &m, std::vector<std::vector<T> > &Df) const{

	if(f.size()!=m)
    	smart_throw("wrong number of saved states in multistep integration"); 

	Df.clear();
	Df.push_back(f[m-1]);

	if(m>1){
		std::vector<std::vector<T> > fp,Dfp;
		for(int j=0; j<m-1; j++)
			fp.push_back(f[j]);
		backward_differences(fp,m-1,Dfp); // recursive call
		for(int j=1; j<m; j++){
			Df.push_back(Df[j-1]);
			for(int i=0; i<f[0].size(); i++){
				Df[j][i]-=Dfp[j-1][i];
			}
		}
	}

	return 0;
}

template < class T >
int AB<T>::integration_step(const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &Df, std::vector<T> &xfinal) const{

	xfinal=x0;
	for(int i=0; i<x0.size(); i++){
		for(int j=0; j<m; j++){		
			xfinal[i]+=h*m_gamma[j]*Df[j][i];
	    }
	}

	return 0;
}

template < class T >
int AB<T>::initialize(const int &m, const double &ti, const double &h, const std::vector<T> &x0, std::vector<std::vector<T> > &f) const{

	f.clear();

	std::vector<T> dx(x0), x(x0), xp(x0);
	std::vector< std::vector<T> > fp;

	integrator::rk4<T> RK(m_dyn); // Runge kutta schemed used for initialization (here RK4)

	/* Computing the initial saved steps */
	m_dyn->evaluate(ti,x,dx);
	fp.push_back(dx);
	double t=ti;
	for(int j=0; j<m-1; j++){
		RK.integration_step(t,-h,x,xp);
		t-=h;
		x=xp;
		m_dyn->evaluate(t,x,dx);
		fp.push_back(dx);
	}

	/* Putting the saved steps in the right order */
	for(int j=0; j<m; j++)
		f.push_back(fp[m-1-j]);

	return 0;
}

template < class T >
int AB<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history) const{
	
	t_history.clear();
	x_history.clear();

	std::vector<T> x(x0),xp(x0),dx(x0);
	std::vector<std::vector<T> > f, fp, Df;
	double t=ti, h = (tend-ti)/nsteps;

	initialize(m_order,ti,h,x0,f);

    for(int k=0; k<nsteps; k++){

		backward_differences(f,m_order,Df);
		integration_step(m_order,h,x,Df,xp);
		
    	/* Updating saved steps */
		fp=f;
		for(int j=0; j<m_order-1; j++)
			f[j]=fp[j+1];
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

template<class T>
int AB<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

	std::vector<std::vector<T> > x_history;
	std::vector<double> t_history;

	integrate(ti,tend,nsteps,x0,x_history,t_history);

	xfinal=x_history.back();
    
	return 0;
}


template class AB<double>;
template class AB<float>;
template class AB<long double>;
#ifdef ENABLE_SMARTUQ
template class AB<smartuq::polynomial::chebyshev_polynomial<double> >;
template class AB<smartuq::polynomial::chebyshev_polynomial<float> >;
template class AB<smartuq::polynomial::chebyshev_polynomial<long double> >;
template class AB<smartuq::polynomial::taylor_polynomial<double> >;
template class AB<smartuq::polynomial::taylor_polynomial<float> >;
template class AB<smartuq::polynomial::taylor_polynomial<long double> >;
#endif
