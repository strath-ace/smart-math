#include "../../include/Integrators/ABM.h"
#include "../../include/Integrators/AB.h"


using namespace smartmath;
using namespace smartmath::integrator;

template < class T >
ABM<T>::ABM(const dynamics::base_dynamics<T> *dyn, const int order) : base_integrator<T>("Adam Bashforth Moulton integration scheme with user-defined order", dyn), m_order(order)
{

	if((order<1)||(order>8))
    	smart_throw("order must be between 1 and 8");    

	double gamma[9]={1.0,-1.0/2.0,-1.0/12.0,-1.0/24.0,-19.0/720.0,-3.0/160.0,-863.0/60480.0,-275.0/24192.0,-33953.0/3628800.0};
	for(int i=0; i<=m_order; i++)
		m_gamma.push_back(gamma[i]);

}

template < class T >
ABM<T>::~ABM(){

}

template < class T >
int ABM<T>::backward_differences(const std::vector<std::vector<T> > &f, const int &m, std::vector<std::vector<T> > &Df) const{

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
int ABM<T>::correction(const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &Df, std::vector<T> &xfinal) const{

	xfinal=x0;
	for(int i=0; i<x0.size(); i++){
		for(int j=0; j<m; j++){		
			xfinal[i]+=h*m_gamma[j]*Df[j][i];
	    }
	}

	return 0;
}

template < class T >
int ABM<T>::integration_step(const double &t, const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal) const{
    	
	if(f.size()!=m)
    	smart_throw("wrong number of saved states in multistep integration"); 

	integrator::AB<T> predictor(m_dyn, m);	    	
	std::vector<T> x(x0), dx(x0);
	std::vector<std::vector<T> > Df;

	/* Prediction */
    predictor.backward_differences(f,m,Df);
    predictor.integration_step(m,h,x0,Df,x);

	/* Updating the saved steps */
   	std::vector<std::vector<T> > fp=f;
	for(int j=0; j<m-1; j++)
		fp[j]=f[j+1];
   	m_dyn->evaluate(t+h, x, dx);
	fp[m-1]=dx;

	/* Correction */
	backward_differences(fp,m,Df);
	correction(m,h,x0,Df,xfinal);	

	return 0;
}

template < class T >
int ABM<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history) const{

	t_history.clear();
	x_history.clear();

	std::vector<T> x(x0),xp(x0),dx(x0);
	std::vector<std::vector<T> > f, fp, Df;
	integrator::AB<T> predictor(m_dyn, m_order);	
	double t=ti, h = (tend-ti)/nsteps;

	predictor.initialize(m_order,ti,h,x0,f); // initializing predictor-corrector by initializing predictor

    for(int k=0; k<nsteps; k++){

    	integration_step(t,m_order,h,x,f,xp);
		x=xp;

		/* Updating the saved steps */
		fp=f;
		for(int j=0; j<m_order-1; j++)
			f[j]=fp[j+1];
   		m_dyn->evaluate(t+h, xp, dx);
   		f[m_order-1]=dx;

   		/* Saving states */
		t+=h;
		t_history.push_back(t);
   		x_history.push_back(x);
	}

	return 0;
}

template < class T >
int ABM<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

	std::vector<std::vector<T> > x_history;
	std::vector<double> t_history;

	integrate(ti,tend,nsteps,x0,x_history,t_history);

	xfinal=x_history.back();

	return 0;
}

template class ABM<double>;
template class ABM<float>;
template class ABM<long double>;
#ifdef ENABLE_SMARTUQ
template class ABM<smartuq::polynomial::chebyshev_polynomial<double> >;
template class ABM<smartuq::polynomial::chebyshev_polynomial<float> >;
template class ABM<smartuq::polynomial::chebyshev_polynomial<long double> >;
template class ABM<smartuq::polynomial::taylor_polynomial<double> >;
template class ABM<smartuq::polynomial::taylor_polynomial<float> >;
template class ABM<smartuq::polynomial::taylor_polynomial<long double> >;
#endif
