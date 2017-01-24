#include "../../include/Integrators/PECEvar.h"
#include "../../include/Integrators/AB.h"


using namespace smartmath;
using namespace smartmath::integrator;

template < class T >
PECEvar<T>::PECEvar(const dynamics::base_dynamics<T> *dyn, const int order_max, const int order_min, const double tol, const double multiplier, const double minstep_events) : base_integrator<T>("Adam Bashforth Moulton integration scheme with variable order", dyn), m_order_max(order_max), m_order_min(order_min), m_tol(tol), m_multiplier(multiplier), m_minstep_events(minstep_events)
{

	if(order_min>order_max)
    	smart_throw("minimum order must be smaller or equal to maximum order");
	if((order_min<1)||(order_min>8))
    	smart_throw("minimum order must be between 1 and 8");    
	if((order_max<1)||(order_max>8))
    	smart_throw("maximum order must be between 1 and 8");
	if(tol<=0.0)
       smart_throw("tolerance for estimated error must be non negative");
    if((multiplier>5.0)||(multiplier<2.0))
       smart_throw("maximum step-multiplier must be between 2 and 5");
	if(minstep_events<=0.0)
       smart_throw("minimum step for events must be non negative");

	double gamma[9]={1.0,-1.0/2.0,-1.0/12.0,-1.0/24.0,-19.0/720.0,-3.0/160.0,-863.0/60480.0,-275.0/24192.0,-33953.0/3628800.0};
	for(int i=0; i<=m_order_max; i++)
		m_gamma.push_back(gamma[i]);

}

template < class T >
PECEvar<T>::~PECEvar(){

}

template < class T >
int PECEvar<T>::backward_differences(const std::vector<std::vector<T> > &f, const int &m, std::vector<std::vector<T> > &Df) const{

	if(f.size()!=m)
    	smart_throw("wrong number of saved states in multistep integration"); 

	Df.clear();
	Df.push_back(f[m-1]);

	if(m>1){
		std::vector<std::vector<T> > fp,Dfp;
		for(int j=0; j<m-1; j++)
			fp.push_back(f[j]);
		backward_differences(fp,m-1,Dfp);
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
int PECEvar<T>::correction(const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &Df, std::vector<T> &xfinal) const{

	xfinal=x0;
	for(int i=0; i<x0.size(); i++){
		for(int j=0; j<m; j++){		
			xfinal[i]+=h*m_gamma[j]*Df[j][i];
	    }
	}

	return 0;
}

template < class T >
int PECEvar<T>::integration_step(const double &t, const int &m, const double &h, const std::vector<T> &x0, const std::vector<std::vector<T> > &f, std::vector<T> &xfinal, T &er) const{
    	
	if(f.size()!=m)
    	smart_throw("wrong number of saved states in multistep integration"); 

	integrator::AB<T> predictor(m_dyn, m);	    	
	std::vector<T> x(x0), dx(x0);
	std::vector<std::vector<T> > Df;

    predictor.backward_differences(f,m,Df);
    predictor.integration_step(m,h,x0,Df,x);
   	std::vector<std::vector<T> > fp=f;
	for(int j=0; j<m-1; j++)
		fp[j]=f[j+1];
   	m_dyn->evaluate(t+h, x, dx);
	fp[m-1]=dx;
	backward_differences(fp,m,Df);
	correction(m,h,x0,Df,xfinal);

	unsigned int l=x.size();
	er=0.0;
	for(unsigned int i=0; i<l; i++)
			er+=pow(xfinal[i]-x[i],2);
	er=sqrt(er);	

	return 0;
}

template < class T >
int PECEvar<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history) const{

	double t0=ti, tf=tend, n=nsteps;
	std::vector<T> x(x0);

	integrate(t0,tf,n,x,x_history,t_history,dummy_event);

	return 0;
}

template < class T >
int PECEvar<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

	std::vector<std::vector<T> > x_history;
	std::vector<double> t_history;

	integrate(ti,tend,nsteps,x0,x_history,t_history);

	xfinal=x_history.back();

	return 0;
}

template<class T>
int PECEvar<T>::integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> t_history, std::vector<int> (*g)(std::vector<T> x, double d)) const{

	t_history.clear();
	x_history.clear();

    int i;
    int check=0;
    double factor;	

	std::vector<T> x(x0),xp(x0),dx(x0);
	std::vector<std::vector<T> > f_max, f, Df;
	T er;
	integrator::AB<T> predictor(m_dyn, m_order_max);	
	int m=m_order_min, mold=m_order_min;
	double t=ti, h = (tend-ti)/nsteps;
	double value;

	std::vector<int> events=g(x0,ti), events2=events;
	int q=events.size();

	predictor.initialize(m_order_max,ti,h,x0,f_max);

	int k=0;
    while(sqrt(pow(t-ti,2))<sqrt(pow(tend-ti,2))){

		if(sqrt(pow(tend-t,2))<sqrt(h*h)){
			h=tend-t;
			predictor.initialize(m_order_max,t,h,x,f_max);	
		}
			
		f.clear();
		for(int j=0; j<m; j++)
			f.push_back(f_max[m_order_max-m+j]);  	

    	integration_step(t,m,h,x,f,xp,er);
		
		/* Step-size and order control */
		error(er,value);
		factor=pow(m_tol/value,1.0/(double(m)+1.0));
		if(value>m_tol){ // unsucessful step
			if(m<m_order_max){
				mold=m;
				m++;
				//std::cout << "increasing order" << std::endl;
			}
			else{				
				h*=0.9*factor;			
				predictor.initialize(m_order_max,t,h,x,f_max);
				//std::cout << "decreasing time-step" << std::endl;	
			}
		}
		else{ // sucessful step
			/* Checking for the events */
			events2=g(xp,t+h);
			i=0;
			check=0;		
			while(i<q){
				if(events2[i]-events[i]!=0){
					check=1;
					i=q;
				}
				else{
					i++;
				}
			}
			if(check==1){
				if(sqrt(h*h)>m_minstep_events){
					h*=0.5;
					predictor.initialize(m_order_max,t,h,x,f_max);
					//std::cout << "decreasing time-step for event detection" << std::endl;
				}
				else{
					tend=t+h; // saving the termination time	
					t=tend; // trick to get out of the while loop	
					x_history.push_back(xp);
					t_history.push_back(t);
					std::cout << "Propagation interrupted by terminal event at time " << tend << " after " << k << " steps" << std::endl;
				}
			}
			else{
				k++; // counting the number of steps
				f=f_max;
				for(int j=0; j<m_order_max-1; j++)
					f_max[j]=f[j+1];
	    		m_dyn->evaluate(t+h, xp, dx);
	    		f_max[m_order_max-1]=dx;
				x=xp; // updating state
				t+=h; // updating current time			
				events=events2;				
				x_history.push_back(x);
				t_history.push_back(t);	
				/* Step-size and order control */
				if((m>m_order_min)&&(mold!=m-1)){
					mold=m;
					m--;
					//std::cout << "decreasing order" << std::endl;
				}
				else if((m==m_order_min)&&(factor>=2.0)){
					if(factor>m_multiplier){
						factor=m_multiplier;
					}	
					h*=factor; // updating step-size
					predictor.initialize(m_order_max,t,h,x,f_max);
					//std::cout << "increasing time-step" << std::endl;
				}
				else{
					//std::cout << "keeping both time-step and order constant" << std::endl;
				}
				//std::cout << k << std::endl;
			}	
		}	
	}

	return 0;
}

template<class T>
int PECEvar<T>::integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal, std::vector<int> (*g)(std::vector<T> x, double d)) const{

	std::vector<std::vector<T> > x_history;
	std::vector<double> t_history;

	integrate(ti, tend, nsteps, x0, x_history, t_history, *g);

	xfinal=x_history.back();

	return 0;
}

template<class T>
std::vector<int> PECEvar<T>::dummy_event(std::vector<T> x, double d){

	std::vector<int> output(1,0);

	return output;

}

template<class T>
int PECEvar<T>::error(const float &x, double &val) const{
	val=x;
return 0;
}

template<class T>
int PECEvar<T>::error(const double &x, double &val) const{
	val=x;
return 0;
}

template<class T>
int PECEvar<T>::error(const long double &x, double &val) const{
	val=x;
return 0;
}

#ifdef ENABLE_SMARTUQ

template<class T>
int PECEvar<T>::error(const smartuq::polynomial::chebyshev_polynomial<double> &x, double &val) const{
	val=x.get_range()[1];
return 0;
}

template<class T>
int PECEvar<T>::error(const smartuq::polynomial::chebyshev_polynomial<float> &x, double &val) const{
	val=x.get_range()[1];
return 0;
}

template<class T>
int PECEvar<T>::error(const smartuq::polynomial::chebyshev_polynomial<long double> &x, double &val) const{
	val=x.get_range()[1];
return 0;
}

template<class T>
int PECEvar<T>::error(const smartuq::polynomial::taylor_polynomial<double> &x, double &val) const{
	val=x.get_coeffs()[0];
return 0;
}

template<class T>
int PECEvar<T>::error(const smartuq::polynomial::taylor_polynomial<float> &x, double &val) const{
	val=x.get_coeffs()[0];
return 0;
}

template<class T>
int PECEvar<T>::error(const smartuq::polynomial::taylor_polynomial<long double> &x, double &val) const{
	val=x.get_coeffs()[0];
return 0;
}

#endif

template class PECEvar<double>;
template class PECEvar<float>;
template class PECEvar<long double>;
#ifdef ENABLE_SMARTUQ
template class PECEvar<smartuq::polynomial::chebyshev_polynomial<double> >;
template class PECEvar<smartuq::polynomial::chebyshev_polynomial<float> >;
template class PECEvar<smartuq::polynomial::chebyshev_polynomial<long double> >;
template class PECEvar<smartuq::polynomial::taylor_polynomial<double> >;
template class PECEvar<smartuq::polynomial::taylor_polynomial<float> >;
template class PECEvar<smartuq::polynomial::taylor_polynomial<long double> >;
#endif
