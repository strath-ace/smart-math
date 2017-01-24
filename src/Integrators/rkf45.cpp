#include "../../include/Integrators/rkf45.h"
#include <iomanip>  

using namespace smartmath;
using namespace smartmath::integrator;

template<class T>
rkf45<T>::rkf45(const dynamics::base_dynamics<T> *dyn, const double tol, const double multiplier, const double minstep_events) : base_integrator<T>("Runge Kutta 4-5 variable step time", dyn),m_tol(tol),m_multiplier(multiplier),	m_minstep_events(minstep_events)
{
	/** sanity checks **/
	if(tol<=0.0)
       smart_throw("tolerance for estimated error must be non negative");
    if((multiplier>5.0)||(multiplier<2.0))
       smart_throw("maximum step-multiplier must be between 2 and 5");
	if(minstep_events<=0.0)
       smart_throw("minimum step for events must be non negative");

}

template<class T>
rkf45<T>::~rkf45(){

}

template<class T>
int rkf45<T>::integration_step(const double &t, const double &h, const std::vector<T> &x, std::vector<T> &xtemp, T &er) const{
		
		int n = x.size();
		std::vector<T> k1(x), k2(x), k3(x), k4(x), k5(x), k6(x), xbar(x);
		double t1, t2, t3, t4, t5, t6;

		t1 = t;
		t2 = t + h/4.0;
		t3 = t + h*3.0/8.0;
		t4 = t + h*12.0/13.0;
		t5 = t + h;
		t6 = t + h/2.0;

		//* Evaluate k1 
		m_dyn->evaluate(t1, x, k1);

		//* Evaluate k2 
		for(int j=0; j<n; j++)
		    xtemp[j] = x[j]+k1[j]*h/4.0;
		m_dyn->evaluate(t2, xtemp, k2);

		//* Evaluate k3 
		for(int j=0; j<n; j++)
		    xtemp[j] = x[j]+k1[j]*h*3.0/32.0+k2[j]*h*9.0/32.0;
		m_dyn->evaluate(t3, xtemp, k3);

		//* Evaluate k4 
		for(int j=0; j<n; j++)
		    xtemp[j] = x[j]+k1[j]*h*1932.0/2197.0-k2[j]*h*7200.0/2197.0+k3[j]*h*7296.0/2197.0;
		m_dyn->evaluate(t4, xtemp, k4);

		//* Evaluate k5
		for(int j=0; j<n; j++)
		    xtemp[j] = x[j]+k1[j]*h*439.0/216.0-k2[j]*h*8.0+k3[j]*h*3680.0/513.0-k4[j]*h*845.0/4104.0;
		m_dyn->evaluate(t5, xtemp, k5);

		//* Evaluate k6
		for(int j=0; j<n; j++)
		    xtemp[j] = x[j]-k1[j]*h*8.0/27.0+k2[j]*h*2.0-k3[j]*h*3544.0/2565.0+k4[j]*h*1859.0/4104.0-k4[j]*h*11.0/40.0;
		m_dyn->evaluate(t6, xtemp, k6);

		//* Return x(t+h) computed from fourth-order Runge Kutta.
		er=0.0*x[0];
		for(int j=0; j<n; j++){
		    xbar[j] = x[j]+ (k1[j]*16.0/135.0+k3[j]*6656.0/12825.0+k4[j]*28561.0/56430.0-k5[j]*9.0/50.0+k6[j]*2.0/55.0)*h;
		    xtemp[j] = x[j] + (k1[j]*25.0/216.0+k3[j]*1408.0/2565.0+k4[j]*2197.0/4104.0-k5[j]/5.0)*h;
		    er+=pow(xbar[j]-xtemp[j],2);
		}
		er=sqrt(er);

	return 0;
}


template<class T>
std::vector<int> rkf45<T>::dummy_event(std::vector<T> x, double d){

	std::vector<int> output(1,0);

	return output;

}

template<class T>
int rkf45<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> &t_history) const{

	double t0=ti, tf=tend, n=nsteps;
	std::vector<T> x(x0);

	integrate(t0,tf,n,x,x_history,t_history,dummy_event);
	
	return 0;
}

template<class T>
int rkf45<T>::integrate(const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal) const{

	std::vector<std::vector<T> > x_history;
	std::vector<double> t_history;

	integrate(ti, tend, nsteps, x0, x_history, t_history);

	xfinal=x_history.back();

	return 0;
}


template<class T>
int rkf45<T>::integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<std::vector<T> > &x_history, std::vector<double> t_history, std::vector<int> (*g)(std::vector<T> x, double d)) const{

	x_history.clear();
	t_history.clear();

    int k;
    int check=0;
    std::vector<T> x(x0), xtemp(x0);

	double factor=1.0, value=0.0, t=ti, h = (tend-ti)/nsteps;
	T er=0.0*x0[0];

	std::vector<int> events=g(x0,ti), events2=events;
	int m=events.size();

	int i=0;
    while(sqrt(pow(t-ti,2))<sqrt(pow(tend-ti,2))){

		if(sqrt(pow(tend-t,2))<sqrt(h*h))
			h=tend-t;

		integration_step(t,h,x,xtemp,er);

		/* Step-size control */
		error(er,value);
		factor=pow(m_tol/value,1.0/(4.0+1.0));
		if(value>m_tol){ // unsucessful step
			h*=0.9*factor; 	
		}
		else{ // sucessful step
			/* Checking for the events */
			events2=g(xtemp,t+h);
			k=0;
			check=0;		
			while(k<m){
				if(events2[k]-events[k]!=0){
					check=1;
					k=m;
				}
				else{
					k++;
				}
			}

			if(check==1){
				if(sqrt(h*h)>m_minstep_events){
					h*=0.5;
				}
				else{
					tend=t+h; // saving the termination time	
					t=tend; // trick to get out of the while loop	
					x_history.push_back(xtemp);
					t_history.push_back(t);
					std::cout << "Propagation interrupted by terminal event at time " << tend << " after " << i << " steps" << std::endl;
				}
			}
			else{
				x=xtemp; // updating state
				t+=h; // updating current time			
				events=events2;				
				x_history.push_back(x);
				t_history.push_back(t);	
				/* Step-size control */
				if(factor>m_multiplier){
					factor=m_multiplier;
				}	
				h*=factor; // updating step-size
				i++; // counting the number of steps
			}			
		
		}
	}

	return 0;
}

template<class T>
int rkf45<T>::integrate(const double &ti, double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal, std::vector<int> (*g)(std::vector<T> x, double d)) const{

	std::vector<std::vector<T> > x_history;
	std::vector<double> t_history;

	integrate(ti, tend, nsteps, x0, x_history, t_history, *g);

	xfinal=x_history.back();

	return 0;
}

template<class T>
int rkf45<T>::error(const float &x, double &val) const{
	val=x;
return 0;
}

template<class T>
int rkf45<T>::error(const double &x, double &val) const{
	val=x;
return 0;
}

template<class T>
int rkf45<T>::error(const long double &x, double &val) const{
	val=x;
return 0;
}

#ifdef ENABLE_SMARTUQ

template<class T>
int rkf45<T>::error(const smartuq::polynomial::chebyshev_polynomial<double> &x, double &val) const{
	val=x.get_range()[1];
return 0;
}

template<class T>
int rkf45<T>::error(const smartuq::polynomial::chebyshev_polynomial<float> &x, double &val) const{
	val=x.get_range()[1];
return 0;
}

template<class T>
int rkf45<T>::error(const smartuq::polynomial::chebyshev_polynomial<long double> &x, double &val) const{
	val=x.get_range()[1];
return 0;
}

template<class T>
int rkf45<T>::error(const smartuq::polynomial::taylor_polynomial<double> &x, double &val) const{
	val=x.get_coeffs()[0];
return 0;
}

template<class T>
int rkf45<T>::error(const smartuq::polynomial::taylor_polynomial<float> &x, double &val) const{
	val=x.get_coeffs()[0];
return 0;
}

template<class T>
int rkf45<T>::error(const smartuq::polynomial::taylor_polynomial<long double> &x, double &val) const{
	val=x.get_coeffs()[0];
return 0;
}

#endif

template class rkf45<double>;
template class rkf45<float>;
template class rkf45<long double>;
#ifdef ENABLE_SMARTUQ
template class rkf45<smartuq::polynomial::chebyshev_polynomial<double> >;
template class rkf45<smartuq::polynomial::chebyshev_polynomial<float> >;
template class rkf45<smartuq::polynomial::chebyshev_polynomial<long double> >;
template class rkf45<smartuq::polynomial::taylor_polynomial<double> >;
template class rkf45<smartuq::polynomial::taylor_polynomial<float> >;
template class rkf45<smartuq::polynomial::taylor_polynomial<long double> >;
#endif
