#include "../../include/Dynamics/oscul.h"

using namespace smartmath;
using namespace dynamics;

template < class T >
oscul<T>::oscul(const std::vector<T> &param) : base_dynamics<T>("drag on orbital elements"),
    m_param(param)
{
    
    if(m_param.size()!=3)
        smart_throw(m_name+": the parameters list need to be of size 3");

}

template < class T >
oscul<T>::~oscul()
{

}

template < class T >
int oscul<T>::evaluate(const double &E, const std::vector<T> &state, std::vector<T> &dstate) const
{
    //sanity checks
    if(state.size()!=4)
        smart_throw(m_name+": the state dimension needs to be 4");

    dstate.clear();

    double omega_earth = 7.2921150e-5;// * m_t_scale;
    double radius_earth = 6378.137e3;
    double h0 = 120.0e3/radius_earth;
    double ht = 1000.0e3/radius_earth;
    double mu=398600.4415e9/pow(radius_earth,3);

    T a = state[0];
    T e = state[1];
    T i = state[2];
    T w = state[3];

    T delta = m_param[0]; //CD*S/m!!
    T rho_0 = 2.438e-8*m_param[1];
    T rho_t = 3.019e-15*m_param[2];

    T H = (h0-ht)/log(rho_t/rho_0);

    T Fminus = 1.0-e*cos(E);
    T Fplus = 2.0-Fminus;
    T eta = sqrt(1.0-e*e);
    T eta2 = 1.0-e*e;
    T fcos = (cos(E)-e)/Fminus;
    T fsin = (eta*sin(E))/Fminus;
    T apow=pow(a,1.5);
    T tau = omega_earth* apow* eta * cos(i) / sqrt(mu);

    T B = 1.0 - tau * Fminus / Fplus ;

    T inter =sqrt(Fminus/Fplus);

    T Sa = B*B*Fplus/inter;
    T Se = eta2*cos(E)*Sa/Fplus + 0.5*e*tau*B*Fminus*inter * sin(E)*sin(E);
    T Zcos = cos(w)*fcos - sin(w)*fsin;
    T Zsin = cos(w)*fsin + sin(w)*fcos;
    T Si = B*Fminus*Fminus*Fplus*inter*Zcos *Zcos;
    T Sw = sin(E)*(eta*Sa/Fplus-(e+cos(E))*(0.5*e*tau/eta)*inter*B*Fminus)-0.5*e*tau*(1.0/eta2)*sin(i)*B*Fminus*Fminus*Fplus*inter*Zcos*Zsin;

    T h = a*Fminus-1.0; // a*Fminus-(1.0-(1.0/298.0)*sin(i)*sin(i)*Zsin*Zsin);
    T factor = a*radius_earth*delta* rho_0 * exp(-(h-h0)/H);

    dstate.push_back( - a * factor * Sa  );
    dstate.push_back( - factor * Se  );
    dstate.push_back( - factor * omega_earth * apow * sin(i) * Si / (2.0*sqrt(mu)*eta)  );
    dstate.push_back( - factor * Sw / e  );

    return 0;

}


template class oscul<double>;
template class oscul<float>;
template class oscul<long double>;
#ifdef ENABLE_SMARTUQ
template class oscul<polynomial::chebyshev_polynomial<double> >;
template class oscul<polynomial::chebyshev_polynomial<float> >;
template class oscul<polynomial::chebyshev_polynomial<long double> >;
template class oscul<polynomial::taylor_polynomial<double> >;
template class oscul<polynomial::taylor_polynomial<float> >;
template class oscul<polynomial::taylor_polynomial<long double> >;
#endif
