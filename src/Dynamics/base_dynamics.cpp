#include "../../include/Dynamics/base_dynamics.h"

using namespace smartmath;
using namespace smartmath::dynamics;

template < class T >
base_dynamics<T>::base_dynamics(const std::string &name,
                                const double &r_scale,
                                const double &t_scale)                           
{
    m_name    = name;
    m_r_scale = r_scale;
    m_t_scale = t_scale;
}

template < class T >
base_dynamics<T>::~base_dynamics()
{
    // NOTHING YET    
}

// Methods for name
template < class T >
std::string base_dynamics<T>::get_name() const
{
    return m_name;
}

template class base_dynamics<double>;
template class base_dynamics<float>;
template class base_dynamics<long double>;
#ifdef ENABLE_SMARTUQ
template class base_dynamics<smartuq::polynomial::chebyshev_polynomial<double> >;
template class base_dynamics<smartuq::polynomial::chebyshev_polynomial<float> >;
template class base_dynamics<smartuq::polynomial::chebyshev_polynomial<long double> >;
template class base_dynamics<smartuq::polynomial::taylor_polynomial<double> >;
template class base_dynamics<smartuq::polynomial::taylor_polynomial<float> >;
template class base_dynamics<smartuq::polynomial::taylor_polynomial<long double> >;
#endif
