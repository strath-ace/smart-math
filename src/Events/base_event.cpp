/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
-------------- e-mail: francesco.torre@strath.ac.uk ------------------
------------------- Author: Francesco Torre --------------------------
*/

#include "../../include/Events/base_event.h"

using namespace smartmath;
using namespace smartmath::events;

// Destructor
template < class T >
base_event<T>::~base_event()
{
    // NOTHING YET
}


// Methods for trigger
template < class T >
bool base_event<T>::get_trigger()
{
    return m_trigger;
}

template < class T >
void base_event<T>::switch_trigger_on(const double &t)
{
    m_trigger = true;
    m_last_time = t;
}

template < class T >
void base_event<T>::switch_trigger_off()
{
    m_trigger = false;
}


// Methods for last time
template < class T >
double base_event<T>::get_last_time()
{
    return m_last_time;
}


template class base_event<double>;
template class base_event<float>;
template class base_event<long double>;
#ifdef ENABLE_SMARTUQ
template class base_event<smartuq::polynomial::chebyshev_polynomial<double> >;
template class base_event<smartuq::polynomial::chebyshev_polynomial<float> >;
template class base_event<smartuq::polynomial::chebyshev_polynomial<long double> >;
template class base_event<smartuq::polynomial::taylor_polynomial<double> >;
template class base_event<smartuq::polynomial::taylor_polynomial<float> >;
template class base_event<smartuq::polynomial::taylor_polynomial<long double> >;
#endif
