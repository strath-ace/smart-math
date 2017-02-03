/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
-------------------- Author: Francesco Torre -------------------------
-------------- e-mail: francesco.torre@strath.ac.uk ------------------
*/

#include "../../include/Dynamics/spring.h"

using namespace smartmath;
using namespace dynamics;

spring::spring(const double &k, const double &t_scale, const double &r_scale) : base_dynamics<double>("Spring Problem"),
    m_k(k), m_t_scale(t_scale), m_r_scale(r_scale)
{
    if(m_k < 0.0)
        smart_throw(this->m_name+": the elastic coefficient must be >= 0");

}

spring::~spring()
{

}

int spring::evaluate(const double &t, const std::vector<double> &state, std::vector<double> &dstate) const
{
    //sanity checks
    if(t<0)
        smart_throw(this->m_name+": negative time supplied in evaluation of the dynamical system");
    if(state.size() > 6)
        smart_throw(this->m_name+": it's difficult to imagine a spring acting ot time or higher dimensions...");
    if(state.size()%2 != 0)
        smart_throw(this->m_name+": half dimensions are not contemplated in this version. Wait for update #42.");

    unsigned int dimension = state.size()/2;

    dstate.clear();

    // velocities
    for(unsigned int index = 0; index < dimension; ++index)
    {
        dstate.push_back(state[dimension+index]);
    }

    // acceerations
    for(unsigned int index = 0; index < dimension; ++index)
    {
        dstate.push_back(-m_k*state[index]);
    }

    if(dstate.size() != state.size())
        smart_throw(this->m_name+": state and dstate dimensions mismatch.");

    return 0;

}
