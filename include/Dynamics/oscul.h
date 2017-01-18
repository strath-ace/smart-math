/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef SMARTMATH_OSCUL_H
#define SMARTMATH_OSCUL_H

#include "base_dynamics.h"
#include "../exception.h"

namespace smartmath
{
    namespace dynamics {

        template < class T >
        class oscul: public base_dynamics<T>
        {

        private:
            using base_dynamics<T>::m_name;

        public:
            oscul(const std::vector<T> &param = std::vector<T>(4));

            ~oscul();

            int evaluate(const double &E, const std::vector<T> &state, std::vector<T> &dstate) const;

        private:
            mutable std::vector<T> m_param;

        };

    }
}


#endif // SMARTMATH_OSCUL_H
