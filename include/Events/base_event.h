/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
-------------- e-mail: francesco.torre@strath.ac.uk ------------------
------------------- Author: Francesco Torre --------------------------
*/

#ifndef SMARTMATH_BASE_EVENT_H
#define SMARTMATH_BASE_EVENT_H

#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>

namespace smartmath
{
    namespace events
    {
        template < class T >
        class base_event
        {
        protected:
            int m_value;
            bool m_trigger;
            double m_last_time;
            
        public:
            /**
            *
            */
            virtual ~base_event(){}

            /**
             *
             */
            bool get_trigger(){return m_trigger;}

            /**
             *
             */
            void switch_trigger_on(const double &t){
                m_trigger = true;
                m_last_time = t;
            }

            /**
             *
             */
            void switch_trigger_off(){m_trigger = false;}

            /**
             *
             */
            double get_last_time(){return m_last_time;}

            /**
            *
            */
            virtual int evaluate(const double &t, const std::vector<T> &state) = 0;
        };
    }
}

#endif // SMARTMATH_BASE_EVENT_H
