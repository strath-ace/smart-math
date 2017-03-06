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
            std::string m_name;
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
            std::string get_name(){return m_name;}

            /**
             *
             */
            int get_value(){return m_value;}

            /**
             *
             */
            bool get_trigger(){return m_trigger;}

            /**
             *
             */
            virtual void switch_trigger_on(const double &t) = 0;

            /**
             *
             */
            virtual void switch_trigger_off() = 0;

            /**
             *
             */
            double get_last_time(){return m_last_time;}

            /**
             *
             */
            void set_last_time(const double &t){m_last_time = t;}

            /**
            *
            */
            virtual int evaluate(const double &t, const std::vector<T> &state) = 0;
        };
    }
}

#endif // SMARTMATH_BASE_EVENT_H
