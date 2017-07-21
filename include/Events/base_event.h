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
            /**
             * @brief m_name Event name.
             */
            std::string m_name;


            /**
             * @brief m_value Event value.
             *
             * Output value of the event after the last evaluation.
             */
            int m_value;


            /**
             * @brief m_last_time Time of the last evaluation.
             */
            double m_last_time;


            /**
             * @brief m_status_scope Evaluation status scope.
             *
             * This attribute determines whether an event needs to be evaluated or not without
             * removing it from the integration method or modifying the list of events.
             * The evaluation of a method might be unnecessary for external reasons (e.g. the
             * sensor associated to this event might be off).
             */
            bool* m_status_scope = NULL;


            /**
             * @brief m_trigger Event trigger.
             *
             * It determines whether the event has been triggered or not.
             * This value is changed by the integrator according to its event handling method.
             */
            bool m_trigger = false;


            /**
             * @brief m_comments Print comments status.
             */
            bool m_comments = false;


        public:
            /**
             * @brief ~base_event base_event destructor (virtual).
             */
            virtual ~base_event(){}


            /**
             * @brief get_name returns event name.
             */
            std::string get_name(){return m_name;}


            /**
             * @brief get_value returns event value.
             */
            int get_value(){return m_value;}


            /**
             * @brief get_last_time returns the last evaluation time.
             */
            double get_last_time(){return m_last_time;}

            /**
             * @brief set_last_time sets the last evaluation time.
             *
             * @param t time.
             */
            void set_last_time(const double &t){m_last_time = t;}


            /**
             * @brief get_status_scope returns pointer to the status scope.
             */
            bool* get_status_scope(){return m_status_scope;}


            /**
             * @brief set_status_scope sets the status scope.
             *
             * @param scope pointer to the status scope.
             */
            void set_status_scope(bool* scope){m_status_scope = scope;}


            /**
             * @brief get_trigger returns trigger.
             */
            bool get_trigger(){return m_trigger;}

            /**
             * @brief switch_trigger_on switches on the trigger.
             *
             * This method is virtual. All the classes that inherit from this must
             * provide an implementation.
             * @param t time.
             */
            virtual void switch_trigger_on(const double &t) = 0;

            
            /**
             * @brief switch_trigger_off switches off the trigger.
             *
             * This method is virtual. All the classes that inherit from this must
             * provide an implementation.
             */
            virtual void switch_trigger_off() = 0;


            /**
             * @brief set_comments sets the status of the comments switch.
             */
            void set_comments(const bool &status){m_comments = status;}

            /**
             * @brief get_comments returns the status of the comments switch.
             */
            bool get_comments(){return m_comments;}


            /**
             * @brief evaluate evaluates the event.
             *
             * This method is virtual. All the classes that inherit from this must
             * provide an implementation.
             * @param t time.
             * @param state state vector.
             */
            virtual int evaluate(const double &t, const std::vector<T> &state) = 0;
        };
    }
}

#endif // SMARTMATH_BASE_EVENT_H
