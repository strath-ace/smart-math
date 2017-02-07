/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
--------- Author: Annalisa Riccardi ----------------------------------
*/

#ifndef SMARTMATH_EXCEPTIONS_H
#define SMARTMATH_EXCEPTIONS_H

#include <exception>
#include <cassert>
#include <iostream>
#include <string>

#define _SMARTMATH_EXCEPTION_QUOTEME(x) #x
#define SMARTMATH_EXCEPTION_QUOTEME(x) _SMARTMATH_EXCEPTION_QUOTEME(x)
#define SMARTMATH_EXCEPTION_EXCTOR(s) ((std::string(__FILE__ "," SMARTMATH_EXCEPTION_QUOTEME(__LINE__) ": ") + s) + ".")
#define SMARTMATH_EX_THROW(s) (throw smartmath_exception(SMARTMATH_EXCEPTION_EXCTOR(s)))

#define smartmath_throw(s) SMARTMATH_EX_THROW(s)

namespace smartmath{

class smartmath_exception: public std::exception {
	public:
		smartmath_exception(const std::string &s):m_what(s) {}
		virtual const char *what() const throw() {
			return m_what.c_str();
		}
		virtual ~smartmath_exception() throw() {}
	protected:
		std::string m_what;
};
}
#endif
