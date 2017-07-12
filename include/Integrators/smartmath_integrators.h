/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/


#ifndef SMARTMATH_INTEGRATORS_H
#define SMARTMATH_INTEGRATORS_H

#include "base_integrator.h"
#include "base_rungekutta.h"
#include "euler.h"
#include "midpoint.h"
#include "heun.h"
#include "rk4.h"
#include "base_multistep.h"
#include "AB.h"
#include "ABM.h"
#include "base_integrationwevent.h"
#include "PECEvar.h"
#include "base_stepsizecontrol.h"
#include "rkf45.h"
#include "rk87.h"
#include "base_bulirschstoer.h"
#include "base_symplectic.h"
#include "euler_symplectic.h"
#include "leapfrog.h"
#include "forest.h"
#include "yoshida6.h"
#include "symplectic_mixedvar.h"
#include "euler_mixedvar.h"
#include "leapfrog_mixedvar.h"
#include "forest_mixedvar.h"
#include "yoshida6_mixedvar.h"

#endif // SMARTMATH_INTEGRATORS_H
