/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/

#ifndef SMARTMATH_INLINEFUNCTIONS_H
#define SMARTMATH_INLINEFUNCTIONS_H

#include <vector>
#include <cmath>
#include <iostream>
#include "../LinearAlgebra/Eigen/Eigen"
#include "../LinearAlgebra/Eigen/Core"
#include "../LinearAlgebra/EigenMultivariateNormal.h"
#include "../exception.h"
#include "chebyshev_inequality.h"
#include <time.h>

namespace smartmath{

const double ZERO = 1.0e-15;

template <class T>
T inverse(T x){
    if(fabs(x)<=ZERO){
        std::cout<<"ERROR: Division by zero."<<std::endl;
        throw std::exception();
    }
    return 1.0/x;
}

//MATH STUFFS
int factorial(int n);
int combination(int n, int k);
void rep(std::vector<std::vector<int> > &res, const std::vector<int> &values, std::vector<int> &item, unsigned int count);
void variations(const std::vector<int> values, const int k, std::vector<std::vector<int> > &res);

typedef double (*fun)(double);
double bisection_method(fun f, double lb, double ub, double prec);

double Legendre(int l, int m, double x);
double Legendre_derivative(int l, int m, double x);

int sample_truncated_normal_distribution(const double &lower_bound,
                                                         const double &upper_bound,
                                                         const double &mean,
                                                         const double &sd,
                                                         const unsigned int &N_samples,
                                                         std::vector<double> &result);

Eigen::MatrixXd sample_multivariate_normal_distribution(const Eigen::VectorXd &mean,
                                                        const Eigen::MatrixXd &covar,
                                                        const int &N_samples);

Eigen::MatrixXd sample_truncated_multivariate_normal_distribution(const Eigen::VectorXd &lower_bounds,
                                                         const Eigen::VectorXd &upper_bounds,
                                                         const Eigen::VectorXd &mean,
                                                         const Eigen::MatrixXd &covar,
                                                         const unsigned int &N_samples, double &proportion_valid_samples);

Eigen::MatrixXd sample_truncated_multivariate_normal_distribution(const Eigen::VectorXd &mean,
                                                                  const Eigen::MatrixXd &covar,
                                                                  const double &min_pr_valid_samples,
                                                                  const unsigned int &N_samples,
                                                                  double &pr_valid_samples);

}

#endif // SMARTMATH_INLINEFUNCTIONS_H
