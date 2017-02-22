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

int find_PC_bounds_normal_distribution(const Eigen::VectorXd &mean,
                                       const Eigen::MatrixXd &covar,
                                       const double &min_pr_valid_samples,
                                       Eigen::VectorXd &lower_bounds,
                                       Eigen::VectorXd &upper_bounds);

/**
 * @brief sample_multivariate_normal_distribution Sampling method for multivariate
 * normal distributions of any dimensionality
 * @param mean Vector containing the mean of the distribution
 * @param covar Covariance Matrix
 * @param N_samples Number of samples needed
 * @return Matrix with N_samples columns. Each column represents a sample
 */
Eigen::MatrixXd sample_multivariate_normal_distribution(const Eigen::VectorXd &mean,
                                                        const Eigen::MatrixXd &covar,
                                                        const int &N_samples);

/**
 * @brief sample_truncated_multivariate_normal_distribution Create samples (x1,x2,...xd)
 * of a multivariate normal distribution, bounded to a given multidimensional
 * rectangular region, such as (a1<=x1<=b1, a2<=x2<=b2,...,ad<=xd<=bd)
 * @param lower_bounds Vector of lower bounds (a1,a2,...,ad)
 * @param upper_bounds Vector of upper bounds (b1,b2,...,bd)
 * @param mean Vector containing the mean of the distribution
 * @param covar Covariance Matrix
 * @param N_samples Number of samples needed
 * @param pr_valid_samples Output argument. Gives the probability that a sample has been marked as
 * valid during the sampling.
 * @return Matrix with N_samples columns. Each column represents a sample
 */
Eigen::MatrixXd sample_truncated_multivariate_normal_distribution(const Eigen::VectorXd &lower_bounds,
                                                         const Eigen::VectorXd &upper_bounds,
                                                         const Eigen::VectorXd &mean,
                                                         const Eigen::MatrixXd &covar,
                                                         const unsigned int &N_samples,
                                                        double &pr_valid_samples);

/**
 * @brief sample_truncated_multivariate_normal_distribution Create samples (x1,x2,...xd)
 * of a multivariate normal distribution, bounded to a multidimensional
 * rectangular region, such as (a1<=x1<=b1, a2<=x2<=b2,...,ad<=xd<=bd). This rectangular region
 * is not given a priori, but it is automatically computed based on the parameter min_pr_valid_samples,
 * which sets the minimum probability that a sample is valid (not rejected) during the sampling process.
 * The higher is this probability, the broader will be the rectangular bounds. Internally, this method
 * creates the samples in the eigenspace corresponding to the given covariance matrix, so that the variables
 * there are independent and we can easily compute the bounds for each dimension independently using the
 * Chebyshev inequality
 * @param mean Vector containing the mean of the distribution
 * @param covar Covariance Matrix
 * @param min_pr_valid_samples Minimum probability that a sample is valid (not rejected) during the sampling process
 * @param N_samples Number of samples needed
 * @param pr_valid_samples Output argument. Actual probability that a sample was valid during the sampling process
 * @return Matrix with N_samples columns. Each column represents a sample
 */
Eigen::MatrixXd sample_truncated_multivariate_normal_distribution(const Eigen::VectorXd &mean,
                                                                  const Eigen::MatrixXd &covar,
                                                                  const double &min_pr_valid_samples,
                                                                  const unsigned int &N_samples,
                                                                  double &pr_valid_samples);

}

#endif // SMARTMATH_INLINEFUNCTIONS_H
