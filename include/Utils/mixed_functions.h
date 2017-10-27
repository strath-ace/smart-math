/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2016 University of Strathclyde--------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
--------- Author: Annalisa Riccardi and Carlos Ortega Absil ----------
*/

/*! \file mixed_functions.h
  \brief A selection of functions for mathematical purposes
 */ 

/*! \fn double Legendre(int l, int m, double x)
  \brief Legendre evaluation of associated Legendre polynomials
  \param l order
  \param m degree
  \param x evaluation point
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
#include <functional>
#include <utility>


template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
{
  for (unsigned int index = 0; index < v.size(); ++index)
  {
      os << " " << v[index];
  }
  return os;
}


namespace smartmath
{
    //Multinormal sampler. Always use this object to create samples, using the methods setMean, setCovar and samples
    extern Eigen::EigenMultivariateNormal<double> multinormal_sampler;


  /**
   * @brief ZERO numerical value that can be used to check if a double can be considered as zero or not
   */  
  const double ZERO = 1.0e-15;

  /**
   * @brief templated function returning the inverse
   * @param[in] x point where to evaluate the inverse
   * @return value of the inverse of x
   */ 
  template <class T>
  T inverse(T x){
      if(fabs(x)<=ZERO)
      {
          std::cout<<"ERROR: Division by zero."<<std::endl;
          throw std::exception();
      }
      return 1.0/x;
  }

  //MATH STUFFS

  /**
   * @brief method computing the factorial recursively
   * @param[in] n integer whose factorial needs to be computed
   * @return n!
   */   
  int factorial(int n);

  /**
   * @brief method computing (n+k)!/(n!k!)
   * @param[in] n first input integer 
   * @param[in] k second input integer 
   * @return (n+k)!/(n!k!)
   */   
  int combination(int n, int k);

  void rep(std::vector<std::vector<int> > &res, const std::vector<int> &values, std::vector<int> &item, unsigned int count);

  void variations(const std::vector<int> values, const int k, std::vector<std::vector<int> > &res);

  template < class T >
  void remove_numerical_zeros(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix, const double &threshold = ZERO)
  {
      matrix = (((matrix.array()) < threshold) &&
                ((matrix.array()) > -threshold) &&
                ((matrix.array()) != 0.0)).select(Eigen::MatrixXd::Constant(matrix.rows(), matrix.cols(), 0), matrix);
  }

  typedef double (*fun)(double);

  /**
   * @brief bisection_method implementation of the bisection method
   * @param[in] f monotonic function whose zero is to be found
   * @param[in] lb0 initial lower bound for root
   * @param[in] ub0 initial upper bound for root
   * @param[in] prec tolerance on root finding
   * @param[in] iter maximum number of iterations
   * @param[out] value of root
   * @return flag: 0 if method converged, 1 if max. number of iterations reached, -2 if initialization is wrong sign-wise, -1 if extremal values are not of opposite signs during iterations
   */  
  int bisection_method(fun f, const double &lb0, const double &ub0, const double &prec, const int &iter, double &root);

  int bisection_method_2(std::function<double(double)> f, const double &lb0, const double &ub0, const double &prec, const int &iter, double &root);

  /**
   * @brief Legendre evaluation of associated Legendre functions
   * @param[in] l order
   * @param[in] m degree
   * @param[in] x evaluation point
   * @return value of Plm at x
   */  
  double Legendre(int l, int m, double x);

  /**
   * @brief Legendre_derivative evaluation of derivative of associated Legendre functions
   * @param[in] l order
   * @param[in] m degree
   * @param[in] x evaluation point
   * @return value of d_Plm / d_x at x
   */  
  double Legendre_derivative(int l, int m, double x);

  /**
   * @brief Lagrange1d evaluation of Lagrange univariate polynomial
   * @param[in] times vector of interpolation points
   * @param[in] values vector of interpolated values
   * @param[in] t point where to evaluate polynomial
   * @return value at t
   */  
  double Lagrange1d(std::vector<double> times, std::vector<double> values, double t);

  /**
   * @brief LagrangeNd evaluation of Lagrange multivariate polynomial
   * @param[in] times vector of interpolation points
   * @param[in] values vector of interpolated vectors
   * @param[in] t point where to evaluate polynomial
   * @return vector of values at t
   */ 
  std::vector<double>  LagrangeNd(std::vector<double> times, std::vector<std::vector<double> > values, double t);

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

  template <class T>
  inline T quantile(const T* d, const unsigned int size, const double q)
  {
    if (size == 0) return T(0);
    if (size == 1) return d[0];
    if (q <= 0) return *std::min_element(d, d + size);
    if (q >= 1) return *std::max_element(d, d + size);

    double pos = (size - 1) * q;
    unsigned int ind = static_cast<unsigned int>(pos);
    double delta = pos - ind;
    std::vector<T> w(size); std::copy(d, d + size, w.begin());
    std::nth_element(w.begin(), w.begin() + ind, w.end());
    T i1 = *(w.begin() + ind);
    T i2 = *std::min_element(w.begin() + ind + 1, w.end());
    return i1 * (1.0 - delta) + i2 * delta;
  }

  template <class T>
  inline T quantile(const std::vector<T>& v, const double q)
    { return quantile(&v[0], v.size(), q); }



  /**
    * Function to compute coefficients of a polynomial given its roots
    *
    * @param   roots: Vector containing roots of polynomial
    *                 Ex. y = (x-x0)(x-x1) -> +x0,+x1 are input roots
    * @return coeffs: Coefficient of the polynomial written in explicit form
    *                 and normalized with a_n = 1.0, sorted from order n to
    *                 order 0.
    */
  std::vector<double> vieta_root2coef(const std::vector<double> &roots);

}

#endif // SMARTMATH_INLINEFUNCTIONS_H
