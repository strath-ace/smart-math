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
#include <time.h>
#include <functional>
#include <utility>

namespace smartmath
{

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
      if(fabs(x) <= ZERO)
      {
          std::cout<<"ERROR: Division by zero."<<std::endl;
          throw std::exception();
      }
      return 1.0 / x;
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

}

#endif // SMARTMATH_INLINEFUNCTIONS_H
