#include "../../include/Utils/chebyshev_inequality.h"

using namespace smartmath;

double chebyshev_inequality::get_k_for_independent_variables(const int &N,
                                                             const double &pr_threshold)
{
    //Sanity check
    if (N <= 0)
        smartmath_throw("N must be a positive integer");
    if (pr_threshold < 0 || pr_threshold > 1)
        smartmath_throw("pr_threshold must be a value in [0,1]");

    return 1/sqrt(1-exp(log(pr_threshold)/N));
}
