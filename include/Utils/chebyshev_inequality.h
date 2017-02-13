#ifndef SMARTMATH_CHEBYSHEV_INEQUALITY_H
#define SMARTMATH_CHEBYSHEV_INEQUALITY_H

#include <cmath>
#include "../exception.h"

namespace smartmath
{
     /**
     * @brief The chebyshev_inequality class contains methods and utils to deal with the
     * chebyshev inequality, which guarantees that, for a wide class of probability
     * distributions, "nearly all" values are close to the mean—the precise statement
     * being that no more than 1/k^2 of the distribution's values can be more than k
     * standard deviations away from the mean (or equivalently, at least 1−1/k2 of the
     * distribution's values are within k standard deviations of the mean)
     * @author Victor Rodriguez
     */
    class chebyshev_inequality
    {
    public:
        /**
         * @brief get_k_for_independent_variables Special case of the chebyshev inequality
         * for the multivariate case, in which we consider that the variables are independent
         * and also, that we want to obtain the same k for each of the variables.
         * @param Number of variables
         * @param pr_threshold Threshold value for the inequality. The probability of finding
         * values in the distribution within k standard deviations from the mean is greater or equal
         * than this value, by the chevyshev inequality.
         * @return
         */
        static double get_k_for_independent_variables(const int &N, const double &pr_threshold);
    };
}

#endif // SMARTMATH_CHEBYSHEV_INEQUALITY_H
