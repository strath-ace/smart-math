#include "../../include/Utils/mixed_functions.h"

int smartmath::factorial(int n)
{
    if(n<0)
        smartmath_throw("FACT: factorial of non positive integer does not exist");

    return (n == 1 || n == 0) ? 1 : smartmath::factorial(n - 1) * n;
}

int smartmath::combination(int n, int k)
{
    int max = std::max(n, k);
    int min = std::min(n, k);
    int res = 1;
    int j = 1;

    while(j <= min)
    {
        res *= max + j;
        j++;
    }

    return res / smartmath::factorial(min);
}

int smartmath::bisection_method(fun f, const double &lb0, const double &ub0, const double &prec, const int &iter, double &root){

    if(lb0 > ub0)
        smartmath_throw("BISECTION_METHOD: lower bound must be smaller than upper one");
    if(prec <= 0.0)
        smartmath_throw("BISECTION_METHOD: required precision must be non-negative");    
    if(iter < 1)
        smartmath_throw("BISECTION_METHOD: maximum number of iterations must be non-negative");    

    double f_low = f(lb0);
    double f_up  = f(ub0);
    root = (ub0 + lb0) / 2.0;

    if( f_low * f_up > 0.0 )
        return -2;    
    else if( ub0 - lb0 <= prec )
        return 0;

    int i = 0;
    double ub = ub0, lb = lb0;
    double f_temp;
    while( (ub - lb > prec) && (i < iter) )
    {
        f_temp = f(root);

        if (f_temp == 0.0)
            return 0;

        if( f_temp * f_low > 0.0 )
        {
            lb = root;
            f_low = f_temp;
        }

        if( f_temp * f_up > 0.0 )
        {
            ub = root;
            f_up = f_temp;
        }

        if( f_low * f_up > 0.0 )
            return -1;

        root = (ub + lb) / 2.0;
        i++;
    }
     
    if(i == iter)
        return 1;
    else
        return 0;
}

double smartmath::Legendre(int l, int m, double x)
{
    if(x * x > 1.0)
        smartmath_throw("LEGENDRE: real number must be in [-1,1]");

    if(l < 0)
        return smartmath::Legendre(-l - 1, m, x);

    if(m < 0)
        return pow(-1.0, -m) * double(factorial(l + m)) * Legendre(l, -m, x) / double(factorial(l - m));

    if(m > l)
        return 0.0;        

    /* l>=m>=0 */
    double out = 1.0; // default value (l=0)

    if(l == 1)
    {
        if(m == 1)
            out = -sqrt(1.0 - x * x);
        if(m == 0)
            out = x;    
    }

    if(l > 1)
    {
        if(l == m)
            out = -double(2 * l - 1) * sqrt(1.0 - x * x)*smartmath::Legendre(l - 1, l - 1, x);
        else if(m == l - 1)
            out = double(2 * m + 1)* x * smartmath::Legendre(m, m, x);        
        else
            out = (double(2 * l - 1) * x * smartmath::Legendre(l - 1, m, x)-double(l - 1 + m)*smartmath::Legendre(l - 2, m, x)) / double(l - m);     
    }

    return out;
}

double smartmath::Legendre_derivative(int l, int m, double x)
{
    return (double(l) * x * smartmath::Legendre(l, m, x) - double(l + m) * smartmath::Legendre(l - 1, m, x)) / (x * x - 1.0); 
}
