#include "../../include/Utils/mixed_functions.h"

using namespace std;
using namespace smartmath;

int smartmath::factorial(int n)
{
    if(n<0)
        smartmath_throw("FACT: factorial of non positive integer does not exist");

    return (n == 1 || n == 0) ? 1 : smartmath::factorial(n - 1) * n;
}

int smartmath::combination(int n, int k)
{
    int max = std::max(n,k);
    int min = std::min(n,k);
    int res = 1;
    int j = 1;

    while(j<=min){
        res *= max+j;
        j++;
    }

    return res/smartmath::factorial(min);
}


void smartmath::rep(std::vector<std::vector<int> > &res, const std::vector<int> &values, std::vector<int> &item, unsigned int count){
    if (count < item.size()){
        for (unsigned int i = 0; i < values.size(); i++) {
            item[count] = values[i];
            unsigned int tmp_count = count + 1;
            smartmath::rep(res, values, item, tmp_count);
        }
    }else{
        res.push_back(item);
    }
}


void smartmath::variations(const std::vector<int> values, const int k, std::vector<std::vector<int> > &res){
    res.clear();

    std::vector<int> item(k);
    smartmath::rep(res, values, item, 0);
}

double smartmath::bisection_method(fun f, double lb, double ub, double prec){
    double f_low = f(lb);
    double f_up  = f(ub);
    if( (f_low*f_up) < 0 && (ub-lb) <= prec )
    {
        return lb;
    }

    double temp   = (ub+lb)/2;
    double f_temp = f(temp);

    if( (f_temp*f_low) > 0 )
        lb = temp;
    
    if( (f_temp*f_up) > 0 )
        ub = temp;
     
return smartmath::bisection_method(f, lb, ub, prec);
}

double smartmath::Legendre(int l, int m, double x)
{
    if(x*x>1.0)
        smartmath_throw("LEGENDRE: real number must be in [-1,1]");

    if(l<0)
        return smartmath::Legendre(-l-1,m,x);

    if(m<0)
        return pow(-1.0,-m)*double(factorial(l+m))*Legendre(l,-m,x)/double(factorial(l-m));

    if(m>l)
        return 0.0;        

    /* l>=m>=0 */
    double out=1.0; // default value (l=0)

    if(l==1){
        if(m==1){
            out=-sqrt(1.0-x*x);
        }
        if(m==0){
            out=x;
        }       
    }

    if(l>1){
        if(l==m){
            out=-double(2*l-1)*sqrt(1.0-x*x)*smartmath::Legendre(l-1,l-1,x);
        }
        else if(m==l-1){
            out=double(2*m+1)*x*smartmath::Legendre(m,m,x);
        }        
        else{
            out=(double(2*l-1)*x*smartmath::Legendre(l-1,m,x)-double(l-1+m)*smartmath::Legendre(l-2,m,x))/double(l-m); 
        }        
    }

    return out;
}


double smartmath::Legendre_derivative(int l, int m, double x)
{
    if(x*x>=1.0)
        smartmath_throw("LEGENDRE: real number must be in ]-1,1[");  

    return (double(l)*x*smartmath::Legendre(l,m,x)-double(l+m)*smartmath::Legendre(l-1,m,x))/(x*x-1.0); 
}

int smartmath::sample_truncated_normal_distribution(const double &lower_bound,
                                                         const double &upper_bound,
                                                         const double &mean,
                                                         const double &sd,
                                                         const unsigned int &N_samples,
                                                         std::vector<double> &result)
{
    //Sanity checks
    if (lower_bound >= upper_bound)
        smartmath_throw("Lower bound must be less than Upper bound");
    if (N_samples < 0)
        smartmath_throw("N_samples must be positive");
    if (result.size() != N_samples)
        smartmath_throw("The size of the result vector must be N samples");

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(mean,sd);

    //While loop until reaching the desired number of samples
    unsigned int valid_samples_counter = 0;
    while (valid_samples_counter <= N_samples)
    {
        double sample = distribution(generator);

        if (sample >= lower_bound && sample <= upper_bound)
        {
            result[valid_samples_counter] = sample;
            valid_samples_counter++;
        }
    }

    return 0;
}
