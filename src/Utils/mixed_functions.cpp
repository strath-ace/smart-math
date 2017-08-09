#include "../../include/Utils/mixed_functions.h"

// using namespace std;
// using namespace smartmath;

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


// template < class T >
// void smartmath::remove_numerical_zeros(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix, const double &threshold)
// {
//     matrix = (((matrix.array()) < threshold) &&
//               ((matrix.array()) > -threshold) &&
//               ((matrix.array()) != 0.0)).select(Eigen::MatrixXd::Constant(matrix.rows(), matrix.cols(), 0), matrix);
// }


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
        return -1;    
    else if( ub0 - lb0 <= prec )
        return 0;

    int i = 0;
    double ub = ub0, lb = lb0;
    double f_temp;
    while( (ub - lb > prec) && (i < iter) )
    {
        f_temp = f(root);

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

int smartmath::bisection_method_2(
        std::function<double(double)> f,
        const double &lb0,
        const double &ub0,
        const double &prec,
        const int &iter,
        double &root) {
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
        return -1;
    else if( ub0 - lb0 <= prec )
        return 0;

    int i = 0;
    double ub = ub0, lb = lb0;
    double f_temp;
    while( (ub - lb > prec) && (i < iter) )
    {
        f_temp = f(root);

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

double smartmath::Lagrange1d(std::vector<double> times, std::vector<double> values, double t)
{
    if(times.size()!=values.size())
        smartmath_throw("LAGRANGE1D: number of values must be equal to number of interpolation points");

    double output = 0.0;
    double prod1, prod2;

    for(unsigned int i=0;i<times.size();i++){
        prod1 = 1.0;
        prod2 = 1.0;  
        for(unsigned int j=0;j<times.size();j++){  
            if(i!=j){
                if(times[j]==times[i])
                    smartmath_throw("LAGRANGE1D: interpolated points must be different");
                prod1 *= t-times[j];
                prod2 *= times[i]-times[j];
            }
        }
        output += prod1*values[i]/prod2;
    }

    return output; 
}

std::vector<double>  smartmath::LagrangeNd(std::vector<double> times, std::vector<std::vector<double> > values, double t){

    if(times.size()!=values.size())
        smartmath_throw("LAGRANGEND: number of values must be equal to number of interpolation points");

    for(unsigned int k=0;k<values.size();k++){
        if(values[k].size()!=values[0].size())
            smartmath_throw("LAGRANGEND: function evaluations at interpolation points must have same number of components");
    }

    std::vector<double> outputs(values[0].size(),0.0);
    double prod1, prod2;

    for(unsigned int i=0;i<times.size();i++){
        prod1 = 1.0;
        prod2 = 1.0;  
        for(unsigned int j=0;j<times.size();j++){  
            if(i!=j){
                if(times[j]==times[i])
                   smartmath_throw("LAGRANGEND: interpolated points must be different");
                prod1 *= t-times[j];
                prod2 *= times[i]-times[j];
            }
        }
        for(unsigned int k=0;k<values[0].size();k++)
            outputs[k] += prod1*values[i][k]/prod2;

    }
    return outputs;
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

    if(l==1)
    {
        if(m==1)
            out=-sqrt(1.0-x*x);
        if(m==0)
            out=x;    
    }

    if(l>1)
    {
        if(l==m)
            out=-double(2*l-1)*sqrt(1.0-x*x)*smartmath::Legendre(l-1,l-1,x);
        else if(m==l-1)
            out=double(2*m+1)*x*smartmath::Legendre(m,m,x);        
        else
            out=(double(2*l-1)*x*smartmath::Legendre(l-1,m,x)-double(l-1+m)*smartmath::Legendre(l-2,m,x))/double(l-m);     
    }

    return out;
}


double smartmath::Legendre_derivative(int l, int m, double x)
{
    return (double(l)*x*smartmath::Legendre(l,m,x)-double(l+m)*smartmath::Legendre(l-1,m,x))/(x*x-1.0); 
}

Eigen::MatrixXd smartmath::sample_multivariate_normal_distribution(const Eigen::VectorXd &mean,
                                                        const Eigen::MatrixXd &covar,
                                                        const int &N_samples)
{
    Eigen::EigenMultivariateNormal<double> sampler(mean,covar,false,time(NULL));
    return sampler.samples(N_samples);
}

Eigen::MatrixXd smartmath::sample_truncated_multivariate_normal_distribution(const Eigen::VectorXd &lower_bounds,
                                                         const Eigen::VectorXd &upper_bounds,
                                                         const Eigen::VectorXd &mean,
                                                         const Eigen::MatrixXd &covar,
                                                         const unsigned int &N_samples,
                                                            double &pr_valid_samples)
{
    Eigen::EigenMultivariateNormal<double> sampler(mean,
                                                   covar,
                                                   false,
                                                   time(NULL));
    Eigen::MatrixXd samples = sampler.samples_truncated(
                                lower_bounds,
                                upper_bounds,
                                N_samples,
                                pr_valid_samples);

    return samples;
}

Eigen::MatrixXd smartmath::sample_truncated_multivariate_normal_distribution(const Eigen::VectorXd &mean,
                                                                  const Eigen::MatrixXd &covar,
                                                                  const double &min_pr_valid_samples,
                                                                  const unsigned int &N_samples,
                                                                  double &pr_valid_samples)
{
    int d = mean.size(); //Number of dimensions

    //Step 1: Diagonalize the covariance matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(covar);
    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
    Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();

    //Step 2: Calculate the range for acceptance in the eigenspace, using the chebyshev inequality
    // for independent variables
    double k = smartmath::chebyshev_inequality::get_k_for_independent_variables(
                d,
                min_pr_valid_samples);

    //Step 3: Transform mean and covariance to the eigenspace
    Eigen::VectorXd mean_eigenspace = (eigenvectors.transpose())*mean;
    Eigen::MatrixXd covar_eigenspace = eigenvalues.asDiagonal();

    // Step 4: Get the lowe and upper bounds in the eigen space
    Eigen::VectorXd sd_eigenspace = eigenvalues.array().sqrt();
    Eigen::VectorXd lower_bounds = mean_eigenspace - k*sd_eigenspace;
    Eigen::VectorXd upper_bounds = mean_eigenspace + k*sd_eigenspace;

    //Step 5: Truncated Sampling in the eigenspace
    Eigen::EigenMultivariateNormal<double> sampler(mean_eigenspace,
                                                   covar_eigenspace,
                                                   false,
                                                   time(NULL));
    Eigen::MatrixXd samples_eigenspace = sampler.samples_truncated(lower_bounds,
                              upper_bounds,
                              N_samples,
                              pr_valid_samples);

    //Step 6: Rotate the samples back to the original space
    Eigen::MatrixXd samples = (eigenvectors)*samples_eigenspace;

    return samples;
}

int smartmath::find_PC_bounds_normal_distribution(const Eigen::VectorXd &mean,
                                       const Eigen::MatrixXd &covar,
                                       const double &min_pr_valid_samples,
                                       Eigen::VectorXd &lower_bounds,
                                       Eigen::VectorXd &upper_bounds)
{
    int d = mean.size(); //Number of dimensions

    //Step 1: Diagonalize the covariance matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(covar);
    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
    Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();

    //Step 2: Calculate the range for acceptance in the eigenspace, using the chebyshev inequality
    // for independent variables
    double k = smartmath::chebyshev_inequality::get_k_for_independent_variables(
                d,
                min_pr_valid_samples);

    //Step 3: Transform mean and covariance to the eigenspace
    Eigen::VectorXd mean_eigenspace = (eigenvectors.transpose())*mean;

    // Step 4: Get the lowe and upper bounds in the eigen space
    Eigen::VectorXd sd_eigenspace = eigenvalues.array().sqrt();
    lower_bounds = mean_eigenspace - k*sd_eigenspace;
    upper_bounds = mean_eigenspace + k*sd_eigenspace;

    return 0;
}

/**
  * Function to compute coefficients of a polynomial given its roots
  *
  * @param   roots: Vector containing roots of polynomial
  *                 Ex. y = (x-x0)(x-x1) -> +x0,+x1 are input roots
  * @return coeffs: Coefficient of the polynomial written in explicit form
  *                 and normalized with a_n = 1.0, sorted from order n to
  *                 order 0.
  */
std::vector<double> smartmath::vieta_root2coef(const std::vector<double> &roots)
{
    if ( roots.empty() || roots.size() == 0 )
        smartmath_throw("vieta_root2coef: Input roots vector has no element");

    // Polynomial degree
    unsigned int deg = roots.size() ;

    // Initialize vector of coefficients to 1.0 (scale vector for a_n)
    std::vector<double> coeffs(deg+1,1.0), coeffs2, roots2;

    if(deg == 1) // linear case
    {
        coeffs[1] = -roots[0];
    }
    else
    {
        for(unsigned int i = 0; i < deg - 1; i++)
            roots2.push_back(roots[i]);

        coeffs2 = vieta_root2coef(roots2);

        /* first and last */
        coeffs[0] = 1.0;
        coeffs[deg] = - roots[deg-1] * coeffs2[deg-1];

        /* recursive formula from 2 to deg */
        for(unsigned int n = 1; n < deg; n++)
            coeffs[n] = coeffs2[n] - roots[deg-1] * coeffs2[n-1];

    }

    return coeffs;
}
