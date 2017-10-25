#include "../include/smartmath.h"

using namespace std;

int main(){

cout << "This is an example to illustrate symplectic integration with mixed variables." << endl;

/* Creating the dynamics */
smartmath::dynamics::spring<double> *dyn = new smartmath::dynamics::spring<double>();

/* Creating integrator */
smartmath::integrator::forest_mixedvar<double> prop(dyn);

/* Setting initial conditions */
std::vector<double> x(2), xf;
x[0] = 0.1; 
x[1] = 0.01; 

/* Integration */
double t_0 = 5.0; // initial time
double t_f = 10.0; // final time
prop.integrate(t_0, t_f, 1e2, x, xf); // integration with 100 steps

cout << "The analytical solution to the harmonic oscillator is compared to the numerical one:" << endl;

cout << sqrt(x[0] * x[0] + x[1] * x[1]) * sin(t_f - (t_0 - atan2(x[0], x[1]))) << " VS " << xf[0] << endl;

delete dyn;
}
