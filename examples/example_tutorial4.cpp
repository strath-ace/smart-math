#include "../include/smartmath.h"

using namespace std;

const double t_0 = 0.0;

std::vector<int> event_test(std::vector<double> x, double t){

std::vector<int> v(1, 0); // default returned value is 0. Several events can be modelled by adding more components to the vector v.

if((x[0] < 0.0) && (t >= t_0))
	v[0] = 1; // if the event is detected, returned value is 1

return v;
}

int main(){

cout << "This example shows you how to use events." << endl;

/* Creating the dynamics */
smartmath::dynamics::pendulum<double> *dyn = new smartmath::dynamics::pendulum<double>();

/* Creating integrators */
double tolerance = 1.0e-8; // integration tolerance for stepsize control
double max_factor = 2.0; // maximum multiplying factor used with stepsize control
double min_step = 1.0e-6; // minimum stepsize 
double max_step = 1.0; // maximum stepsize 
smartmath::integrator::rkf45<double> prop(dyn, tolerance, max_factor, min_step, max_step);

/* Setting initial conditions */

std::vector<double> x(2), xf;
x[0] = 0.1; // angle
x[1] = 0.01; // angular velocity

cout << "The propagation of a pendulum is set until the angle becomes negative" << endl;

/* Integration */
double t_f = 10.0; // projected final time
prop.integrate(t_0, t_f, 100, x, xf, event_test); // the variable t_f is erased and set to the termination time if an event is detected before its original value

delete dyn;

}
