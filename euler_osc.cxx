/* solve x'' + x = 0 with forward and backward Euler methods
exact solution: cos(t) */

#include<fstream>
#include<cmath>
#include<string>
#include<sstream>

using namespace std;

int main(){
  
  const double t0 = 0.0;
  const double tend = 20.0*M_PI;
  const int dim = 2;
  
  const double dt = M_PI/10;
  ofstream out("forward_dt_pi_over_10.dat");
  ofstream out2("backward_dt_pi_over_10.dat");
  
  // const double dt = M_PI/100;
  // ofstream out("forward_dt_pi_over_100.dat");
  // ofstream out2("backward_dt_pi_over_100.dat");
  
  double yold[dim];
  double y[dim];

  // initial conditions
  const double x0 = 1.0;
  const double xp0 = 0.0;
    
  yold[0] = x0;
  yold[1] = xp0;

  out << t0 << "\t"  << yold[0] << "\t" << yold[1] << endl;

  // forward Euler
  for(double t = t0 + dt; t <= tend; t += dt)
    {
      y[0] = yold[0] + dt*yold[1];
      y[1] = yold[1] - dt*yold[0];
    
     out << t << "\t" << y[0] << "\t" << y[1] << endl;
      
      yold[0] = y[0];
      yold[1] = y[1];
    }

  
  // backward Euler
  yold[0] = x0;
  yold[1] = xp0;
  
  out2 << t0 << "\t"  << yold[0] << "\t" << yold[1] << endl;
  
  for(double t = t0 + dt; t <= tend; t += dt)
    {
      y[0] = (1.0/(1.0 + dt*dt))*(yold[0] + dt*yold[1]);
      y[1] = (1.0/(1.0 + dt*dt))*(yold[1] - dt*yold[0]);
    
     out2 << t << "\t" << y[0] << "\t" << y[1] << endl;
      
      yold[0] = y[0];
      yold[1] = y[1];
    }

  out.close();
  out2.close();

  return(0);
}
