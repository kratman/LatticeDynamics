#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;

//A program to calculate the Cv based on the Debye model

//Constants
const double hbar = 1.054571726e-34;
const double pi = 4.0*atan(1.0);
const double k = 1.3806488e-23;
const double R = 8.3145;

int main()
{
  //Set up variables
  double Cv = 0, T = 1, kappa = 0, tauU, tauV;
  double theta = 1, tau = 1, v = 1, vp = 1;
  double Xd, x, y, f = 0, Na = 1, vs = 1, Nm = 1;
  double alpha = 1, mu = 1, comp = 1, V = 1;
  fstream inp;
  string filename = "data.txt";
  inp.open(filename.c_str(), ios_base::in);
  string dummy;
  //Loop over data range
  int ct = 0;
  double dx = 0.0000001;
  inp >> Na >> f; //Note: f = N_defect/N_mol^2
  cout << "0.0  0.0  0.0" << '\n';
  while (ct <= 19)
  {
    //Read Input
    ct += 1;
    inp >> T >> V >> vp >> vs >> dummy;
    inp >> mu >> dummy >> comp;
    inp >> alpha >> theta >> Nm;
    V *= 1e-10;
    V *= V*V;
    vp *= 1000;
    vs *= 1000;
    v = (1/(3*vp*vp*vp))+(2/(3*vs*vs*vs));
    v = pow(v,(-1.0/3.0));
    mu *= 1e9;
    comp /= 1e9;
    Cv = 0;
    x = dx;
    Xd = theta/T;
    kappa = 0;
    //Integrate Cv
    while (x <= Xd)
    {
      y = exp(x);
      Cv += dx*x*x*x*x*y/((y-1)*(y-1));
      x += dx;
    }
    Cv *= 3*Nm*k*T*T*T/(theta*theta*theta);
    cout << T << " " << R*Cv/k;
    //Integrate thermal conductivity
    kappa = 0;
    x = dx;
    while (x <= Xd)
    {
      tauU = 2*alpha*alpha*Na*V*k*k*T*T*T*x*x;
      tauU /= (Cv*Cv*mu*comp*comp*theta*hbar);
      tauV = 9*V*f*k*k*k*k*T*T*T*T*x*x*x*x;
      tauV /= (4*pi*hbar*hbar*hbar*hbar*v*v*v);
      tau = 1/(tauU+tauV);
      tau *= k*k*k*k*T*T*T;
      tau /= hbar*hbar*hbar*2*pi*pi*v;
      y = exp(x);
      kappa += tau*dx*x*x*x*x*y/((y-1)*(y-1));
      x += dx;
    }
    cout << " " << kappa << '\n';
  }
  inp.close();
  return 0;
}

