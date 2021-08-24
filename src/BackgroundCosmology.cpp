#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{
  // Constants
  const double pi = 3.14159265;

  // TODO: Compute OmegaR, OmegaNu, OmegaLambda, H0
  this->H0 = h * Constants.H0_over_h;                   // Hubble parameter
  this->OmegaR = 2 * pow(pi, 2) / 30 * pow(Constants.k_b * TCMB, 4) / (pow(Constants.hbar, 3) * pow(Constants.c, 5)) \
                 * 8 * pi * Constants.G / (3 * pow(H0, 2));      // Present radiation density parameter
  this->OmegaNu = Neff * 7 / 8 * pow(4 / 11, 4 / 3) * OmegaR;                         // Present nutrino density parameter
  this->OmegaLambda = 1 - (OmegaK + OmegaB + OmegaCDM + OmegaR + OmegaNu);            // Present cosmological constant density parameter
}

// Do all the solving. Compute eta(x) and t(x)

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
  Utils::StartTiming("Time");

  // TODO: Set the range of x and the number of points for the splines
  double x_start = -20;
  double x_end = 0;
  double n = 100;
  Vector x_array = Utils::linspace(x_start, x_end, n);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx) {

    // TODO: Set the rhs of the detadx ODE
    detadx[0] = Constants.c / (Hp_of_x(x));

    return GSL_SUCCESS;
  };

  // The ODE for dtdx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx) {
    // TODO: Set the rhs of the dtdx ODE
    dtdx[0] = 1 / H_of_x(x);

    return GSL_SUCCESS;
  };

  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline and t_of_x_spline
  
  // Set Initial condition
  double eta_0 = Constants.c / (exp(x_start) * Hp_of_x(x_start));
  Vector eta_i{eta_0};
  double t_0 = 1 / (2 * H_of_x(x_start));
  Vector t_i{t_0};

  // Solve eta ODE
  ODESolver ode;
  ode.solve(detadx, x_array, eta_i);
  Vector y_array = ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, y_array, "eta");

  Utils::EndTiming("Eta");

  // Solve t ODE
  ode.solve(dtdx, x_array, t_i);
  y_array = ode.get_data_by_component(0);

  t_of_x_spline.create(x_array, y_array, "t");

  Utils::EndTiming("Time");
}

// Get methods

double BackgroundCosmology::H_of_x(double x) const{
  // TODO: Implement...
  double H = H0 * sqrt((OmegaR + OmegaNu) * exp(-4 * x) + (OmegaB + OmegaCDM) * exp(-3 * x) + OmegaK * exp(-2 * x) + OmegaLambda);
  
  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  // TODO: Implement...
  double Hp = exp(x) * H_of_x(x);

  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  // TODO: Implement...
  double dHpdx = exp(x) * (2 * pow(H_of_x(x), 2) + (-4 * (OmegaR + OmegaNu) * exp(-4 * x) - 3 * (OmegaB + OmegaCDM) * exp(-3 * x) - 2 * OmegaK * exp(-2 * x)) / (2 * H_of_x(x)));

  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  // TODO: Implement...
  double A = (OmegaR + OmegaNu);
  double B = (OmegaB + OmegaCDM);
  double C = OmegaK;
  double ABC = 2 * pow(H0, 2) * (-4 * A * exp(-4 * x) - 3 * B * exp(-3 * x) - 2 * C * exp(-2 * x));
  double fraction = (ABC + 16 * A * exp(-4 * x) + 9 * B * exp(-3 * x) + 4 * C * exp(-2 * x)) * (2 * pow(H_of_x(x), 2) - ABC) / (2 * pow(H_of_x(x), 3));
  double ddHddx = exp(x) * (dHpdx_of_x(x) + fraction);

  return ddHddx;
}

double BackgroundCosmology::get_OmegaB(double x) const{
  // TODO: Implement...
  if (x == 0.0) {return OmegaB;}
  else {return OmegaB * pow(H0, 2) / (exp(3 * x) * pow(H_of_x(x), 2));}
}

double BackgroundCosmology::get_OmegaR(double x) const{
  // TODO: Implement...
  if (x == 0.0) {return OmegaR;}
  else {return OmegaR * pow(H0, 2) / (exp(4 * x) * pow(H_of_x(x), 2));}
}

double BackgroundCosmology::get_OmegaNu(double x) const{
  // TODO: Implement...
  if (x == 0.0) {return OmegaNu;}
  else {return OmegaNu * pow(H0, 2) / (exp(4 * x) * pow(H_of_x(x), 2));}
}

double BackgroundCosmology::get_OmegaCDM(double x) const{
  // TODO: Implement...
  if (x == 0.0) {return OmegaCDM;}
  else {return OmegaCDM * pow(H0, 2) / (exp(3 * x) * pow(H_of_x(x), 2));}
}

double BackgroundCosmology::get_OmegaLambda(double x) const{
  // TODO: Implement...
  if (x == 0.0) {return OmegaLambda;}
  else {return OmegaLambda * pow(H0, 2) / (pow(H_of_x(x), 2));}
}

double BackgroundCosmology::get_OmegaK(double x) const{
  // TODO: Implement...
  if (x == 0.0) {return OmegaK;}
  else {return OmegaK * pow(H0, 2) / (exp(2 * x) * pow(H_of_x(x), 2));}
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << "Time:        " << t_of_x(0) / (60 * 60 * 24 * 365.25 * 1e9) << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
double x_start = -20;
double x_end = 0;
double n = 1000;
Vector x_array = Utils::linspace(x_start, x_end, n);

void BackgroundCosmology::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << H_of_x(x)          << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp << t_of_x(x)          << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

