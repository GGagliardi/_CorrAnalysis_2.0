#include "../include/ChPT_form_factors.h"


using namespace std;

using namespace std::complex_literals;

const double M_k = 0.493677;
const double M_Rho = 0.7754 ;
const double M_pi = 0.139570;
const double f_k =  0.155; //0.1136*sqrt(2); //0.155;
const double f_pi = 0.1304; //0.0932*sqrt(2);  //0.1304;
const double f_chir = f_k;
const double l10 = -0.0055;
const double l9 = 0.0069;
const double tol= 1.0e-12;

double H(double xk2, double m, double mu) {


  double H_val;
  
  
  double k = (1.0/(32.0*M_PI*M_PI))*(2.0*log( m/mu) + 1.0);
  
  auto F= [&](double x) {

	    complex<double> z{1.0- (xk2/pow(m,2))*x*(1.0-x),0.0};

	    return real(log(z));

	  };

    double J = -(1.0/(16.0*M_PI*M_PI))*boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F, 0.0, 1.0, 5, tol);

    complex<double> kr2{1.0 - 4.0*m*m/xk2,0.0};

    complex<double> kr = sqrt(kr2);

    complex<double> zdual = (1.0+kr)*(1.0+kr)/( (1.0-kr)*(1.0-kr)); 

 
    double JIII =  (1.0/(32.0*M_PI*M_PI))*(4.0 - real(kr*log(zdual)));

    /* if( J/JIII > 1 + 1e-6 || J/JIII < 1 + 1e-6) {

    cout.precision(15);
    
    cout<<"xk2: "<< xk2<<" xk2/m^2: "<<xk2/(m*m)<<" J: "<<J<<" JIII: "<<JIII<<endl;
    
    //crash("J and JIII in ChPT form factors do not coincide");

    } */
    

   H_val = (2.0/3.0)*l9 + (1.0/12.0)*(1.0 - 4.0*pow(m,2)/xk2)*J + (1.0/(288.0*M_PI*M_PI))  - k/6.0;


  return H_val;


}

void Compute_ChPT_form_factors( function<double(double,double)>& H1, function<double(double,double)>& H2, function<double(double,double)>& FA, function<double(double,double)>& FV) {

  
  auto def_FA = [&](double xk, double xq) -> double { return (8*M_k/f_chir)*(l9 + l10);   };
  auto def_FV = [&](double xk, double xq) -> double { return M_k/(4.0*(pow(M_PI,2))*f_chir);};
  auto def_H =  [&](double xk, double xq) -> double { return (2*M_k*f_k)*(    (2.0/(f_chir*f_chir))*H(pow(M_k*xk,2),M_pi, M_Rho) + 2.0*(2.0/(f_chir*f_chir))*H(pow(M_k*xk,2),M_k,M_Rho));};

  FA= def_FA;
  FV= def_FV;
  H1= def_H;
  H2= def_H;

  return;
}
