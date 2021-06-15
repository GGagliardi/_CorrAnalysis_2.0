#include "../include/g_minus_2_utilities.h"


using namespace std;

const double m_muon= 0.10565837;
const double relerr= 1e-12;
const double Nresonances=8;


double kernel_K(double t) {

  auto F = [=](double x) -> double {

    if (x<eps(16)) return 0;
    return (4.0/pow(m_muon,2))*(1.0/sqrt(4.0 + pow(x,2)))*pow(  (sqrt(4.0+pow(x,2))-x)/(sqrt(4.0+pow(x,2))+x),2)*( (cos(m_muon*t*x)-1)/pow(x,2) + (1.0/2.0)*pow(t*m_muon,2));
  };


  double err;
  return boost::math::quadrature::gauss_kronrod<double,20>::integrate(F, 0.0, numeric_limits<double>::infinity() , 5, relerr, &err);


}

double Vdual(double t, double m_rho, double Edual, double Rdual) {


  return (5.0/(18.0*pow(M_PI,2)))*(Rdual/pow(t,3))*exp(-(m_rho+Edual)*t)*( 1.0 + (m_rho+Edual)*t + 0.5*pow( (m_rho+Edual)*t,2));

}


double Gamma_rpp(double omega, double g_rho_pipi, double Mpi) {

  double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));
  return pow(g_rho_pipi,2)*pow(k,3)/(6.0*M_PI*pow(omega,2));  
}

double Gamma_rpp_der(double omega, double g_rho_pipi, double Mpi) {

  double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));
  return 2*k*( (8*pow(Mpi,2) + pow(omega,2))/(4*pow(omega,3)))*pow(g_rho_pipi,2)/(6.0*M_PI)
  
    }

double h(double omega, double g_rho_pipi, double Mpi) {

  double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));
  return pow(g_rho_pipi,2)*pow(k,3)*(2.0/(omega*6.0*pow(M_PI,2)))*log( (omega + 2*k)/(2.0*Mpi));

}

double h_prime(double omega, double g_rho_pipi, double Mpi) {

  double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));
  return pow(g_rho_pipi,2)*pow(k,2)*(1.0/(6.0*pow(M_PI,2)*omega))*( 1.0 + (1.0+ 2.0*pow(Mpi/omega,2))*(omega/k)*log( (omega+2*k)/(2.0*Mpi)));

}



double A_pipi_0(double m_rho, double g_rho_pipi, double Mpi) {

  return h(m_rho,g_rho_pipi, Mpi) - (m_rho/2.0)*h_prime(m_rho, g_rho_pipi, Mpi) + pow(g_rho_pipi,2)*pow(Mpi,2)/(6.0*pow(M_PI,2));  

}


double Re_A_pipi(double omega, double m_rho, double g_rho_pipi, double Mpi) {

  return h(m_rho, g_rho_pipi, Mpi) + (pow(omega,2)- pow(m_rho,2))*(h_prime(m_rho, g_rho_pipi,Mpi)/(2.0*m_rho)) - h(omega, g_rho_pipi,Mpi);
}

double Im_A_pipi(double omega, double g_rho_pipi, double Mpi) {

  return omega*Gamma_rpp(omega,g_rho_pipi,Mpi);

}

Pfloat A_pipi(double omega, double m_rho, double g_rho_pipi, double Mpi) {

  return make_pair(Re_A_pipi(omega, m_rho, g_rho_pipi,Mpi), Im_A_pipi(omega, g_rho_pipi,Mpi));

}

double F_pi_GS_mod(double omega, double m_rho, double g_rho_pipi, double Mpi) {

  double Re = pow(m_rho,2) -pow(omega,2) - Re_A_pipi(omega,m_rho, g_rho_pipi,Mpi);
  double Im = -1.0*Im_A_pipi(omega, g_rho_pipi, Mpi);
  return fabs((pow(m_rho,2) - A_pipi_0(m_rho, g_rho_pipi,Mpi)))/( sqrt( pow(Re,2)+ pow(Im,2))); 

}

double cot_d_11(double k, double m_rho, double g_rho_pipi, double Mpi) {


  double omega= 2.0*sqrt( pow(Mpi,2) + pow(k,2));

  return  ( pow(m_rho,2) -pow(omega,2) - h(m_rho, g_rho_pipi,Mpi) -(pow(omega,2) - pow(m_rho,2))*(h_prime(m_rho,g_rho_pipi,Mpi)/(2.0*m_rho)) + h(omega, g_rho_pipi,Mpi))/Im_A_pipi(omega, g_rho_pipi, Mpi);

}

double d_11(double k, double m_rho, double g_rho_pipi, double Mpi) {


  return acot( cot_d_11(k, m_rho, g_rho_pipi, Mpi));

}


double cot_d_11_der(double k, double m_rho, double g_rho_pipi, double Mpi) {


  double omega = 2.0*sqrt( pow(Mpi,2) + pow(k,2));

  double d_omega_d_k = 4.0*k/omega;

  double num= ( pow(m_rho,2) -pow(omega,2) - h(m_rho, g_rho_pipi,Mpi) -(pow(omega,2) - pow(m_rho,2))*(h_prime(m_rho,g_rho_pipi,Mpi)/(2.0*m_rho)) + h(omega, g_rho_pipi,Mpi));

  double den = Im_A_pipi(omega, g_rho_pipi, Mpi);

  double num_der = (-2.0*omega -2.0*omega*(h_prime(m_rho, g_rho_pipi,Mpi)/(2.0*m_rho)) + h_prime(m_rho,g_rho_pipi,Mpi))*d_omega_d_k;

  
  double den_der = Gamma_rpp(omega, g_rho_pipi, Mpi) + omega*Gamma_rpp_der(omega, g_rho_pipi, Mpi);

  return d_omega_d_k*( num_der/den - num*den_der/pow(den,2));



}


double d_11_der(double k, double m_rho, double g_rho_pipi, double Mpi) {


  return -1.0*pow(sin(d_11(k,m_rho, g_rho_pipi,Mpi)),2)*cot_d_11_der(k, m_rho, g_rho_pipi, Mpi);

}

int degeneracy(int m) {

  int deg_val=0;

  //return number of lattice sites with distance m>0 from origin.
  assert(m >= 0);

  //assume x >= y >= z and moreover x, y ,z >= 0. We take into account degeneracy after
  int xmin= (int)(m/sqrt(3))  +1;
  int xmax= m;

  for(int x=xmin;x<=xmax;x++) {
    int dif = ipow(m,2) - ipow(x,2);
    int ymax;
    if(ceil((double)sqrt(dif)) == floor((double)sqrt(dif))) {
      //dif is a perfect square
      ymax= floor((double)sqrt(df));
    }
    else  ymax= (int)sqrt( pow(m,2) -pow(x,2));

    if( (m-x) % 2 == 0) {  int m_f = (m-x)/2; int prod = m_f*(m+x);
      if(ceil((double)sqrt(prod)) == floor((double)sqrt(prod))) ymin= floor((double)sqrt(prod));
      else ymin= (int)sqrt(prod) +1;
    }
    else ymin= (int)sqrt( (pow(m,2)-pow(x,2))/2.0) +1;

    
    int ymin;

    for(int y=ymin;y<=ymax;y++) {

      if(y<=x) {

	int diff = ipow(m,2) -ipow(x,2) - ipow(y,2);

	if (ceil((double)sqrt(diff)) == floor((double)sqrt(diff))) {
	  //number is a perfect square
	  //compute degeneracy
	  int z= floor((double)sqrt(diff));

	  if(z<=y) {
	    if(z==0 && y==0) deg_val += 6;
	    else if(z==0 && y!= 0) {
	      if(x==y) deg_val += 3*4;
	      else deg_val+=3*4*2;

	    }
	    else {  // x, y, z > 0
	      if(z==y && y==x) deg_val += 2*2*2;
	      else if(z==y && y != x) deg_val += 3*2*4;
	      else if(x==y && x != y) deg_val += 3*2*4;
	      else deg_val += 6*2*2*2;
	  
	    }
	  }
	}
      }
      
    }
  }

  
  return deg_val;

}

double Zeta_function_laplacian(double z) {

  double val= -M_PI -(1.0/sqrt(4*M_PI))*exp(pow(z,2))/pow(z,2);

  auto F1 = [&](double t) -> double { return (M_PI/2.0)*(exp(t*pow(z,2))-1.0)/(pow(t,3.0/2.0));};

  val += boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F1, 0, 1, 5, relerr);

  auto F2 = [&](int m) {

    double x=z;
    int y = m;
    auto F3 = [&](double t) -> double { return (M_PI/2.0)*exp(t*x*x -pow(M_PI,2)*y/t)/(pow(t,3.0/2.0));};

    double ret_val= boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F3, 0, 1, 5, relerr)  ;
    ret_val += (1.0/sqrt(4*M_PI))*exp(-(y-pow(x,2)))/(y-pow(x,2));
    return ret_val;
  };


  //compute sum m   nu(m)F2(m)

  double series_val=0;

  int max_m= 1000000;
  double new_val;

  for(int m=1;m<=max_m;m++) {

    new_val = degeneracy(m)*F3(m);
    if(new_val<= relerr*series_val) return series_val+val;
    else series_val += new_val;
  }

  crash("Unable to evaluate Zeta function laplacian to required accuracy. m_max: "+to_string(max_m));

  return 0.0;
}


double tan_phi(double z) {

  return -1.0*pow(M_PI, 3.0/2.0)*z/Zeta_function_laplacian(z);

}

double phi_der(double z) {

  auto f = [&](double t) { return phi(t);}

  double x=z;

  return boost::math::tools::finite_difference_derivative(f,x);


}

double Amplitude(double k, int L, double m_rho, double g_rho_pipi, double Mpi) {


  double omega= 2*sqrt( pow(k,2)+pow(Mpi,2));

  return ((2.0*pow(k,5)/(3.0*M_PI*pow(omega,2))))*pow(F_pi_GS_mod(omega, m_rho,g_rho_pipi,Mpi),2)/(k*d11_der(k, m_rho, g_rho_pipi,Mpi) + (k*L/(2.0*M_PI))*phi_der(k*L/(2.0*M_PI)));

}



double Find_pipi_energy_lev(int L, int n, double m_rho, double g_rho_pipi, double Mpi) {

  //find energy level n of pi-pi bound state

  assert(n>0);

  double Precision = 1e-9;
  double delta = 0.001;
  
 
  double a,b;

  double STEP_SEARCH= (1/10)*2.0*M_PI/L;

  auto F = [&](double k) -> double { return d_11(k, m_rho, g_rho_pipi, Mpi) + phi(k*L/(2.0*M_PI))-n*M_PI;};
  a= STEP_SEARCH;
  b=STEP_SEARCH/10.0;

  while (F(a)*F(b) > 0) {a  +=  STEP_SEARCH;}
 
 



  if(fabs(F(a)) < fabs(F(b))) {double atemp=a; a=b; b=atemp;}

  double c=a;
  bool FLAG = true;
  double s=b;
  double d=0;

  
  
  while(F(s) !=0 && fabs(b-a)>= Precision*fabs((double)(b+a)/2.0) ) {

  
    if((F(a) != F(c)) && (F(b) != F(c))) {//inverse quadratic interpolation
      s= a*F(b)*F(c)/((F(a)-F(b))*(F(a)-F(c))) + b*F(a)*F(c)/((F(b)-F(a))*(F(b)-F(c))) + c*F(a)*F(b)/((F(c)-F(a))*(F(c)-F(b)));
    }
    else s= b-(F(b)*(b-a)/(F(b)-F(a)));

    double s1= (double)(3*a+b/4.0);
    if( (s < s1 || s> b) || (FLAG==true && fabs(s-b) >= (double)fabs(b-c)/2) || (FLAG==false && fabs(s-b) >= (double)fabs(c-d)/2) || (FLAG==true && fabs(b-c) < delta) || (FLAG==false && fabs(c-d) < delta)) {
      
      FLAG=true;
      s= (a+b)/2.0;
      
    }
    
    else FLAG= false;

    d= c;
    c= b;
    if (F(a)*F(s)<0) b=s;
    else a=s;

    if(fabs(F(a)) < fabs(F(b))) {double atemp=a; a=b; b=atemp;}

   
  }

  
  return s;
  




}


double V_pipi(double t, int L, double m_rho, double g_rho_pipi, double Mpi) {

  double ret_val=0;

  for(int i_lev=1;i_lev<=Nresonances;i_lev++) {

    double k_n= Find_pipi_energy_lev(int L, int i_lev, double m_rho, double g_rho_pipi, double Mpi);

    double omega_n = 2.0*sqrt( pow(Mpi,2) + pow(k_n,2));
    ret_val += Amplitude(k_n, L, m_rho, g_rho_pipi, Mpi)*exp(-omega_n*t);

  }

  return ret_val;

}


double V_pipi_infL(double t, double m_rho_infL, double g_rho_pipi_infL, double Mpi_infL) {



  auto Integrand = [&](double omega) -> double {
    return (1.0/(48.0*pow(M_PI,2)))*pow(omega,2)*pow(1.0- pow(2.0*Mpi_infL/omega,2), 3.0/2.0)*exp(-omega*t)*F_pi_GS_mod(omega, m_rho_infL, g_rho_pipi_infL,Mpi_infl);};


  return boost::math::quadrature::gauss_kronrod<double, 15>::integrate(Integrand, 2*Mpi_infL, numeric_limits<double>::infinity(), 5, relerr)  ;


  
}

