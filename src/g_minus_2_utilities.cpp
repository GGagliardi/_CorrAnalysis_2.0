#include "../include/g_minus_2_utilities.h"


using namespace std;

const double m_muon= 0.10565837;
const double relerr= 1e-12;
const double Nresonances=4;



bool Is_perfect_square(int x) {

  int sqrt_x = sqrt(x);
  if (sqrt_x*sqrt_x == x) return true;
  else return false;
}


double kernel_K(double t, double MV, string channel) {

  double m_muon_MV;
  if(channel=="l") {
    m_muon_MV= m_muon;
  }
  else if(channel == "s") m_muon_MV= m_muon/1.0195;
  else if(channel == "c") m_muon_MV= m_muon/3.0969;
  else crash("In kernel_K channel: "+channel+" is not implemented");

  auto F = [=](double x) -> double {

    if (x<eps(16)) return 0;
    return (4.0/pow(m_muon_MV*MV,2))*(1.0/sqrt(4.0 + pow(x,2)))*pow(  (sqrt(4.0+pow(x,2))-x)/(sqrt(4.0+pow(x,2))+x),2)*( (cos(m_muon_MV*MV*t*x)-1)/pow(x,2) + (1.0/2.0)*pow(t*m_muon_MV*MV,2));
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
  return 2*k*( (8*pow(Mpi,2) + pow(omega,2))/(4*pow(omega,3)))*pow(g_rho_pipi,2)/(6.0*M_PI);
  
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


  return atan( 1.0/cot_d_11(k, m_rho, g_rho_pipi, Mpi));

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

  if(m==0) return 1;
  int deg_val=0;

  //return number of lattice sites with distance m>0 from origin.
  assert(m > 0);

  int x2_max= m;
  int x2_min;
  if(m%3==0) x2_min=m/3;
  else x2_min= floor( m/3 + 1) ;
  

  

  for(int x2=x2_min;x2<=x2_max;x2++) {
    if(Is_perfect_square(x2)) {  //x2 is a perfect square
      
      int y2_max= m-x2;
      int y2_min;
      if( (m-x2)%2 == 0) y2_min= (m-x2)/2;
      else y2_min = floor( (m-x2)/2 +1);
      
      
      for(int y2=y2_min;y2<=y2_max;y2++) {
	
	if(y2<=x2 && Is_perfect_square(y2))  { //y2 is a perfect square
	
	  int z2 = m -x2 - y2;
	
	  if (Is_perfect_square(z2) && z2<=y2) { // z2 is a perfect square
	      if(z2==0 && y2==0) deg_val += 6;
	      else if(z2==0 && y2!= 0) {
		if(x2==y2) deg_val += 3*4;
		else deg_val+=3*4*2;
	      
	      }
	      else {  // x, y, z > 0
		if(z2==y2 && y2==x2) deg_val += 2*2*2;
		else if(z2==y2 && y2 != x2) deg_val += 3*2*4;
		else if(x2==y2 && z2 != y2) deg_val += 3*2*4;
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

  assert(z> 0);

  double val= -M_PI -(1.0/sqrt(4.0*M_PI))*exp(pow(z,2))/pow(z,2);

  cout<<"z: "<<z<<endl;
  cout<<endl<<"val0: "<<val<<endl;

  auto F1 = [&](double t) -> double {
	      if(t<eps(10)) return 0.0;
	      else return (M_PI/2.0)*(exp(t*pow(z,2))-1.0)/(pow(t,3.0/2.0));
	    };

  val += boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F1, 0, 1, 5, relerr);

  cout<<endl<<"val1: "<<val<<endl;

  auto F2 = [&](int m) {

    double x=z;
    int y = m;
    auto F3 = [&](double t) -> double {
		if(t<eps(10)) return 0.0;
		else return (M_PI/2.0)*exp(t*x*x -pow(M_PI,2)*((double)y)/t)/(pow(t,3.0/2.0));
	      };

    double ret_val= boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F3, 0, 1, 5, relerr)  ;
    ret_val += (1.0/sqrt(4*M_PI))*exp(-(y-pow(x,2)))/(y-pow(x,2));
    return ret_val;
  };


  //compute sum m   nu(m)F2(m)

  double relerr_laplacian= 1.0e-7;
  double series_val=0;

  int max_m= 1000;
  double new_val;

  for(int m=1;m<=max_m;m++) {

    new_val = degeneracy(m)*F2(m);
    if(new_val<= relerr_laplacian*series_val) return series_val+val;
    else series_val += new_val;
  }

  crash("Unable to evaluate Zeta function laplacian to required accuracy. m_max: "+to_string(max_m)+", Series_val: "+to_string_with_precision(series_val,6)+" , last term: "+to_string_with_precision(new_val,6));

  return 0.0;
}

double Zeta_function_laplacian_Luscher(double z) {


  
  double inf = 1e-20;

  assert(z>= 0) ;

  //compute Zeta_function using Luscher parametrization


  //l2 is cutoff
  int l2= (int)pow(z,2) + 1;
  int thresh = l2+ 10;

  double t0= pow(z,2)>2.0?M_PI/pow(z,2):M_PI/2.0;
  

  //sum first terms
  double val=0;
  for(int m=0; m < l2; m++) val += (1.0/sqrt(4*M_PI))*degeneracy(m)*(1.0/(m-pow(z,2)));

  auto F = [&](double t) -> double {

	     bool mode = t>t0;
	     double heat_kernel_red=0;
	     if(mode) {
	       for(int m=l2;m<=thresh; m++)

		 heat_kernel_red += (1.0/pow(2.0*M_PI,3))*degeneracy(m)*exp(-t*m);
	     }
	     else {
	       for(int m=0; m<= thresh;m++) {
		 if(m<l2) {
		   heat_kernel_red += degeneracy(m)*(( 1.0/pow(4.0*M_PI*t,3.0/2.0))*exp(-(pow(M_PI,2)/t)*m)-(1.0/pow(2.0*M_PI,3))*exp(-t*m)  );
		   
		 }
		 else {
		   heat_kernel_red += degeneracy(m)*( 1.0/pow(4.0*M_PI*t,3.0/2.0))*exp(-(pow(M_PI,2)/t)*m);
		 }

	       }
	       
	     }
	     
	     double res= pow(2.0*M_PI,3)*( exp(t*pow(z,2))*(1.0/sqrt(4*M_PI))*heat_kernel_red - 1.0/(pow(4.0*M_PI,2)*pow(t,3.0/2.0)));

	     
	     if( isnan(res)) crash("In integral repr. Zeta function t: "+to_string_with_precision(t,3)+" m: "+to_string_with_precision(pow(z,2), 3)+" , integrand is nan");
	    
	     return res;
	   };


  double int_val= boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F, inf, 20.0 , 5, relerr);

  if( isnan( int_val)) crash("integral in Generalized zeta is nan");
  return val + int_val;
}


double tan_phi(double z) {

  int n2= (int)pow(z,2); //check whether is integer.
  int n2_p= n2+1;
  if( pow(z,2) < 1e-8) return 0.0;
  if( fabs( pow(z,2) - n2) < 1e-8 ||  fabs(pow(z,2) - n2_p) < 1e-8 ) return 0.0;

  return -1.0*pow(M_PI, 3.0/2.0)*z/Zeta_function_laplacian_Luscher(z);

}

double phi_der(double z) {

  auto tf = [&](double t) { return tan_phi(t);};

  double x=z;
  double err;

  double result= boost::math::differentiation::finite_difference_derivative(tf,x, &err);

  if(err/result > 1.0e-6) crash("boost::math::differentiation::finite_difference_derivative has reached a low level of accuracy in phi_der(double). rel_err: "+to_string_with_precision(err/result,4));

  return result/(1.0+ pow(tan_phi(z),2));

}

double phi(double z) {

  int n2= (int)pow(z,2); //check whether is integer.
  int n2_p= n2+1;
  if( pow(z,2) < 1e-8) return 0.0;
  if( fabs( pow(z,2) - n2) < 1e-8 ||  fabs(pow(z,2) - n2_p) < 1e-8 ) return 0.0;
  double sum= Zeta_function_laplacian_Luscher(z);
  if( fabs(sum) < 1e-8) return 0.0;

  return atan( -1.0*pow(M_PI, 3.0/2.0)*z/sum);

}

double Amplitude(double k, int L, double m_rho, double g_rho_pipi, double Mpi) {

  
  double omega= 2*sqrt( pow(k,2)+pow(Mpi,2));

  //cout<<endl<<"d_11_der: "<<d_11_der(k, m_rho, g_rho_pipi,Mpi)<<endl;
  //cout<<"phi_der("<<k*L/(2.0*M_PI)<<"): "<<phi_der(k*L/(2.0*M_PI))<<endl;

  return ((2.0*pow(k,5)/(3.0*M_PI*pow(omega,2))))*pow(F_pi_GS_mod(omega, m_rho,g_rho_pipi,Mpi),2)/(k*d_11_der(k, m_rho, g_rho_pipi,Mpi) + (k*L/(2.0*M_PI))*phi_der(k*L/(2.0*M_PI)));

}



double Find_pipi_energy_lev(int L, int n, double m_rho, double g_rho_pipi, double Mpi) {

  //find energy level n of pi-pi bound state

  assert(n>=0);

  double Precision = 1e-8;
  double delta = 0.001;
 
 
  double a,b;

  double b_min= (n-0.15)*2.0*M_PI/L;
  b_min= (b_min >=0)?b_min:0.;
  double a_min= (n+0.15)*2.0*M_PI/L;

  auto F = [&](double k) -> double { return d_11(k, m_rho, g_rho_pipi, Mpi) + phi(k*L/(2.0*M_PI));};
  a=a_min;
  b=b_min;

  /*
  cout<<"n: "<<n<<endl;
  cout<<"m_rho: "<<m_rho<<endl;
  cout<<"g_rho_pipi: "<<g_rho_pipi<<endl;
  cout<<"M_pi: "<<Mpi<<endl;
  cout<<"d11(a): "<<d_11(a,m_rho,g_rho_pipi,Mpi)<<endl;
  cout<<"d11(b): "<<d_11(b,m_rho,g_rho_pipi, Mpi)<<endl;
  cout<<"phi(a): "<<phi(a*L/(2.0*M_PI))<<endl;
  cout<<"F(a): "<<F(a)<<endl;
  cout<<"F(b): "<<F(b)<<endl;
  */

  while (F(a)*F(b) > 0) { a+=0.15*2.0*M_PI/L;  b -= 0.15*2.0*M_PI/L; }
 



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

  cout<<"Entering V_pipi"<<endl;
  cout<<"########parameters#########"<<endl;
  cout<<"L: "<<L<<endl;
  cout<<"mrho: "<<m_rho<<endl;
  cout<<"Mpi: "<<Mpi<<endl;
  cout<<"g_rho: "<<g_rho_pipi<<endl;
  cout<<"#################"<<endl;
  cout<<"n        k        omega    ampl"<<endl;
  
  for(int i_lev=1;i_lev<=Nresonances;i_lev++) {
    double k_n= Find_pipi_energy_lev(L,i_lev, m_rho, g_rho_pipi, Mpi);
    
    double omega_n = 2.0*sqrt( pow(Mpi,2) + pow(k_n,2));
    double Ampl= Amplitude(k_n, L, m_rho, g_rho_pipi, Mpi);
    ret_val += Ampl*exp(-omega_n*t);
    cout<<i_lev<<"   "<<k_n<<"    "<<omega_n<<"     "<<Ampl<<endl;

  }
  cout<<"Exiting V_pipi"<<endl;

  return ret_val;

}


double V_pipi_infL(double t, double m_rho_infL, double g_rho_pipi_infL, double Mpi_infL) {



  auto Integrand = [&](double omega) -> double {
    return (1.0/(48.0*pow(M_PI,2)))*pow(omega,2)*pow(1.0- pow(2.0*Mpi_infL/omega,2), 3.0/2.0)*exp(-omega*t)*F_pi_GS_mod(omega, m_rho_infL, g_rho_pipi_infL,Mpi_infL);};


  return boost::math::quadrature::gauss_kronrod<double, 15>::integrate(Integrand, 2*Mpi_infL, numeric_limits<double>::infinity(), 5, relerr)  ;


  
}

