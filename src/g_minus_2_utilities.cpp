#include "../include/g_minus_2_utilities.h"


using namespace std;

const double m_muon= 0.10565837;
const double relerr= 1e-12;









double kernel_K(double t, double MV) {

  double m_mu= m_muon;

  auto F = [&](double x) -> double {

    if (x<1e-40) return 0;
    return (4.0/pow(m_mu*MV,2))*(1.0/sqrt(4.0 + pow(x,2)))*pow(  (sqrt(4.0+pow(x,2))-x)/(sqrt(4.0+pow(x,2))+x),2)*( (cos(m_mu*MV*t*x)-1)/pow(x,2) + (1.0/2.0)*pow(t*m_mu*MV,2));
  };


  double err;
  double tol_kernel= 1e-14;
  return boost::math::quadrature::gauss_kronrod<double,61>::integrate(F, 0.0, numeric_limits<double>::infinity() , 5, tol_kernel, &err);


}


double Zeta_function_laplacian_Luscher(double z) {


  
  double inf = 1e-20;
  double tol_zeta= 1e-10;

  assert(z>= 0) ;

  //compute Zeta_function using Luscher parametrization

  double z2= pow(z,2);

  //l2 is cutoff
  int l2= z2 + 2;
  // if(degeneracy(l2) == 0) l2++;

  
  int thresh = l2+ 10;

  double t0= z2>2.0?M_PI/z2:0.5*M_PI;

  //t0=1.0;
  

  //sum first terms
  double val=0;
  double sum0=0.0;
  for(int m=0; m < l2; m++) {
    if(degeneracy(m) > 0 ) {
      sum0 += (1.0/sqrt(4*M_PI))*degeneracy(m)*(1.0/(m-z2));
    }
  }


  auto F2 = [&](double t) -> double {

	      double sum_t=0.0;
	      double exp_t= exp(-t);
	      double exp_qt= exp(-pow(M_PI,2)/t);
	      for(int m=0; m<l2;m++) sum_t += (1.0/sqrt(4.0*M_PI))*degeneracy(m)*(1.0/pow(2.0*M_PI,3))*pow(exp_t,m);
	      
	      double heat_kernel_red= (exp(t*z2)-1.0)/(pow(4.0*M_PI,2)*pow(t,1.5)); //m=0 term
	      for(int m=1; m<= thresh;m++) {
		heat_kernel_red += ( exp(t*z2)*(1.0/sqrt(4*M_PI)))*degeneracy(m)*( 1.0/pow(4.0*M_PI*t,3.0/2.0))*pow(exp_qt,m);
	      }
	      
	      double res= pow(2.0*M_PI,3)*(heat_kernel_red -exp(t*z2 + log(sum_t)));

	     
	      if( isnan(res)) crash("In integral repr. Zeta function_1 t: "+to_string_with_precision(t,3)+" m: "+to_string_with_precision(z2, 3)+" , integrand is nan");
	    
	      return res;

	      
	    };
    
  auto F = [&](double t) -> double {
	     
	     double heat_kernel_red=0;
	     double exp_t= exp(-t);
	     for(int m=l2;m<=thresh; m++) {
	       heat_kernel_red += (1.0/pow(2.0*M_PI,3))*degeneracy(m)*pow(exp_t,m);
	     }
	     
	     
	     double res= pow(2.0*M_PI,3)*(1.0/sqrt(4*M_PI))*exp(t*z2 + log(heat_kernel_red));

	     
	     if( isnan(res)) crash("In integral repr. Zeta function_2 t: "+to_string_with_precision(t,3)+" m: "+to_string_with_precision(z2, 3)+" , integrand is nan");
	    
	     return res;
	   };


  double int_val_2= boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F, t0, 30.0 , 5, tol_zeta);
  double int_val_1= boost::math::quadrature::gauss_kronrod<double, 15>::integrate(F2, inf, t0 , 5, tol_zeta);

  double int_val = int_val_2 + int_val_1 + sum0 -M_PI/sqrt(t0);

  if( isnan( int_val)) crash("integral in Generalized zeta is nan");
  return int_val;
}


double tan_phi(double z) {

  int n2= (int)pow(z,2); //check whether is integer.
  int n2_p= n2+1;
  if(  (fabs( pow(z,2) - n2) < 1e-18 && degeneracy(n2*n2) > 0) ||  (fabs(pow(z,2) - n2_p) < 1e-18 && degeneracy(n2_p*n2_p) > 0) ) return 0.0;

  return -1.0*pow(M_PI, 1.5)*z/Zeta_function_laplacian_Luscher(z);

}

double tan_phi_der(double z) {

  auto tf = [&](double t) { return tan_phi(t);};

  //double x=z;

  double result, err;

  gsl_function_pp<decltype(tf)> Fp(tf);

  gsl_function *G = static_cast<gsl_function*>(&Fp);
  if(z < 2e-3) {
     gsl_deriv_forward(G, z, 1e-3, &result, &err);
  }
  else {
  gsl_deriv_central(G, z, 1e-3, &result, &err);
  }

  //double result= boost::math::differentiation::finite_difference_derivative(tf,x, &err);

  if(err/result > 1.0e-1) crash("numerical derivative in tan_phi'(z) has reached a low level of accuracy in phi_der. z: "+to_string_with_precision(z,10)+" rel_err: "+to_string_with_precision(err/result,10));


  if(isnan(result)) crash("Function tan_phi_der(z) returned nan");

  return result;

}

double phi_der(double z) {

  auto tf = [&](double t) { return tan_phi(t);};

  //double x=z;

  double result, err;

  gsl_function_pp<decltype(tf)> Fp(tf);

  gsl_function *G = static_cast<gsl_function*>(&Fp);
  if(z < 2e-3) {
     gsl_deriv_forward(G, z, 1e-3, &result, &err);
  }
  else {
  gsl_deriv_central(G, z, 1e-3, &result, &err);
  }

  //double result= boost::math::differentiation::finite_difference_derivative(tf,x, &err);

  if(err/result > 1.0e-1) crash("numerical derivative in phi'(z) has reached a low level of accuracy in phi_der. z: "+to_string_with_precision(z,10)+" rel_err: "+to_string_with_precision(err/result,10));

  double res= result/(1.0+ pow(tan_phi(z),2));

  if(isnan(res)) crash("Function phi_der(z) returned nan");

  return res;

}

double phi_der_for_back(double z, int mode) {

  auto tf = [&](double t) { return phi(t);};

  //double x=z;

  double result, err;

  gsl_function_pp<decltype(tf)> Fp(tf);

  gsl_function *G = static_cast<gsl_function*>(&Fp);
  if(mode==1)  gsl_deriv_forward(G, z, 1e-3, &result, &err);
  else if(mode==-1)  gsl_deriv_backward(G, z, 1e-3, &result, &err);
  else crash("mode: "+to_string( mode)+" not yet implemented");

  //double result= boost::math::differentiation::finite_difference_derivative(tf,x, &err);

  if(err/result > 1.0e-2) crash("numerical for_back derivative in phi'(z) has reached a low level of accuracy in phi_der. z: "+to_string_with_precision(z,10)+" rel_err: "+to_string_with_precision(err/result,10));


  if(isnan(result)) crash("Function phi_der_for_back(z) returned nan");

  return result;

}

void Zeta_function_zeroes(int Nzeros, Vfloat &res) { //using Newton method

  res.clear();

  auto F = [&](double z) -> double {

	     return Zeta_function_laplacian_Luscher(sqrt(z));

	   };

  auto dF = [&](double z) -> double {

	      auto tf = [&](double t) {return Zeta_function_laplacian_Luscher(sqrt(t));};
	      double der, err_der;
	      gsl_function_pp<decltype(tf)> Fp2(tf);

	      gsl_function *G = static_cast<gsl_function*>(&Fp2);
	      gsl_deriv_forward(G, z, 1e-5, &der, &err_der);
	      return der;
	     

	   };

  auto fdF = [&](double z, double* f, double *df) -> void {
	       *f = F(z);
	       *df = dF(z);
	       return;
	     };


  double Precision= 1e-10;
  double offset = 1e-10;

  int N2old=0;
  
  for(int izero=0; izero<Nzeros;izero++) {

    int status;   
    int iter = 0, max_iter = 100;
    const gsl_root_fdfsolver_type *TT;
    gsl_root_fdfsolver *ss;
    double x0,x;

    if(izero==0) x= offset;
    else {
      int temp=N2old+1;
      while(degeneracy(temp) == 0) temp++;
      x= temp + offset;
      N2old=temp;
    }

    
    gsl_function_fdf_pp<decltype(F), decltype(dF), decltype(fdF)> Fp(F, dF, fdF);
    
    gsl_function_fdf *FDF = static_cast<gsl_function_fdf*>(&Fp);

    

    
    TT = gsl_root_fdfsolver_newton;
    ss = gsl_root_fdfsolver_alloc (TT);
    gsl_root_fdfsolver_set(ss, FDF, x);

    
    do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate(ss);
      x0 = x;
      x = gsl_root_fdfsolver_root (ss);
      status = gsl_root_test_delta (x, x0, 0, Precision);

           
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    if(status != GSL_SUCCESS) {

      cout<<"gsl_root_fdf_solver_newton was unable to find zero of the Luscher zeta function for  izero: "<<izero<<endl;
      crash("Aborting...");
    }

    gsl_root_fdfsolver_free(ss);

    res.push_back(x);
  }

  

  return;
}

double phi(double z) {

  int n2= (int)pow(z,2); //check whether is integer.
  int n2_p= n2+1;
  if( (fabs( pow(z,2) - n2) < 1e-16 && degeneracy(n2*n2) > 0) || ( fabs(pow(z,2) - n2_p) < 1e-16 && degeneracy(n2_p*n2_p) > 0) ) return 0.0;
  double sum= Zeta_function_laplacian_Luscher(z);

  double res= atan( -1.0*pow(M_PI,1.5)*z/sum);

  if(isnan(res)) crash("Function phi(z) returned nan");
  
  return res;

}

double LL_functions::Vdual(double t, double m_rho, double Edual, double Rdual) {


  return (5.0/(18.0*pow(M_PI,2)))*(Rdual/pow(t,3))*exp(-(m_rho+Edual)*t)*( 1.0 + (m_rho+Edual)*t + 0.5*pow( (m_rho+Edual)*t,2));

}


double LL_functions::Gamma_rpp(double omega, double g_rho_pipi, double Mpi) {

  double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));
  return pow(g_rho_pipi,2)*pow(k,3)/(6.0*M_PI*pow(omega,2));  
}

double LL_functions::Gamma_rpp_der(double omega, double g_rho_pipi, double Mpi) {

  double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));

  double dk_domega = omega/(4.0*k);

  return (pow(g_rho_pipi,2)/(6.0*M_PI))*( -2.0*pow(k,3)/pow(omega,3) + 3.0*pow(k,2)*dk_domega/pow(omega,2));
  //return 2*k*( (8*pow(Mpi,2) + pow(omega,2))/(4*pow(omega,3)))*pow(g_rho_pipi,2)/(6.0*M_PI);
  
}

double LL_functions::h(double omega, double g_rho_pipi, double Mpi) {
  double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));
  return pow(g_rho_pipi,2)*pow(k,3)*(2.0/(omega*6.0*pow(M_PI,2)))*log( (omega + 2*k)/(2.0*Mpi));

}

double LL_functions::h_prime(double omega, double g_rho_pipi, double Mpi) {

  double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));
  return pow(g_rho_pipi,2)*pow(k,2)*(1.0/(6.0*pow(M_PI,2)*omega))*( 1.0 + (1.0+ 2.0*pow(Mpi/omega,2))*(omega/k)*log( (omega+2*k)/(2.0*Mpi)));

}

double LL_functions::h_second(double omega, double g_rho_pipi, double Mpi) {

   double k= sqrt( pow(omega,2)/4.0 - pow(Mpi,2));

   return (1.0/6.0)*pow(g_rho_pipi/M_PI,2)*(3.0/4.0 + 2.0*pow(Mpi/omega,2) + (1.0-2.0*pow(Mpi/omega,2))*(k/omega)*log( (omega + 2.0*k)/(2.0*Mpi)));

}



double LL_functions::A_pipi_0(double m_rho, double g_rho_pipi, double Mpi, double kappa) {

  return h(m_rho,g_rho_pipi, Mpi) - (m_rho/2.0)*h_prime(m_rho, g_rho_pipi, Mpi) + pow(g_rho_pipi,2)*pow(Mpi,2)/(6.0*pow(M_PI,2)) + kappa*(pow(m_rho,2)/8.0)*(h_second(m_rho, g_rho_pipi, Mpi) - h_prime(m_rho, g_rho_pipi, Mpi)/m_rho) ;  

}


double LL_functions::Re_A_pipi(double omega, double m_rho, double g_rho_pipi, double Mpi, double kappa) {

  return h(m_rho, g_rho_pipi, Mpi) + (pow(omega,2)- pow(m_rho,2))*(h_prime(m_rho, g_rho_pipi,Mpi)/(2.0*m_rho)) - h(omega, g_rho_pipi,Mpi) + (kappa/2.0)*pow((pow(omega,2)-pow(m_rho,2)),2)*(1.0/(pow(2.0*m_rho,2)))*( h_second(m_rho,g_rho_pipi,Mpi) - (1.0/m_rho)*h_prime(m_rho, g_rho_pipi,Mpi));

  
}

double LL_functions::Im_A_pipi(double omega, double g_rho_pipi, double Mpi) {

  return omega*Gamma_rpp(omega,g_rho_pipi,Mpi);

}

Pfloat LL_functions::A_pipi(double omega, double m_rho, double g_rho_pipi, double Mpi, double kappa) {

  return make_pair(Re_A_pipi(omega, m_rho, g_rho_pipi,Mpi, kappa), Im_A_pipi(omega, g_rho_pipi,Mpi));

}

double LL_functions::F_pi_GS_mod(double omega, double m_rho, double g_rho_pipi, double Mpi, double kappa) {

  double Re = pow(m_rho,2) -pow(omega,2) - Re_A_pipi(omega,m_rho, g_rho_pipi,Mpi, kappa);
  double Im = -1.0*Im_A_pipi(omega, g_rho_pipi, Mpi);
  return fabs((pow(m_rho,2) - A_pipi_0(m_rho, g_rho_pipi,Mpi, kappa)))/( sqrt( pow(Re,2)+ pow(Im,2))); 

}

double LL_functions::cot_d_11(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa) {


  double omega= 2.0*sqrt( pow(Mpi,2) + pow(k,2));

  double Im_a_pipi= Im_A_pipi(omega, g_rho_pipi, Mpi);

  double res = (pow(m_rho,2) - pow(omega,2) - Re_A_pipi(omega, m_rho, g_rho_pipi,Mpi, kappa))/Im_a_pipi;

  return res;

}

double LL_functions::d_11(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa) {


  double d_11=  cot_d_11(k, m_rho, g_rho_pipi, Mpi, kappa);
  
  //double sgn= 1.0;

  //if(d_11 < 0.0) sgn = -1.0;

  double ret_val = atan(1.0/d_11);
  
  //ret_val=  asin( sgn/sqrt(1.0 + pow(d_11,2)));

  return ret_val;

}


double LL_functions::cot_d_11_der(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa) {


  double omega = 2.0*sqrt( pow(Mpi,2) + pow(k,2));

  double d_omega_d_k = 4.0*k/omega;

  double num= ( pow(m_rho,2) -pow(omega,2) - Re_A_pipi(omega, m_rho, g_rho_pipi, Mpi, kappa));


   double den = Im_A_pipi(omega, g_rho_pipi, Mpi);

   double num_der = (-2.0*omega -2.0*omega*(h_prime(m_rho, g_rho_pipi,Mpi)/(2.0*m_rho)) + h_prime(omega,g_rho_pipi,Mpi)    -(2.0*kappa*omega*( pow(omega,2) -pow(m_rho,2)   ))*(1.0/(pow(2.0*m_rho,2)))*( h_second(m_rho,g_rho_pipi,Mpi) - (1.0/m_rho)*h_prime(m_rho, g_rho_pipi,Mpi)))*d_omega_d_k;

  
  double den_der = (Gamma_rpp(omega, g_rho_pipi, Mpi) + omega*Gamma_rpp_der(omega, g_rho_pipi, Mpi))*d_omega_d_k;

 
  double res_2 = cot_d_11(k, m_rho, g_rho_pipi,Mpi, kappa)*( num_der/num - den_der/den);

 

  return res_2;


}


double LL_functions::d_11_der(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa) {



  return -1.0*(1.0/(1.0 + pow(cot_d_11(k,m_rho,g_rho_pipi,Mpi, kappa),2)))*cot_d_11_der(k, m_rho, g_rho_pipi, Mpi, kappa);

}

double LL_functions::d_11_der_num(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa) {

  auto F = [&](double x) { return d_11(x, m_rho, g_rho_pipi,Mpi, kappa);};

  //double x=z;

  double result, err;

  gsl_function_pp<decltype(F)> Fp(F);

  gsl_function *G = static_cast<gsl_function*>(&Fp);

  gsl_deriv_central(G, k, 1e-3, &result, &err);


  if(err/result > 1.0e-5) crash("numerical derivative in d11'(z) has reached a low level of accuracy in d_11_der_num. rel_err: "+to_string_with_precision(err/result,4));

  return result;


}





double LL_functions::Amplitude(double k, double L, double m_rho, double g_rho_pipi, double Mpi, double kappa) {

  
  double omega= 2.0*sqrt( pow(k,2)+pow(Mpi,2));


  return ((2.0*pow(k,5)/(3.0*M_PI*pow(omega,2))))*pow(F_pi_GS_mod(omega, m_rho,g_rho_pipi,Mpi, kappa),2)/(k*d_11_der(k, m_rho, g_rho_pipi,Mpi, kappa) + (k*L/(2.0*M_PI))*phi_der_spline(k*L/(2.0*M_PI)));

}



void LL_functions::Find_pipi_energy_lev(double L, double m_rho, double g_rho_pipi, double Mpi, double kappa,  Vfloat &res) {

  //find energy level n of pi-pi bound state

  Vfloat divergent_levels= Luscher_zeroes;

  double z_crit = (L/(2.0*M_PI))*sqrt( pow(m_rho,2)/4 - pow(Mpi,2));
  double z2_crit= z_crit*z_crit;


  //add z2_crit to divergent_levels
  divergent_levels.push_back(z2_crit);
  sort(divergent_levels.begin(), divergent_levels.end());


  res.clear();

  auto F = [&](double k) -> double {
	     return d_11(k, m_rho, g_rho_pipi, Mpi, kappa) + phi_spline(k*L/(2.0*M_PI));
	   };

  auto F_exact = [&](double k) -> double {
		   return d_11(k, m_rho, g_rho_pipi, Mpi,kappa) + phi(k*L/(2.0*M_PI));
	   };

  auto dF = [&](double k) -> double {

	      return -cot_d_11_der(k, m_rho, g_rho_pipi, Mpi,kappa)/pow(cot_d_11(k, m_rho, g_rho_pipi, Mpi,kappa),2) + (L/(2.0*M_PI))*tan_phi_der(k*L/(2.0*M_PI));
	   };

  auto fdF = [&](double k, double* f, double *df) -> void {
	       *f = F(k);
	       *df = dF(k);
	       return;
	     };


  double Precision = 1e-4;

  
  for(int ilev=1; ilev<=Nres;ilev++) {
  
    int status;
    int iter = 0, max_iter = 100;
    // const gsl_root_fdfsolver_type *T;
    const gsl_root_fsolver_type *T_Br;
    //gsl_root_fdfsolver *s;
    gsl_root_fsolver *s_Br;
    double x0,x;
  
    
    // gsl_function_fdf_pp<decltype(F), decltype(dF), decltype(fdF)> Fp(F, dF, fdF);
    gsl_function_pp<decltype(F)> Fp_Brent(F);
    gsl_function_pp<decltype(F_exact)> Fp_Brent_exact(F_exact);
    
    // gsl_function_fdf *FDF = static_cast<gsl_function_fdf*>(&Fp);
    gsl_function *F_Brent = static_cast<gsl_function*>(&Fp_Brent);
    gsl_function *F_Brent_exact = static_cast<gsl_function*>(&Fp_Brent_exact);
    
 
    
    //   T = gsl_root_fdfsolver_newton;
    // s = gsl_root_fdfsolver_alloc (T);
    T_Br= gsl_root_fsolver_brent;
    s_Br= gsl_root_fsolver_alloc(T_Br);
   

    double XL= (2.0*M_PI/L)*sqrt(divergent_levels[ilev-1]);
    double XH= (2.0*M_PI/L)*sqrt(divergent_levels[ilev]);
    double XL_start= XL;
    double XH_start = XH;
    double DX= (XH -XL)*1e-3;
    XL+=(XH-XL)*1e-6;
    XH-=(XH-XL)*1e-6;

  
    gsl_set_error_handler_off();
    int Brent_status=gsl_root_fsolver_set(s_Br,F_Brent, XL, XH);

    //gsl_root_fdfsolver_set(s, FDF, x);
    
    while(Brent_status != GSL_SUCCESS) {
      XL+=DX;
      XH-=DX;
      Brent_status= gsl_root_fsolver_set(s_Br,F_Brent, XL, XH);
     

      
      if(XL > XH){
	cout.precision(16);
	cout<<"Brent algo is not able to bracket the root"<<endl;
	cout<<"status: "<<Brent_status<<endl;
	cout<<"PRINTING INFO: "<<endl;
	cout<<"Mrho: "<<m_rho<<endl;
	cout<<"Mpi: "<<Mpi<<endl;
	cout<<"g: "<<g_rho_pipi<<endl;
	cout<<"kappa: "<<kappa<<endl;
	cout<<"L: "<<L<<endl;
	cout<<"ilev: "<<ilev<<endl;
	cout<<"Bracketing values: "<<XL<<" "<<XH<<endl;
	cout<<"Starting bracketing values: "<<XL_start<<" "<<XH_start<<endl;
	cout<<"Values of F at bracketing values (spline) : "<<F(XL)<<" "<<F(XH)<<endl;
	cout<<"Values of F at bracketing values (exact) : "<<F_exact(XL)<<" "<<F_exact(XH)<<endl;
	cout<<"Values of F at starting bracketing values (spline) : "<<F(XL_start)<<" "<<F(XH_start)<<endl;
	cout<<"Values of F at starting bracketing values (exact) : "<<F_exact(XL_start)<<" "<<F_exact(XH_start)<<endl;
	printV(divergent_levels, "Printing diverging levels....", 1);
	cout<<"rho resonance divergence is at z^2: "<<z2_crit<<endl;
	crash("Exiting....");
      }
    }

    double x_lo,x_hi;
    
    do
    {
      iter++;
      status = gsl_root_fsolver_iterate(s_Br);
      //x0 = x;
      x = gsl_root_fsolver_root (s_Br);
      //status = gsl_root_test_delta (x, x0, 0, Precision);

      x_lo = gsl_root_fsolver_x_lower(s_Br);
      x_hi = gsl_root_fsolver_x_upper(s_Br);
      status = gsl_root_test_interval (x_lo, x_hi, 0.0, Precision);
      
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    if(status != GSL_SUCCESS) {

      cout<<"gsl_root_solver was unable to find a solution to the Luscher quantization equation for ilev: "<<ilev<<endl;
      cout<<"Printing info: "<<endl;
      cout<<"L: "<<L<<endl;
      cout<<"Mrho: "<<m_rho<<endl;
      cout<<"g: "<<g_rho_pipi<<endl;
      cout<<"Mpi: "<<Mpi<<endl;
      cout<<"kappa: "<<kappa<<endl;
      cout<<"Last root position: "<<x<<endl;
      cout<<"divergent level position 1: "<<(2.0*M_PI/L)*sqrt(divergent_levels[ilev-1]+2e-2)<<endl;
      cout<<"divergent level position 2: "<<(2.0*M_PI/L)*sqrt(divergent_levels[ilev]+2e-2)<<endl;
      crash("Aborting...");
    }


    //if(switch_to_exact_function) {
    //cout.precision(12);
    // cout<<"ilev: "<<ilev<<" switched to exact function evaluation"<<endl;
    // cout<<"root(approximate): "<<z_failed<<" root exact: "<<x<<endl;
    //}

    //if(!switch_to_exact_function) { root_found= F_exact(x_lo)*F_exact(x_hi) <= 0; z_failed=x;}
    //if(!root_found) switch_to_exact_function=true;
    
    // gsl_root_fdfsolver_free (s);
    gsl_root_fsolver_free(s_Br);

    res.push_back(x);

  }

  

  return;
  
}


double LL_functions::V_pipi(double t, double L, double m_rho, double g_rho_pipi, double Mpi, double kappa,  Vfloat &Knpp) {

  double ret_val=0;

 
  /*
  cout<<"Entering V_pipi"<<endl;
  cout<<"########parameters#########"<<endl;
  cout<<"L: "<<L<<endl;
  cout<<"mrho: "<<m_rho<<endl;
  cout<<"Mpi: "<<Mpi<<endl;
  cout<<"g_rho: "<<g_rho_pipi<<endl;
  cout<<"#################"<<endl;
  cout<<"n        k        omega"<<endl;
  for(int i_lev=0;i_lev<Nres;i_lev++)  cout<<i_lev+1<<"   "<<Knpp[i_lev]<<"    "<<2.0*sqrt(pow(Mpi,2)+pow(Knpp[i_lev],2))<<endl;

  */
  
  for(int i_lev=0;i_lev<Nres;i_lev++) {
    double k_n= Knpp[i_lev];  
    double omega_n = 2.0*sqrt( pow(Mpi,2) + pow(k_n,2));
    double Ampl= Amplitude(k_n, L, m_rho, g_rho_pipi, Mpi,kappa);
    ret_val += Ampl*exp(-omega_n*t);
   

  }
  // cout<<"Exiting V_pipi"<<endl;

  return ret_val;

}


double LL_functions::V_pipi_infL(double t, double m_rho_infL, double g_rho_pipi_infL, double Mpi_infL, double kappa_infL) {

  double tol_infL= 1e-14;


  auto Integrand = [&](double omega) -> double {
		     return (1.0/(48.0*pow(M_PI,2)))*pow(omega,2)*pow(1.0- pow(2.0*Mpi_infL/omega,2), 3.0/2.0)*exp(-omega*t)*pow(F_pi_GS_mod(omega, m_rho_infL, g_rho_pipi_infL,Mpi_infL,kappa_infL),2);};


  return boost::math::quadrature::gauss_kronrod<double, 61>::integrate(Integrand, 2*Mpi_infL, numeric_limits<double>::infinity(), 5, tol_infL)  ;


  
}

