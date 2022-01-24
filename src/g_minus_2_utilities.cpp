#include "../include/g_minus_2_utilities.h"


using namespace std;

const double m_muon= 0.10565837;  //1.7768 ;
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

void Plot_kernel_K(int npoints) {

  double tmin=0.0;
  double tmax= 50*5.0; //time is in fermi

  Vfloat pts;
  Vfloat val;

  for(int i=0; i<=npoints;i++) {
    double pt = tmin + i*(tmax-tmin)/npoints;
    pts.push_back(pt);
    val.push_back( kernel_K( pt/0.1975, 1.0));
  }

  Print_To_File({}, {pts, val}, "../data/gm2/light/kernel.dat", "", "");

  
  return;

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

void Compute_clogSD_free(int scale_max) {

  auto C_cont = [](double t) -> double { return 1.0/(2.0*M_PI*M_PI*t*t*t);};
  double t0= 0.4; //fm
  double Q= pow(2.0/3.0,2)+ pow(1.0/3.0,2);
  double Delta = 0.15; //fm
  double Alpha=1.0/137.04;
  auto SD_th = [&t0, &Delta](double t) ->double { return 1.0- 1.0/(1.0 + exp(-2.0*(t-t0)/Delta));};
  double Nc=3; //three colors
  double r= 1.0;
  double mu2=0.0;
  double a=0.2;
  double resc=0.2;
 

  auto VkVk_b = [&](double z0, double q0, double p1, double p2, double p3) -> double {

		  double VkVk= (4*Nc*(mu2*a*a -r*(4 - cos(z0) - cos(p1) - cos(p2) - cos(p3))*(4 - cos(p1) - cos(p2) - cos(p3) - cos(q0)) + 
					(pow(sin(p1),2) + pow(sin(p2),2) + pow(sin(p3),2))/3.0 + sin(z0)*sin(q0)))/
		     ((mu2*a*a + pow(4 - cos(z0) - cos(p1) - cos(p2) - cos(p3),2) + pow(sin(z0),2) + pow(sin(p1),2) + pow(sin(p2),2) + pow(sin(p3),2))*
		      (mu2*a*a + pow(4 - cos(p1) - cos(p2) - cos(p3) - cos(q0),2) + pow(sin(p1),2) + pow(sin(p2),2) + pow(sin(p3),2) + pow(sin(q0),2)));

		  /* double VkVk_p0 = (4*Nc*(mu2*a*a -r*(4 - cos(q0) - cos(p1) - cos(p2) - cos(p3))*(4 - cos(p1) - cos(p2) - cos(p3) - cos(q0)) + 
					(pow(sin(p1),2) + pow(sin(p2),2) + pow(sin(p3),2))/3.0 + sin(q0)*sin(q0)))/
		     ((mu2*a*a + pow(4 - cos(q0) - cos(p1) - cos(p2) - cos(p3),2) + pow(sin(q0),2) + pow(sin(p1),2) + pow(sin(p2),2) + pow(sin(p3),2))*
		      (mu2*a*a + pow(4 - cos(p1) - cos(p2) - cos(p3) - cos(q0),2) + pow(sin(p1),2) + pow(sin(p2),2) + pow(sin(p3),2) + pow(sin(q0),2)));


		      if(VkVk + VkVk_p0< 0) crash("VkVk > VkVk_p0"); */
		  return VkVk;
		};
  
  auto VkVk_FT = [&](vector<double> const &par) -> double {

		   double pt = par[0];
		   double p0 = par[1];
		   double p1 = par[2];
		   double p2 = par[3];
		   double p3 = par[4];
		   double fact= pow(1.0/(2.0*M_PI),5);
		   double VkVk=  2.0*(VkVk_b(p0+pt,p0,p1,p2,p3) + VkVk_b(p0-pt,p0,p1,p2,p3));
		   VkVk += 2.0*( VkVk_b(M_PI-p0+pt,M_PI-p0,p1,p2,p3) + VkVk_b(M_PI-p0-pt, M_PI-p0, p1,p2,p3));
		   double VkVk_ref = 2.0*( VkVk_b(p0 +(M_PI-pt), p0,p1,p2,p3) + VkVk_b(p0- (M_PI-pt), p0,p1,p2,p3));
		   VkVk_ref += 2.0*( VkVk_b(M_PI-p0+(M_PI-pt),M_PI-p0,p1,p2,p3) + VkVk_b(M_PI-p0-(M_PI-pt), M_PI-p0, p1,p2,p3));
		   double res_prod= pow(2.0,3)*fact*(1.0/a/a/a)*VkVk;
		   double res_prod_ref = pow(2.0,3)*fact*(1.0/a/a/a)*VkVk_ref;
		   double SD_fact_1 = 0.0;
		   double SD_fact_2 =0.0;
		   for(int it=1; it<(int)2.0*t0/a;it++) {
		     double tt=it*a;
		     SD_fact_1+= cos(tt*pt/a)*SD_th(tt)*pow(tt,4)*pow(m_muon,2);
		     SD_fact_2+= cos(tt*(M_PI-pt)/a)*SD_th(tt)*pow(tt,4)*pow(m_muon,2);
		     } 
		   return a*4.0*pow(Alpha,2)*(res_prod*SD_fact_1 + res_prod_ref*SD_fact_2);
	      };

  auto VkVk_FT_summed= [&](vector<double> const & par) -> double {
			  double pt = par[0];
			  double p0 = par[1];
			  double p1 = par[2];
			  double p2 = par[3];
			  double p3 = par[4];

			 
			  Vfloat par2({ M_PI/2 -pt, p0,p1,p2,p3});
			  Vfloat par3({ pt,M_PI/2- p0,p1,p2,p3});
			  Vfloat par4({ M_PI/2 -pt, M_PI/2 -p0,p1,p2,p3});
	        


			  double res=  VkVk_FT(par) + VkVk_FT(par2) + VkVk_FT(par3)+ VkVk_FT(par4);


			  return res;


		       };

  auto VkVk_FT_summed_p0 = [&](vector<double> const & par) -> double {


			     double Vk_p0 = VkVk_b(par[1],par[1], par[2], par[3], par[4]);
			     if(Vk_p0<0) crash("Vk_pt0: "+to_string_with_precision(Vk_p0,10));
			     return VkVk_FT_summed(par);

			   };

  auto Int_SD_cont = [&](double t) ->double {

		       return 4.0*pow(Alpha,2)*C_cont(t)*SD_th(t)*pow(t,4)*pow(m_muon,2);

		     };

  //compute SD window in the continuum
  double error_amu_cont;
  double eps_cont = 1.0e-10;
  double a_mu_SD_cont= boost::math::quadrature::gauss_kronrod<double,15>::integrate(Int_SD_cont, 0.0, numeric_limits<double>::infinity() , 5, eps_cont, &error_amu_cont);


 


   //VEGAS GLS INTEGRATION

  auto display_results = [](char *title, double result, double error) ->void
			 {
			   printf ("%s ==================\n", title);
			   printf ("result = % .16f\n", result);
			   printf ("relative error  =  %.16f  %\n ", 100*error/result);
			   return;
			 };

  
    size_t calls = 5000000; //old was 1000000;
    int warmup_calls= 50000; //old was 100000;

    double bound_l[5];
    bound_l[0] = 0.0;
    bound_l[1] = 0.0;
    bound_l[2] = 0.0;
    bound_l[3] = 0.0;
    bound_l[4] = 0.0;
    double bound_u[5];
    bound_u[0] = M_PI/4;
    bound_u[1] = M_PI/4;
    bound_u[2] = M_PI;
    bound_u[3] = M_PI;
    bound_u[4] = M_PI;


    Vfloat res_r0, res_r1;
    Vfloat err_r0, err_r1;

    double a0=0.2; //starting value of the lattice spacing (does not make any sense absolute scale with massless fermions)

    for(int ir=0;ir<2;ir++) {

      if(ir==0) r=-1;
      else r= 1;
      
     
      for(int scale=0;scale<scale_max;scale++) {
	a= a0/(1.0+resc*scale);
	double a_mu_SD_lat=0.0;
	double a_mu_SD_lat_err=0.0;
	string r_tag= (ir==0)?"same_r":"opposite_r";
	cout<<"Computing SD windows for a: "<<a<<" "<<r_tag<<endl;

	double res_vegas, err_vegas;

	
	const gsl_rng_type *T;
	gsl_rng *rr;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	rr = gsl_rng_alloc(T);

   
	//to enter a generic seed mySeed
	gsl_rng_set(rr, 534543);

	gsl_monte_function_pp<decltype(VkVk_FT_summed_p0)> Fp(VkVk_FT_summed_p0);
    
	gsl_monte_function *G = static_cast<gsl_monte_function*>(&Fp);



	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(5);

 

	gsl_monte_vegas_integrate (G, bound_l, bound_u, 5, warmup_calls, rr, s,
                               &res_vegas, &err_vegas);



	printf ("converging...\n");
	int k=1;
	do
	  {
	    int new_calls=calls*k;
	    gsl_monte_vegas_integrate (G,bound_l, bound_u, 5, new_calls/5, rr, s,
                                   &res_vegas, &err_vegas);
	   
	    cout<<"k: "<<k<<" err_vegas: "<<100*err_vegas/res_vegas<<endl;
	    k*=3;
       
      }
	while (fabs (100*err_vegas/res_vegas) > 0.1);

	//display_results("vegas final", res_vegas, err_vegas);

	gsl_monte_vegas_free(s);

	gsl_rng_free(rr);

	a_mu_SD_lat = res_vegas;
        a_mu_SD_lat_err= Q*a_mu_SD_lat/(a*a*log(1/a));
      
	if(ir==0)  {res_r0.push_back( Q*(a_mu_SD_lat - a_mu_SD_cont)/(a*a*log(1/a))); err_r0.push_back( a_mu_SD_lat_err);}
	else {res_r1.push_back( Q*(a_mu_SD_lat - a_mu_SD_cont)/(a*a*log(1/a))); err_r1.push_back( a_mu_SD_lat_err);}
	cout<<"SD(lat): "<<Q*a_mu_SD_lat<<" SD(cont): "<<Q*a_mu_SD_cont<<endl;
	cout<<"Done!"<<endl;
      }
    }
  
  
    //print the result
    cout<<"printing result (same r)"<<endl;
    for(int i=0;i<(signed)res_r0.size();i++) {
      cout<<"a: "<<a0/(1.0+resc*i)<<" C^log_SD: "<< res_r0[i]<<" +- "<<err_r0[i]<<endl;
    }

    cout<<"printing result (opposite r)"<<endl;
    for(int i=0;i<(signed)res_r1.size();i++) {
      cout<<"a: "<<a0/(1.0+resc*i)<<" C^log_SD: "<< res_r1[i]<<" +- "<<err_r1[i]<<endl;
    }




    return;
}


void Compute_SD_window_Free() {

  VVfloat C_samer;
  VVfloat C_opper;
  VVfloat C_cons;

   double t0= 0.4; //fm
   double a0=1.0;
   double Q= pow(2.0/3.0,2)+ pow(1.0/3.0,2);
   double Delta = 0.15; //fm
   double Alpha=1.0/137.04;
   auto SD_th = [&t0, &Delta](double t) ->double { return 1.0- 1.0/(1.0 + exp(-2.0*(t-t0)/Delta));};
   auto int_SD_cont = [&](double t) -> double { return 4.0*pow(Alpha,2)*(1.0/(2.0*M_PI*M_PI*t*t*t))*SD_th(t)*kernel_K(t,1.0);};

   double error_amu_cont;
   double eps_cont = 1.0e-10;
   double a_mu_SD_cont= boost::math::quadrature::gauss_kronrod<double,15>::integrate(int_SD_cont, 0.0, numeric_limits<double>::infinity() , 5, eps_cont, &error_amu_cont);

  for(int k=1;k<=15;k++) {
    Vfloat rr= Read_From_File("../theo_RM123/VkVk_free_data/corr_r1_1_r2_1_L_48_T_96_mu_0.00000_scale_"+to_string(k), 1, 2);
    Vfloat rmr= Read_From_File("../theo_RM123/VkVk_free_data/corr_r1_1_r2_0_L_48_T_96_mu_0.00000_scale_"+to_string(k), 1, 2);
    Vfloat rcons= Read_From_File("../theo_RM123/VkVk_free_data/corr_r1_1_r2_2_L_48_T_96_mu_0.00000_scale_"+to_string(k), 1, 2);
    C_samer.push_back(rr);
    C_opper.push_back(rmr);
    C_cons.push_back(rcons);
  }

  Vfloat a2_Corr_48_96_samer;
  Vfloat a2_Corr_48_96_opper;

  a2_Corr_48_96_samer= Read_From_File("../theo_RM123/VkVk_free_data/48_SAMER", 1, 2);
  a2_Corr_48_96_opper= Read_From_File("../theo_RM123/VkVk_free_data/48_OPPOR", 1, 2);

  Vfloat cont_lim_samer = Sum_vectors<double>(a2_Corr_48_96_samer, C_opper[0]);
  Vfloat cont_lim_opper = Sum_vectors<double>(a2_Corr_48_96_opper, C_samer[0]);


  for(int t=1; t< (signed)cont_lim_samer.size();t++) cont_lim_samer[t]*=2.0*M_PI*M_PI*t*t*t;
  for(int t=1; t< (signed)cont_lim_opper.size();t++) cont_lim_opper[t]*=2.0*M_PI*M_PI*t*t*t;



  Vfloat a2_coeff_samer, a2_coeff_opper, a2_coeff_cons;
  Vfloat cont_samer, cont_opper, cont_cons;

  //fit correlator for each t
  for(int r=0;r<3;r++) {
  for(int t=1;t<10;t++) {
    

    Vfloat X,Y,Z;
    for(int i=5;i<15;i++) {X.push_back( a0/(1.0+i)); Y.push_back( (1.0/a0/a0/a0)*C_samer[i][t]); Z.push_back( (1.0/a0/a0/a0)*C_opper[i][t]);}
    //glue Y and Z
    Vfloat V=Y;
    Vfloat Xd = X;
    Vfloat Y_cons;
    for(int i=5;i<15;i++) {Y_cons.push_back( (1.0/a0/a0/a0)*C_cons[i][t]);}
    V.insert(V.end(), Z.begin(), Z.end());
    Xd.insert(Xd.end(), Xd.begin(), Xd.end());
    T_fit Fit_f((r==2)?X:Xd, (r==2)?Y_cons:V);
    if(r==2) Fit_f.add_pars(3.0); //cont. lim
    Fit_f.add_pars(2.5); //a^2 samer
    if(r != 2) Fit_f.add_pars(2.5); //a^2 opper
    //Fit_f.add_pars(1.0); //a^4 samer
    //Fit_f.add_pars(1.0); //a^4 opper
  

    Fit_f.ansatz = [&](const Vfloat &ip, double x, int imeas) -> double {

		     

		     double tph = t*a0;
		     double K= 1.0/(2.0*M_PI*M_PI*pow(tph,3));
		     if(r==2) return (1.0/a0/a0/a0)*K*ip[0] + 1.0*K*(1.0/a0/a0/a0)*(x/tph)*(x/tph)*ip[1] + 0.0*K*(1.0/a0*a0*a0)*pow((x/tph),4);
		     else return (1.0/a0/a0/a0)*K*(imeas<Y.size()?cont_lim_samer[t]:cont_lim_opper[t])  + 1.0*(imeas<Y.size()?ip[0]:ip[1])*K*(x/tph)*(x/tph) + 0.0*K*pow((x/tph),4);
		   };

    fit_t_res fit_result = Fit_f.fit();
    if(r==0) {a2_coeff_samer.push_back(fit_result.pars[0]); cont_samer.push_back(fit_result.pars[0]);}
    else if(r==1) {a2_coeff_opper.push_back(fit_result.pars[1]); cont_opper.push_back(fit_result.pars[0]);}
    else {a2_coeff_cons.push_back( fit_result.pars[1]); cont_cons.push_back(fit_result.pars[0]);}

   

  }
  }

  cout<<"Printing a^2 fit coefficients (same r)"<<endl;
  printV(a2_coeff_samer, "", 1);
  cout<<"Printing a^2 fit coefficients (opposite r)"<<endl;
  printV(a2_coeff_opper, "", 1);
  cout<<"Printing a^2 fit coefficients conserved"<<endl;
  printV(a2_coeff_cons, "", 1);
  cout<<"Printing cont fit coefficients (same r)"<<endl;
  printV(cont_samer, "", 1);
  cout<<"Printing cont fit coefficients (opposite r)"<<endl;
  printV(cont_opper, "", 1);
  cout<<"Printing cont fit coefficients conserved"<<endl;
  printV(cont_cons, "", 1);
  

  

  
  Vfloat alist, aSD_samer_list, aSD_opper_list;
  Vfloat Clog_SD_samer_list, Clog_SD_opper_list;
  //Compute_SD_window
  for(int k=0;k<15;k++) {
    double aSD_samer=0.0;
    double aSD_opper=0.0;
    double a= a0/(1.0+k);
    for(int t=1; t< 30;t++) { aSD_samer += a0*(1.0/pow(a0,3))*4.0*pow(Alpha,2)*SD_th(t*a0)*C_samer[k][t]*kernel_K(t*a0,1.0);}
    for(int t=1; t< 30;t++) { aSD_opper += a0*(1.0/pow(a0,3))*4.0*pow(Alpha,2)*SD_th(t*a0)*C_opper[k][t]*kernel_K(t*a0,1.0);}
    alist.push_back(a);
    aSD_samer_list.push_back(aSD_samer);
    aSD_opper_list.push_back(aSD_opper);
    Clog_SD_samer_list.push_back(aSD_samer/(a*a*log(1.0/a)));
    Clog_SD_opper_list.push_back(aSD_opper/(a*a*log(1.0/a)));
      

  }

  
  
  Vfloat cont_val(15,a_mu_SD_cont);
 


   //fit SD window
  Vfloat a2log_samer, a2log_opper, SD_cont_fit_samer, SD_cont_fit_opper;
  for(int r=0;r<2;r++) {
    

    Vfloat X,Y,Z;
    for(int i=3;i<15;i++) {X.push_back( a0/(1.0+i)); Y.push_back( aSD_samer_list[i]); Z.push_back( aSD_opper_list[i]);}
    T_fit Fit_f_SD(X, r==0?Y:Z);
    Fit_f_SD.add_pars(1.3); //cont. lim
    Fit_f_SD.add_pars(0.5); //a^2*log(1/a)
    Fit_f_SD.add_pars(1.0); //a^2

    Fit_f_SD.ansatz = [&](const Vfloat &ip, double x, int imeas) -> double {

			double K= cont_val[0];
			double xd = x/a0;
			return 1.0*ip[0]*K + 1.0*ip[1]*xd*xd*log(1/xd) + ip[2]*xd*xd;
		   };

    fit_t_res fit_result_SD = Fit_f_SD.fit();
    if(r==0) {a2log_samer.push_back(fit_result_SD.pars[1]); SD_cont_fit_samer.push_back(fit_result_SD.pars[0]);}
    else {a2log_opper.push_back(fit_result_SD.pars[1]); SD_cont_fit_opper.push_back(fit_result_SD.pars[0]);}
  }

  //print
  cout<<"Printing a^2*log(1/a) fit coefficients (same r)"<<endl;
  printV(a2log_samer, "", 1);
  cout<<"Printing a^2*log(1/a) fit coefficients (opposite r)"<<endl;
  printV(a2log_opper, "", 1);
  cout<<"Printing cont fit SD (same r)"<<endl;
  printV(SD_cont_fit_samer, "", 1);
  cout<<"Printing cont fit SD (opposite r)"<<endl;
  printV(SD_cont_fit_opper, "", 1);

  VVfloat C_samer_swap= C_samer;
  VVfloat C_opper_swap= C_opper;
  Transpose_VV<double>(C_samer_swap);
  Transpose_VV<double>(C_opper_swap);

  for(int t=1;t<C_samer[0].size()/2;t++) { Print_To_File( {}, {alist, C_samer_swap[t], C_opper_swap[t]}, "../theo_RM123/VkVk_free_data/Corr_"+to_string(t)+".dat", "", "");}


  Print_To_File({}, {alist, aSD_samer_list, aSD_opper_list, cont_val }, "../theo_RM123/VkVk_free_data/SD.out", "", "# a   SD(samer)   SD(opper)   Clog(samer)   Clog(opper)");

  return;

};



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
	cout.precision(15);
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

