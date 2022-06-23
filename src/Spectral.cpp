#include "../include/Spectral.h"


const double step_size = 0.1; //in units of sigma
const bool INCLUDE_ERRORS= true;
double lambda= INCLUDE_ERRORS?0.9:0.0;
bool FIND_OPTIMAL_LAMBDA= true;
string COV_MATRIX_MODE = "";
const int Nmoms=1;
const int alpha=2;
const bool Use_balance_condition = true;


using namespace std;


double Get_exact_gauss(const double &E, const double &m , const double &s, const double &E0) {

 double e= exp(-0.5*(E-m)*(E-m)/(s*s));
 double norm= s*( 2.0 + 0.0*erf( (m-E0)/(s*sqrt(2))))*sqrt(M_PI/2) ;

 return e/norm;

}

PrecFloat Get_exact_gauss(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0) {
  PrecFloat e = exp( -0.5*(E-m)*(E-m)/(s*s));
  PrecFloat norm= s*( 2 + 0.0*erf( (m-E0)/(s*sqrt(PrecFloat(2)))))*sqrt(precPi()/PrecFloat(2)) ;
  return e/norm;
}


PrecFloat Get_exact_gaussE2(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0) {
  return Get_exact_gauss(E, m, s, E0)*E*E;
}

PrecFloat Get_gaussE2_norm(const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0) {

  PrecFloat a1 = exp( -sqr(m-E0)/(s*s*2))*(E0+m)*sqr(s);
  PrecFloat norm= s*( 1.0 + erf( (m-E0)/(s*sqrt(PrecFloat(2)))))*sqrt(precPi()/PrecFloat(2)) ;

  return a1/norm + sqr(m)+sqr(s);
  
}

PrecFloat Get_legoE2_norm(const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0) {

  return m*m+ s*s/12;  

}

PrecFloat Get_exact_lego(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0) {
  PrecFloat a1, a2;

  a1= m - s/2;
  a2= m + s/2;
  if( E >= a1 && E <= a2) return 1/s;
  else return PrecFloat(0);
}

PrecFloat Get_exact_legoE2(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0) {
  PrecFloat a1, a2;

  return E*E*Get_exact_lego(E,m,s,E0);
}

PrecFloat Get_exact_cauchy(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0) {
  PrecFloat Z = 0.5 + (1/precPi())*atan( (m-E0)/s);
  PrecFloat c = (1/precPi())*s/( sqr(E-m) + s*s);
  return c/Z;
}

PrecFloat Get_exact_func(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f ) {

  PrecFloat res;

   if(SMEARING_FUNC == "GAUSSIAN") { res = Get_exact_gauss(E, m, s, E0); }
   else if(SMEARING_FUNC == "GAUSSIAN_E2") {     res = Get_exact_gaussE2(E, m, s, E0);   }
   else if(SMEARING_FUNC == "LEGO") {  res = Get_exact_lego(E, m, s, E0);  } 
   else if(SMEARING_FUNC == "LEGO_E2") {  res = Get_exact_legoE2(E, m, s, E0); }  
   else  res= f(E,m,s,E0);

   return res ;

}
   


PrecFloat BaseFunc(const PrecFloat& E, int t, int T) { return exp(-E*t) + exp( -E*T+ E*t);}

PrecFloat F_gauss(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t) {

  PrecFloat x,y,z;

  y = exp( -( 2*m - s*s*t)*t/2);
  x = erfc( (E0-m)/(s*sqrt(PrecFloat(2))));
  z = erfc(  (E0 -m + s*s*t)/(s*sqrt(PrecFloat(2))));


  return y*z/x;
  
}


PrecFloat F_gaussE2(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t) {

  PrecFloat x,y,z;

  PrecFloat y1, z1;

  PrecFloat y2, z2; 

  y = exp( -( 2*m - s*s*t)*t/2);
  x = erfc( (E0-m)/(s*sqrt(PrecFloat(2))));
  z = erfc(  (E0 -m + s*s*t)/(s*sqrt(PrecFloat(2))));

  y1 =  exp( ( -2*m + s*s*t)*t/2)*(-m + s*s*t);
  y2 =  exp( ( -2*m + s*s*t)*t/2)*( sqr(s*s*t -m) + sqr(s));

  z1 = -exp( -sqr(E0 -m + s*s*t)/(s*s*2))*s*sqrt( PrecFloat(2)/precPi());
  z2 = exp( -sqr(E0 -m + s*s*t)/(s*s*2))*sqrt( PrecFloat(2)/precPi())*s*(E0 - m + s*s*t);


  return (y2*z + y*z2 + +2*y1*z1)/x;
  
}


PrecFloat F_lego(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t) {

  return ( exp( -max(E0, m - s/2)*t) - exp( -(m + s/2)*t))/(s*t);

}


PrecFloat F_legoE2(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t) {

  PrecFloat num = exp( -max(E0, m - s/2)*t) - exp( -(m + s/2)*t);
  PrecFloat den = s*t;
  PrecFloat num1 = -max(E0, m - s/2)*exp( -max(E0, m - s/2)*t) + (m + s/2)*exp( -(m + s/2)*t);
  PrecFloat num2 = sqr(max(E0, m - s/2))*exp( -max(E0, m - s/2)*t) - sqr((m + s/2))*exp( -(m + s/2)*t);
  PrecFloat den1 = s;

  return num2/den - 2*num1*den1/sqr(den) + 2*num*sqr(den1)/(den*den*den);


}



//PrecFloat F_cauchy(const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int t) {

//  return ;

//}

PrecFloat aE0(const PrecFloat &E0, int t, int n) { if(n==0) { return exp(-E0*t)/t;} else return pow(E0, PrecFloat(n))*exp(-E0*t)/t + PrecFloat(n)*aE0(E0,t,n-1)/t;}

void Get_Atr(PrecMatr& Atr, const PrecFloat &E0, int T, int tmin, int tmax)  {

  Atr.resize(tmax-tmin+1, tmax-tmin+1);

  for(int t=tmin;t<= tmax; t++)
    for(int r=tmin; r<= tmax; r++) Atr(t-tmin,r-tmin) = aE0(E0, t+r,alpha) + aE0(E0, T - t +r, alpha) + aE0(E0, T+t -r, alpha) + aE0(E0, 2*T -t -r, alpha);
    
  
  return;
}


void Get_Rt(PrecVect& Rt, const PrecFloat &E0,  int T, int tmin, int tmax) {

  Rt.resize(tmax-tmin+1);

  

  for(int t=tmin;t<=tmax;t++) {
    PrecFloat x= exp(-E0*t)/t;
    PrecFloat y = exp(-E0*(T-t))/(T-t);
    Rt(t-tmin) =  x+y;
  }

  return;

}

void Get_ft(PrecVect& ft, const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int T, int tmin, int tmax, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f) {

  ft.resize(tmax-tmin+1);

  if(SMEARING_FUNC == "GAUSSIAN") {
       for(int t=tmin;t<=tmax;t++) ft(t-tmin) = F_gauss(E0, m, s, t) + F_gauss(E0, m,s, T-t);
  }
  else if(SMEARING_FUNC == "GAUSSIAN_E2") {
       for(int t=tmin;t<=tmax;t++) ft(t-tmin) = F_gaussE2(E0, m, s, t) + F_gaussE2(E0, m,s, T-t);
  }
  else if(SMEARING_FUNC == "LEGO") {
       for(int t=tmin;t<=tmax;t++) ft(t-tmin) = F_lego(E0, m, s, t) + F_lego(E0, m,s, T-t);
  }
  else if(SMEARING_FUNC == "LEGO_E2") {
       for(int t=tmin;t<=tmax;t++) ft(t-tmin) = F_legoE2(E0, m, s, t) + F_legoE2(E0, m,s, T-t);
  }

  else {
    for(int t=tmin;t<=tmax;t++) {

   
    const auto ftT=
      [&f, &m,&s,&E0, &t, &T](const PrecFloat& x) -> PrecFloat
    {
    
      return pow(x,PrecFloat(alpha))*f(x,m,s,E0)*(exp(-x*t) + exp(-x*(T-t))) ;
    };


    
    ft(t-tmin) =   integrateUpToInfinite(ftT, E0.get());
    }
  }

  return;
}

void Get_bt(PrecVect& bt,const PrecFloat &E,  int T, int tmin, int tmax) {

  bt.resize(tmax-tmin+1);

  for(int t=tmin;t<=tmax;t++) bt(t-tmin) = BaseFunc(E, t, T);

  return;


}

PrecFloat Get_norm_constraint(PrecFloat &m, PrecFloat &s, PrecFloat &E0, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f) {

  PrecFloat res=1;

  if(SMEARING_FUNC == "GAUSSIAN_E2") {
     res= Get_gaussE2_norm(m, s, E0);
   }
  else if(SMEARING_FUNC == "LEGO_E2") {
    return Get_legoE2_norm(m,s,E0);
  }
  else {

    const auto f1=
      [&f, &m,&s,&E0](const PrecFloat& x) -> PrecFloat
    {
      return pow(x,PrecFloat(alpha))*f(x, m, s, E0) ;
    };

  return  integrateUpToInfinite(f1, E0.get());
  }
  return res;
}


PrecFloat Get_M2(PrecFloat &m, PrecFloat &s, PrecFloat &E0,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f) {


  const auto f2=
    [&f, &m,&s,&E0](const PrecFloat& x) -> PrecFloat
    {
      return pow(x,PrecFloat(alpha))*f(x, m, s, E0)*f(x,m,s,E0) ;
    };

  return  integrateUpToInfinite(f2, E0.get());

}

void Get_M_N(PrecFloat &m, PrecFloat &s, PrecFloat &E0,  const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f,  PrecVect &M_n) {

  M_n.resize(Nmoms);

  for(int n=0;n<Nmoms;n++) {

     const auto fn=
       [&f, &m,&s,&E0, &n](const PrecFloat& x) -> PrecFloat
    {
      return pow(x,PrecFloat(alpha))*f(x, m, s, E0)*pow(x,PrecFloat(n)) ;
    };

     M_n(n) = integrateUpToInfinite(fn,E0.get());   

  }

  return;

}

void Get_M_tilde_N(const PrecVect &ft, const PrecMatr &Atr_inv, vector<PrecVect> &Rt_n,  PrecVect &M_tilde_n) {

  M_tilde_n.resize(Nmoms);

  for(int n=0;n<Nmoms;n++) {

    M_tilde_n(n) = Rt_n[n].transpose()*Atr_inv*ft;   

  }

  return;

}


void Get_Rt_up_to_N(PrecFloat &E0, int T, int tmin, int tmax, vector<PrecVect> &Rt_n) {

  Rt_n.resize(Nmoms);


  for(int n=0;n<Nmoms;n++) {

    Rt_n[n].resize(tmax-tmin+1);

    for(int t=tmin;t<=tmax;t++) Rt_n[n](t-tmin) = aE0(E0,t,alpha+n) + aE0(E0,T-t,alpha+n);

  }


  return;

}

void Get_G_matrix(PrecMatr &G,const PrecMatr &Atr_inv, vector<PrecVect> &Rt_n) {

  G.resize(Nmoms,Nmoms);

  for(int i=0;i<Nmoms;i++) {
    for(int j=0;j<Nmoms;j++) {

      G(i,j) = Rt_n[i].transpose()*Atr_inv*Rt_n[j];

    }
  }

  return;
  
}

void Compute_covariance_matrix(PrecMatr &B,const PrecMatr &Atr, const distr_t_list &corr, int tmin, int tmax, string MODE) {

  
  PrecFloat norm= sqr(PrecFloat(corr.ave(0)));

  B.resize(tmax-tmin+1, tmax-tmin+1);

  for(int t=tmin;t<=tmax;t++) {
    for(int r=tmin;r<=tmax;r++) {
      if(t==r) {
	if(MODE =="TANT") B(t-tmin,r-tmin) =  (1/norm)*sqr( PrecFloat(corr.err(t)));
	else B(t-tmin,r-tmin) = Atr(t-tmin,r-tmin)*sqr(PrecFloat(corr.err(t)/corr.ave(t)));
      }
      else if( (r==t+1) || (t==r+1)) {
	if(COV_MATRIX_MODE == "NN") {
	   PrecFloat cov_t_t1 = corr.distr_list[t]%corr.distr_list[r];
	   if(MODE == "TANT") B(t-tmin,r-tmin) =  cov_t_t1/norm;
	   else B(t-tmin,r-tmin) = Atr(t-tmin,r-tmin)*PrecFloat(cov_t_t1/(corr.ave(t)*corr.ave(r)));
	}
	else {
	  B(t-tmin, r-tmin) = PrecFloat(0);
	}
      }
      else {
	B(t-tmin, r-tmin) = PrecFloat(0);
      }
    }
  }

}


void Get_optimal_lambda(const PrecMatr &Atr,const PrecMatr &B,const PrecVect &ft,const PrecVect &Rt,const PrecFloat & M2,const double &mean, const double &sigma, const double &Estart,  double& lambda_opt, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f, vector<PrecVect> Rt_n, const PrecVect &M_n ,const distr_t_list & corr,int T, int tmin, int tmax,const double mult,  string MODE, string curr_type, string SMEARING_FUNC, string CORR_NAME, string FLAV) {



  //create print directory
  boost::filesystem::create_directory("../data/spectral_reconstruction/smearing/lambda_stability");

  int MAX_Iters = 1000;

 
  string out_path = MODE+"_"+CORR_NAME+"_"+curr_type+"_"+SMEARING_FUNC+"_E*_"+to_string_with_precision(mean,3)+"_sigma_"+to_string_with_precision(sigma,3)+"_E0_"+to_string_with_precision(Estart,3)+"_T_"+to_string(T)+"_tmax_"+to_string(tmax);

  ofstream Print_R_at_lambda("../data/spectral_reconstruction/smearing/lambda_stability/"+out_path);
  
  cout<<"Finding optimal lambda* ..."<<endl;

   
  const auto A1=
    [&M2, &Atr, &ft](const PrecVect& gmin, PrecFloat lx) -> PrecFloat
    {
      PrecVect Atr_g = Atr*gmin;
      PrecFloat g_Atr_g = gmin.transpose()*Atr_g;
      PrecFloat ft_g = ft.transpose()*gmin;
      return (M2 + g_Atr_g -2*ft_g)/M2  ;
    };

  const auto B1=
    [&B, &M2, &MODE](const PrecVect& gmin, PrecFloat lx) -> PrecFloat
    {
      PrecVect B_g = B*gmin;
      PrecFloat g_B_g = gmin.transpose()*B_g;
      if(MODE=="SANF") g_B_g = g_B_g/M2; 
      return g_B_g ;
    };

  
  //bisection search
  PrecFloat l_low =0;
  PrecFloat l_up = 1;
  PrecFloat diff;
  PrecFloat l_start=1.0;
  int Nit=0;
  int Nit_Ag0=0;
  double lambda_Ag0;
  double lambda_balance;
  bool lambda_found_Ag_A0=false;
  bool lambda_balance_found=false;
  PrecFloat Ag_ov_A0_target=1e-4;
  if(FLAV=="fake") {
    if(MODE=="SANF") {
      if(mean < 0.2) Ag_ov_A0_target = 5e-3 ;
      else if(mean < 0.3) Ag_ov_A0_target = 5e-3   ;
      else if(mean < 0.4) Ag_ov_A0_target = 5e-3  ;
      else Ag_ov_A0_target = 5e-3   ;

    }
    if(MODE=="TANT") {
       if(mean < 0.2) Ag_ov_A0_target = 5e-3 ;
      else if(mean < 0.3) Ag_ov_A0_target = 5e-3   ;
      else if(mean < 0.4) Ag_ov_A0_target = 5e-3  ;
      else Ag_ov_A0_target = 5e-3   ;

    }
    
  }
  if(FLAV=="light") {
    Ag_ov_A0_target = (MODE=="TANT")?1.0e-3:1.0e-6;
  }
  else if (FLAV=="strange") {
    Ag_ov_A0_target = (MODE=="TANT")?1.0e-3:1.0e-6;
  }
  else if(FLAV=="charm") {
    Ag_ov_A0_target = (MODE=="TANT")?1.0e-4:1.0e-9;
  }

  cout<<"Finding lambda corresponding to A[g]/A[0] = "<<Ag_ov_A0_target<<endl;

  vector<PrecFloat> lambdas({0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85});
  double Nfs_L= 50;
  double Nfs_S= 100;
  double Nfs_SSS= 100;
  double Nfs_SSSS= 100;
  double Nfs_SSSSS= 100;
  double Nfs_LL = 100;
  double Nfs_SS=100;
  for(int i=0; i < Nfs_SSSSS; i++) lambdas.push_back( PrecFloat( 0.0000001 + i*0.000001/Nfs_SSSSS));
  for(int i=0; i < Nfs_SSSS; i++) lambdas.push_back( PrecFloat( 0.000001 + i*0.00001/Nfs_SSSS));
  for(int i=0; i < Nfs_SSS; i++) lambdas.push_back( PrecFloat( 0.00001 + i*0.0001/Nfs_SSS));
  for(int i=0; i < Nfs_SS; i++) lambdas.push_back( PrecFloat( 0.0001 + i*0.001/Nfs_S));
  for(int i=0; i < Nfs_S; i++) lambdas.push_back( PrecFloat( 0.001 + i*0.01/Nfs_S));
  for(int i=0;i< Nfs_L;i++) lambdas.push_back( PrecFloat( 0.90 + i*0.1/50));
  for(int i=0;i< Nfs_LL;i++) lambdas.push_back( PrecFloat( 0.99 + i*0.01/100));

  
  
  
  //bisection search for given A[g]/A[0]
  while( !lambda_found_Ag_A0 ) {

    PrecFloat lambda_mid;
    if(Nit_Ag0 < (signed)lambdas.size()) lambda_mid = lambdas[Nit_Ag0];
    else  lambda_mid =  (Nit_Ag0==(signed)lambdas.size())?l_start:(l_up+l_low)/2;
    PrecMatr C = Atr*(1-lambda_mid)/M2 + B*lambda_mid/(MODE=="SANF"?M2:1.0);
    PrecMatr C_inv = C.inverse();
    PrecVect ft_l = ft*(1-lambda_mid)/M2;
    PrecVect M_tilde_n;
    PrecMatr G_n;

      //Lagrangian multipliers
      if(Nmoms > 0) {
	Get_M_tilde_N(ft_l, C_inv, Rt_n, M_tilde_n);
	Get_G_matrix(G_n, C_inv, Rt_n);
      }
      
      
      //get matrix-vector product C_inv * ft_l

      const PrecVect C_inv_ft_l = C_inv*ft_l; 

      //get matrix-vector product C_inv * Rt

      const PrecVect C_inv_Rt = C_inv*Rt;
      
      //get scalar product Rt * ( C_inv * ft_l)

      const PrecFloat Rt_C_inv_ft_l = Rt.transpose()*C_inv_ft_l;

      //get scalar product Rt * ( C_inv * Rt)

      const PrecFloat Rt_C_inv_Rt = Rt.transpose()*C_inv_Rt;


      PrecVect gm = C_inv*ft_l;

      PrecVect pn;

      PrecVect M_n_diff = M_n - M_tilde_n;

      if(Nmoms > 0) pn = G_n.inverse()*M_n_diff;

      for(int n=0;n<Nmoms;n++) {
	PrecVect lmult_n = C_inv*Rt_n[n]; 
	gm = gm + lmult_n*pn(n);
      }

      PrecFloat A1_val = A1(gm, lambda_mid);
      PrecFloat B1_val = B1(gm, lambda_mid);
      PrecFloat W_val = (1-lambda_mid)*A1_val + lambda_mid*B1_val/(MODE=="SANF"?M2:1.0);

      cout.precision(10);
      //cout<<"lambda : "<<lambda_mid.get()<<" A[g]/A[0] : "<<A1_val<<endl;

      if(Nit_Ag0 >= MAX_Iters) crash("After "+to_string(Nit_Ag0)+" iterations, target A[g]/A[0] cannot be obtained for CORR: "+CORR_NAME+" , MODE: "+MODE+", CURR_TYPE: "+curr_type+". Actual A[g]/A[0]: "+to_string_with_precision( A1_val.get(), 10) );


      //##########################################################################################
      //compute anti-Laplace transform corresponding to lambda_mid:
      distr_t R_E_lambda;
      
      for(int ijack=0; ijack< corr.distr_list[1].size(); ijack++) {

	PrecFloat spec_lambda_d_jack=0;

	for(int t=tmin;t<=tmax;t++) spec_lambda_d_jack = spec_lambda_d_jack + gm(t-tmin)*corr.distr_list[t].distr[ijack]; 

	R_E_lambda.distr.push_back( spec_lambda_d_jack.get());
      }

      Print_R_at_lambda<<"lambda: "<<lambda_mid<<" A[g]/A[0]: "<<A1_val<<" B[g]: "<<B1_val<<" W[g]/A[g]/A[0]: "<<W_val/A1_val<<" R: "<<R_E_lambda.ave()<<" +- "<<R_E_lambda.err()<<" 0"<<endl;

      //##########################################################################################

      if(Nit_Ag0 >= (signed)lambdas.size()) {
   

      if(A1_val > Ag_ov_A0_target) { // lambda_mid is new l_low
          l_up =lambda_mid;
      }
      else { //lambda_mid is new l_up
	l_low = lambda_mid;
      }

      }
      
      lambda_Ag0= lambda_mid.get();
      Nit_Ag0++;
      if( (A1_val < 1.1*Ag_ov_A0_target && A1_val > 0.90*Ag_ov_A0_target) &&  (Nit_Ag0 >= (signed)lambdas.size() ) ) lambda_found_Ag_A0=true;
   
  }

  cout.precision(20);
  cout<<"lambda(A[g]/A[0] ="<<Ag_ov_A0_target<<") = : "<<lambda_Ag0<<endl;
  

  cout<<"Finding lambda from balance condition A=mult*B..."<<endl;
  l_up=1.0;
  l_low =0.0;
  diff= l_up-l_low;

 
  //bisection search for condition A = B
  while( !lambda_balance_found ) {

    if(Nit > MAX_Iters) crash("After "+to_string(Nit)+" iterations, balance condition A = mult*B cannot be obtained for CORR: "+CORR_NAME+" , MODE: "+MODE+", CURR_TYPE: "+curr_type+", mult = "+to_string_with_precision( mult, 8));

    //evaluate the minimum at midpoint
    PrecFloat lambda_mid = (Nit==0)?l_start:(l_up+l_low)/2;
    PrecMatr C = Atr*(1-lambda_mid)/M2 + B*lambda_mid/(MODE=="SANF"?M2:1);
    PrecMatr C_inv = C.inverse();
    PrecVect ft_l = ft*(1-lambda_mid)/M2;
    PrecVect M_tilde_n;
    PrecMatr G_n;

      //Lagrangian multipliers
      if(Nmoms > 0) {
	Get_M_tilde_N(ft_l, C_inv, Rt_n, M_tilde_n);
	Get_G_matrix(G_n, C_inv, Rt_n);
      }
      
      
      //get matrix-vector product C_inv * ft_l

      const PrecVect C_inv_ft_l = C_inv*ft_l; 

      //get matrix-vector product C_inv * Rt

      const PrecVect C_inv_Rt = C_inv*Rt;
      
      //get scalar product Rt * ( C_inv * ft_l)

      const PrecFloat Rt_C_inv_ft_l = Rt.transpose()*C_inv_ft_l;

      //get scalar product Rt * ( C_inv * Rt)

      const PrecFloat Rt_C_inv_Rt = Rt.transpose()*C_inv_Rt;


      PrecVect gm = C_inv*ft_l;

      PrecVect pn;

      PrecVect M_n_diff = M_n - M_tilde_n;

      if(Nmoms > 0) pn = G_n.inverse()*M_n_diff;

      for(int n=0;n<Nmoms;n++) {
	PrecVect lmult_n = C_inv*Rt_n[n]; 
	gm = gm + lmult_n*pn(n);
      }

      PrecFloat A1_val = A1(gm, lambda_mid);
      PrecFloat B1_val = B1(gm, lambda_mid);
      PrecFloat W_val = (1-lambda_mid)*A1_val + lambda_mid*B1_val/(MODE=="SANF"?M2:1.0);
      


    
      
 
     
      if(mult*B1_val > A1_val) { // lambda_mid is new l_low
	l_low =lambda_mid;
      }
      else { //lambda_mid is new l_up
	l_up = lambda_mid;
      }

      diff = l_up-l_low;
      lambda_balance= lambda_mid.get();
      Nit++;
      if(diff/(l_up+l_low) < 0.01) lambda_balance_found=true;


      //##########################################################################################
      //compute anti-Laplace transform corresponding to lambda_mid:
      distr_t R_E_lambda;
      
      for(int ijack=0; ijack< corr.distr_list[1].size(); ijack++) {

	PrecFloat spec_lambda_d_jack=0;

	for(int t=tmin;t<=tmax;t++) spec_lambda_d_jack = spec_lambda_d_jack + gm(t-tmin)*corr.distr_list[t].distr[ijack]; 

	R_E_lambda.distr.push_back( spec_lambda_d_jack.get());
      }

      Print_R_at_lambda<<"lambda: "<<lambda_mid<<" A[g]/A[0]: "<<A1_val<<" B[g]: "<<B1_val<<" W[g]/A[g]/A[0]: "<<W_val/A1_val<<" R: "<<R_E_lambda.ave()<<" +- "<<R_E_lambda.err()<<" "<<lambda_balance_found<<endl;

      //##########################################################################################
    
   
  }


  cout<<"lambda_opt = "<<lambda_balance<<endl;


  if(Use_balance_condition) lambda_opt =  lambda_balance;
  else lambda_opt = lambda_Ag0;

  Print_R_at_lambda.close();

  
  return;

}


distr_t Get_Laplace_transfo( double mean, double sigma, double Estart, int T, int tmax, int prec, string SMEARING_FUNC, const function<PrecFloat(const PrecFloat&, const PrecFloat&,const PrecFloat&,const PrecFloat&)> &f, const distr_t_list &corr, double &syst,const double mult, double& lambda_ret, string MODE, string cur_type, string CORR_NAME, string FLAV) {


  if(MODE != "TANT" && MODE != "SANF") crash("MODE: "+MODE+" not recognized");


  //create output directory
  boost::filesystem::create_directory("../data/spectral_reconstruction");
  boost::filesystem::create_directory("../data/spectral_reconstruction/smearing");

  PrecFloat::setDefaultPrecision(prec);

  PrecFloat s = sigma;
  PrecFloat m = mean;
  PrecFloat E0 = Estart;


  PrecMatr Atr,B;
  PrecVect ft, Rt;
  PrecFloat M2, N;


  //vectors for Lagrangian multiplier
  vector<PrecVect> Rt_n;
  PrecVect M_n;
  PrecVect M_tilde_n;
  PrecMatr G_n;

  //get vectors ft, Rt, and matrix Atr

 
  Get_ft(ft, E0, m, s, T, 1, tmax, SMEARING_FUNC, f);

  Get_Rt(Rt,E0, T, 1, tmax);
  
  Get_Atr(Atr, E0, T, 1, tmax);

  M2= Get_M2(m,s,E0, f);

  N= Get_norm_constraint(m,s,E0,SMEARING_FUNC, f);

 
  if(Nmoms > 0) {

     Get_Rt_up_to_N(E0, T, 1, tmax, Rt_n);
     Get_M_N(m,s,E0,f, M_n);
 
  }
  

  if(INCLUDE_ERRORS) Compute_covariance_matrix(B,Atr, corr,1,tmax, MODE);

  double lambda_opt= lambda;

  if(INCLUDE_ERRORS && FIND_OPTIMAL_LAMBDA) Get_optimal_lambda(Atr, B, ft, Rt, M2, mean, sigma, Estart, lambda_opt , f, Rt_n, M_n, corr, T , 1 , tmax, mult,  MODE, cur_type, SMEARING_FUNC,  CORR_NAME, FLAV);


    							         
  if(INCLUDE_ERRORS) {Atr = Atr*(1-lambda_opt)/M2 + B*lambda_opt/((MODE=="SANF")?M2:1);  ft= ft*(1-lambda_opt)/M2;       }


  //invert Atr

  const PrecMatr Atr_inv = Atr.inverse();


  
  //to compute Lagrangian multiplier

  if(Nmoms > 0) {
    Get_M_tilde_N(ft, Atr_inv, Rt_n, M_tilde_n);
    Get_G_matrix(G_n, Atr_inv, Rt_n);
  }
  

  //get matrix-vector product Atr_inv * ft


  const PrecVect Atr_inv_ft = Atr_inv*ft; 

   //get matrix-vector product Atr_inv * Rt

  const PrecVect Atr_inv_Rt = Atr_inv*Rt;

  //get scalar product Rt * ( Atr_inv * ft)

  const PrecFloat Rt_Atr_inv_ft = Rt.transpose()*Atr_inv_ft;


  //get scalar product Rt * ( Atr_inv * Rt)

  const PrecFloat Rt_Atr_inv_Rt = Rt.transpose()*Atr_inv_Rt;


  //get p as in Francesco's note

  // const PrecFloat p = (N- Rt_Atr_inv_ft)/Rt_Atr_inv_Rt;


  //get g(t) coefficient vector
  

  PrecVect g = Atr_inv*ft;

  //add Lagrangian multipliers

  PrecVect pn;

  PrecVect M_n_diff = M_n - M_tilde_n;

  if(Nmoms > 0) pn= G_n.inverse()*M_n_diff;

  for(int n=0;n <Nmoms;n++) {
    PrecVect lmult_n = Atr_inv*Rt_n[n];
    g = g + lmult_n*pn(n);
  }

  // cout<<"##### PRINTING G ######"<<endl;
  //cout<<g<<endl;







  //print reconstructed gaussian, exact Gaussian, diff

  //compute it from E0 up to 10 E* with a step size of sigma * step_size

  vector<double> Reco;
  vector<double> Exact;
  vector<double> Err;
  vector<double> Erg;


  int Npoints=  (int)(((m+10*s -E0)/(s*step_size)).get());

  for(int ip=0; ip<Npoints;ip++) {
    PrecVect bt;
    PrecFloat E = E0 + (ip*s)*step_size;
    Get_bt(bt, E, T, 1, tmax);  
    const PrecFloat reco_result = g.transpose()*bt;
    const PrecFloat exact_result = Get_exact_func(E,m,s,E0, SMEARING_FUNC, f);
    Reco.push_back( reco_result.get());
    Exact.push_back( exact_result.get());
    Err.push_back(  (exact_result -reco_result).get());
    Erg.push_back(E.get());    
  }

  
  //print to file
  Print_To_File({}, { Erg, Reco, Exact, Err}, "../data/spectral_reconstruction/smearing/"+MODE+"_"+CORR_NAME+"_"+cur_type+"_"+SMEARING_FUNC+"_func_E*_"+to_string_with_precision(mean,3)+"_sigma_"+to_string_with_precision(sigma,3)+"_E0_"+to_string_with_precision(E0.get(),3)+"_T_"+to_string(T)+"_tmax_"+to_string(tmax), "", "#id  E   reco    exact    error");

  distr_t Spec_dens_at_E_star; //Uses Jackknife distr by default


  if(INCLUDE_ERRORS) {

    //build the spectral density at E*
    

    for(int ijack=0; ijack< corr.distr_list[1].size(); ijack++) {

      PrecFloat spec_d_jack=0;

      for(int t=1;t<=tmax;t++) spec_d_jack = spec_d_jack + g(t-1)*corr.distr_list[t].distr[ijack]; 

      Spec_dens_at_E_star.distr.push_back( spec_d_jack.get());
    }


    
  //print coefficient
  ofstream PrintCoeff("../data/spectral_reconstruction/smearing/"+MODE+"_"+CORR_NAME+"_"+cur_type+"_"+SMEARING_FUNC+"_coeff_E*"+to_string_with_precision(mean,3)+"_sigma_"+to_string_with_precision(sigma,3)+"_E0_"+to_string_with_precision(E0.get(),3)+"_T_"+to_string(T)+"_tmax_"+to_string(tmax));
  if(!INCLUDE_ERRORS) {
  PrintCoeff<<g<<endl;
  }
  else {
    for(int t=1;t<=tmax;t++) PrintCoeff<<g(t-1)<<"\t"<<g(t-1)*corr.ave(t)<<"\t"<<Spec_dens_at_E_star.ave()<<"\t"<<Spec_dens_at_E_star.err()<<endl;
  }
  PrintCoeff.close();
  


    //evaluate the smeared function at the peak
    PrecVect bt;
    Get_bt(bt, m, T, 1, tmax);
    const PrecFloat reco_at_peak= g.transpose()*bt;
    
    
    syst = 0.68*abs(Spec_dens_at_E_star.ave())*abs(  1 -  reco_at_peak.get()/Get_exact_func(m,m,s,E0, SMEARING_FUNC, f).get());
    
  }

  //set lambda_ret to lambda_optimal

  lambda_ret= lambda_opt;
  
  return Spec_dens_at_E_star;  
  
}
