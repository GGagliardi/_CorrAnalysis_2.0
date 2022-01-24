#include "../include/Spectral.h"


const double step_size = 0.1; //in units of sigma



using namespace std;



PrecFloat Get_exact_gauss(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0) {
  PrecFloat e = exp( -0.5*(E-m)*(E-m)/(s*s));
  PrecFloat norm= s*( 1.0 + erf( (m-E0)/(s*sqrt(PrecFloat(2)))))*sqrt(precPi()/PrecFloat(2)) ;
  return e/norm;
}


PrecFloat Get_exact_gaussE2(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0) {
  return Get_exact_gauss(E, m, s, E0)*E*E;
}

PrecFloat Get_gaussE2_norm(const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0) {

  PrecFloat a1 = exp( -sqr(m-E0)/(s*s*2))*(E0+m)*sqr(s);
  //PrecFloat a2 = sqrt( precPi()/PrecFloat(2))*s*(sqr(m)+sqr(s))*( 1 +erf( (m-E0)/(s*sqrt(PrecFloat(2)))) );
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

PrecFloat Get_exact_func(const PrecFloat &E,const PrecFloat &m,const PrecFloat &s,const PrecFloat &E0, string SMEARING_FUNC) {

  PrecFloat res;

   if(SMEARING_FUNC == "GAUSSIAN") { res = Get_exact_gauss(E, m, s, E0); }
   else if(SMEARING_FUNC == "GAUSSIAN_E2") {     res = Get_exact_gaussE2(E, m, s, E0);   }
   else if(SMEARING_FUNC == "LEGO") {  res = Get_exact_lego(E, m, s, E0);  } 
   else if(SMEARING_FUNC == "LEGO_E2") {  res = Get_exact_legoE2(E, m, s, E0); }  
   else  crash("In Get_exact_func cannot find SMEARING_FUNC: "+SMEARING_FUNC);

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

PrecFloat aE0(const PrecFloat &E0, int t) { return exp(-E0*t)/t;}

void Get_Atr(PrecMatr& Atr, const PrecFloat &E0, int T, int tmin, int tmax)  {

  Atr.resize(tmax-tmin+1, tmax-tmin+1);

  for(int t=tmin;t<= tmax; t++)
    for(int r=tmin; r<= tmax; r++) Atr(t-tmin,r-tmin) = aE0(E0, t+r) + aE0(E0, T - t +r) + aE0(E0, T+t -r) + aE0(E0, 2*T -t -r);
    
  
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

void Get_ft(PrecVect& ft, const PrecFloat &E0, const PrecFloat &m, const PrecFloat &s, int T, int tmin, int tmax, string SMEARING_FUNC) {

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

  else  crash("In Get_ft cannot find SMEARING_FUNC: "+SMEARING_FUNC);


  return;
}

void Get_bt(PrecVect& bt,const PrecFloat &E,  int T, int tmin, int tmax) {

  bt.resize(tmax-tmin+1);

  for(int t=tmin;t<=tmax;t++) bt(t-tmin) = BaseFunc(E, t, T);

  return;


}

PrecFloat Get_norm_constraint(PrecFloat &m, PrecFloat &s, PrecFloat &E0, string SMEARING_FUNC) {

  PrecFloat res=1;

  if(SMEARING_FUNC == "GAUSSIAN_E2") {
     res= Get_gaussE2_norm(m, s, E0);
   }
  else if(SMEARING_FUNC == "LEGO_E2") {
    return Get_legoE2_norm(m,s,E0);
  }

  cout<<"norm: "<<res<<endl;
  return res;
}


void Get_Laplace_transfo( double mean, double sigma, double Estart, int T, int tmax, int prec, string SMEARING_FUNC) {


  PrecFloat::setDefaultPrecision(prec);

  PrecFloat s = sigma;
  PrecFloat m = mean;
  PrecFloat E0 = Estart;

  PrecMatr Atr;
  PrecVect ft, Rt;

  //get vectors ft, Rt, and matrix Atr


  Get_ft(ft, E0, m, s, T, 1, tmax, SMEARING_FUNC);    
  
  Get_Rt(Rt,E0, T, 1, tmax);
  
  Get_Atr(Atr, E0, T, 1, tmax);

  //invert Atr

  const PrecMatr Atr_inv = Atr.inverse();

  //get matrix-vector product Atr_inv * ft


  const PrecVect Atr_inv_ft = Atr_inv*ft; 

   //get matrix-vector product Atr_inv * Rt

  const PrecVect Atr_inv_Rt = Atr_inv*Rt;

  //get scalar product Rt * ( Atr_inv * ft)

  const PrecFloat Rt_Atr_inv_ft = Rt.transpose()*Atr_inv_ft;


  //get scalar product Rt * ( Atr_inv * Rt)

  const PrecFloat Rt_Atr_inv_Rt = Rt.transpose()*Atr_inv_Rt;


  //get p as in Francesco's note

  const PrecFloat p = (Get_norm_constraint(m,s,E0,SMEARING_FUNC) - Rt_Atr_inv_ft)/Rt_Atr_inv_Rt;

  cout<<"p: "<<p<<endl;


  //get g(t) coefficient vector


  const PrecVect g = Atr_inv*( ft + p*Rt);

  cout<<"##### PRINTING G ######"<<endl;
  cout<<g<<endl;


  //create output directory
  boost::filesystem::create_directory("../data/spectral_reconstruction");
  boost::filesystem::create_directory("../data/spectral_reconstruction/smearing");


  //print coefficient
  ofstream PrintCoeff("../data/spectral_reconstruction/smearing/"+SMEARING_FUNC+"_coeff_"+to_string_with_precision(mean,3)+"_sigma_"+to_string_with_precision(sigma,3)+"_E0_"+to_string_with_precision(E0.get(),3)+"_T_"+to_string(T)+"_tmax_"+to_string(tmax));
  PrintCoeff<<g<<endl;
  PrintCoeff.close();


  //print reconstructed gaussian, exact Gaussian, diff

  //compute it from E0 up to 10 E* with a step size of sigma * step_size

  vector<double> Reco;
  vector<double> Exact;
  vector<double> Err;
  vector<double> Erg;


  int Npoints=  (int)(((5*m -E0)/(s*step_size)).get());

  for(int ip=0; ip<Npoints;ip++) {
    PrecVect bt;
    PrecFloat E = E0 + (ip*s)*step_size;
    Get_bt(bt, E, T, 1, tmax);  
    const PrecFloat reco_result = g.transpose()*bt;
    const PrecFloat exact_result = Get_exact_func(E,m,s,E0, SMEARING_FUNC);
    Reco.push_back( reco_result.get());
    Exact.push_back( exact_result.get());
    Err.push_back(  (exact_result -reco_result).get());
    Erg.push_back(E.get());    
  }

  
  //print to file
  Print_To_File({}, { Erg, Reco, Exact, Err}, "../data/spectral_reconstruction/smearing/"+SMEARING_FUNC+"_func_"+to_string_with_precision(mean,3)+"_sigma_"+to_string_with_precision(sigma,3)+"_E0_"+to_string_with_precision(E0.get(),3)+"_T_"+to_string(T)+"_tmax_"+to_string(tmax), "", "#id  E   reco    exact    error");
  
  






  return;  

}
