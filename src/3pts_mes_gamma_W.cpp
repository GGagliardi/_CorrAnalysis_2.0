#include "../include/3pts_mes_gamma_W.h"



using namespace std;
namespace plt = matplotlibcpp;

const double M2PiPhys=pow(0.135,2);
const double alpha = 1/137.04;
const double e2 = alpha*4.0*M_PI;
const int Nbranches = 8;
const int Nboots= 100;
const bool UseJack=1;
const int nboots=150;
const int Njacks=15;
int sm_lev=1;
const double e_f1 = 2.0/3.0; //electric charge of u-type quark
const double e_f2 = -1.0/3.0; //electric charge of d-type quark
bool Include_k0_noise=0;
bool Determine_contaminations=0;
int verbosity=0; //used in Fit Contaminations to print infos about minimization
const string Meson="K";

class ExpFit {

 
public:
  ExpFit(const Vfloat &par)  { if((signed)par.size() != 5) crash("par size != 5");
    a1=par[0];
    a2=par[1];
    a3=par[2];
    a4=par[3];
    a5=par[4];

  }

  double a1;
  double a2;
  double a3;
  double a4;
  double a5;
};

class Val{
public:
  Val() : t(0), meas(0), err(0) {}
  double t;
  double meas;
  double err;
};

void Plot_form_factors(string W, distr_t_list& F, distr_t& F_fit, int Tmin, int Tmax, int Nt, string Ens_tag, double xg, double offsh, int smearing) {

  boost::filesystem::create_directory("../plots");
  boost::filesystem::create_directory("../plots/form_factors");
  boost::filesystem::create_directory("../plots/form_factors/"+Meson);
  plt::clf();
  plt::figure_size(1500,1000);
  plt::xlabel("$t/a$");
  if(W=="V") plt::ylabel("$F_{V}$");
  else plt::ylabel("$F_{A}$");
  plt::xlim(3, Nt/2 -3);
  double y_min=1e10, y_max= -1e10;
  for(int t=3;t<=Nt/2 -3;t++) { y_min = min( y_min, F.ave()[t] - 1.2*F.err()[t]); y_max= max( y_max, F.ave()[t] + 1.2*F.err()[t]);}  
  plt::ylim(y_min, y_max);
  plt::grid(4);
  Vfloat TT;
  Vfloat TT_tot;
  for(int t=0;t<Nt;t++) TT_tot.push_back(t);
  //generate data_point for fit;
  Vfloat F_val(Tmax - Tmin +1, F_fit.ave());
  Vfloat F_err(Tmax - Tmin +1, F_fit.err());
  for(int t=Tmin;t<=Tmax;t++)  TT.push_back(t);
  plt::errorbar(TT_tot, F.ave(), F.err(), { {"c", "black"}, {"marker", "."} , {"ls" , "-"}, {"label", "xg="+to_string_with_precision(xg,4)}});
  plt::errorbar(TT, F_val, F_err, { {"c", "red"}, {"marker", "."}, {"ls", "-"}, {"label", "fitted value"}});
  // Enable legend.
  plt::legend();
  string figure_path;
  if(W=="V") figure_path= "../plots/form_factors/"+Meson+"/FV_"+Ens_tag+"_xg_"+to_string_with_precision(xg,4)+"_offsh_"+to_string_with_precision(offsh,4)+"_smlev_"+to_string(smearing)+"_k0_noise_"+to_string(Include_k0_noise)+".png";
  else figure_path= "../plots/form_factors/"+Meson+"/FA_"+Ens_tag+"_xg_"+to_string_with_precision(xg,4)+"_offsh_"+to_string_with_precision(offsh,4)+"_smlev_"+to_string(smearing)+"_k0_noise_"+to_string(Include_k0_noise)+".png";
  plt::save(figure_path.c_str());
  plt::close();
  return;
}

void Fit_contaminations(string W, distr_t_list& F, distr_t &F_fit, string Ens_tag, int Nt, double xg, double offsh, int smearing) {

  

  double T_min, T_max_sx, T_min_dx, T_max;
  T_min=2;
  if(W=="RV" || W=="FV") T_max= Nt/2 - 3;
  else T_max=Nt/2 -2;
  T_max_sx= (double)min(14, Nt/4);
  T_min_dx= T_max_sx+1;
  int ndof= (T_max -T_min - 5);
  //choose the t-interval
  auto ansatz_sx= [&](const Vfloat& par, double t) -> double { return par[0] +par[1]*exp(-1*par[2]*t);};
  auto ansatz_dx= [&](const Vfloat& par, double t) -> double { return par[0] +par[1]*exp(-1*par[2]*(Nt/2-t));};
  auto ansatz_tot=[&](const Vfloat& par, double t) -> double { return par[0] +par[1]*exp(-1*par[2]*t) + par[3]*exp(-1*par[4]*(Nt/2-t));};
  auto ansatz = [=](const ExpFit& par, const Val& val) -> double {return par.a1 + par.a2*exp(-1*par.a3*val.t) +  par.a4*exp(-1*par.a5*(Nt/2-val.t));};
  Vfloat T_dx, T_sx, Y_dx, Y_sx, Y_err_dx, Y_err_sx, T_dx_sx, Y_dx_sx, Y_err_dx_sx;
  for(int t=T_min;t<=T_max_sx;t++) {T_sx.push_back(t); Y_sx.push_back( F.ave()[t]); Y_err_sx.push_back( F.err()[t]);}
  for(int t=T_min_dx;t<=T_max;t++) {T_dx.push_back(t); Y_dx.push_back( F.ave()[t]); Y_err_dx.push_back( F.err()[t]);}
  for(int t=T_min;t<=T_max;t++) {T_dx_sx.push_back(t); Y_dx_sx.push_back( F.ave()[t]); Y_err_dx_sx.push_back( F.err()[t]);}


  double a1_guess, a1_guess_dx, a2_guess, a3_guess, a4_guess,a5_guess;
  pair<double,double> fit_interval_a3;
  pair<double, double> fit_interval_a5;
  int sign_a2, sign_a4;


  //choose guess intervals depending on the observable
  //###############################################################
  fit_interval_a3= make_pair(0.01, 1.0);
  fit_interval_a5= make_pair(0.01, 1.0);
  if( (W=="RV" || W=="FV")) {
    if(smearing) {
      if(W=="FV") {
	sign_a2=-1;
	sign_a4=1;
      }
      else {
	sign_a2=1;
	sign_a4=-1;
      }
    }
    else {
      if(W=="FV") {
	sign_a2=1;
	sign_a4=1;
      }
      else {
	sign_a2=-1;
	sign_a4=1;
      }

    }
  }
  else {
    if(W=="FA") {
      sign_a2=1;
      sign_a4=-1;
    }
    else {
      if(smearing) {
	sign_a2=1;
	sign_a4=-1;
      }
      else {
	sign_a2=1;
	sign_a4=-1;
      }
    }
  }
  //################################################################

  
  //SX_FIT
  UniformMersenne  A2(6545534,log(0.001), log(10.0));
  UniformMersenne A4(87756756,log(0.001), log(10.0));
  UniformMersenne A3(659878,fit_interval_a3.first, fit_interval_a3.second);
  UniformMersenne A5(546353,fit_interval_a5.first, fit_interval_a5.second);
  int its_sx=0;
  bool sx_min_found=false;
  while(!sx_min_found) {
    its_sx++;
    T_fit fit_sx(T_sx, Y_sx, Y_err_sx);
    fit_sx.ansatz = ansatz_sx;
    int ndof = (T_max_sx-T_min +1 - 3);
    double sx_1 = sign_a2*exp( A2());
    double sx_2 = A3();
    fit_sx.add_pars({a1_guess, sx_1, sx_2});
    fit_sx.add_par_errs({F_fit.err(), sx_1/10, sx_2/10});
    fit_t_res sx_fit_res = fit_sx.fit();
    if(sx_fit_res.chi2/ndof < 1) { cout<<"ch2_sx for OBS: "<<W<<" xg: "<<xg<<" ch2: "<<sx_fit_res.chi2/ndof<<"  ITs: "<<its_sx<<endl; a1_guess= sx_fit_res.pars[0];a2_guess= sx_fit_res.pars[1]; a3_guess=sx_fit_res.pars[2]; sx_min_found=true;}
      if (its_sx % 100 == 0) cout<<"IT(sx): "<<its_sx<<" OBS: "<<W<<" ,chi2/ndof: "<<sx_fit_res.chi2/ndof<<endl;
    if(sx_fit_res.chi2/ndof < 4 && its_sx > 10000) {
      cout<<"ch2_sx for OBS: "<<W<<" xg: "<<xg<<" ch2: "<<sx_fit_res.chi2/ndof<<" ITs: "<<its_sx<<endl;
      cout<<"do you want to continue?"<<endl;
      cin >> sx_min_found;
      if(sx_min_found) { a1_guess= sx_fit_res.pars[0]; a2_guess=sx_fit_res.pars[1]; a3_guess=sx_fit_res.pars[2];}
    }
    if(its_sx > 20000) {sign_a2 *= -1;}
    if(its_sx > 30000) return;  //not able to perform the fit
  }

 

  //DX_FIT
  bool dx_min_found=false;
  int its_dx=0;
  while(!dx_min_found) {
    its_dx++;
    T_fit fit_dx(T_dx, Y_dx, Y_err_dx);
    fit_dx.ansatz= ansatz_dx;
    int ndof = (T_max-T_min_dx +1 - 3);
    double dx_1= sign_a4*exp(A4());
    double dx_2= A5();
    fit_dx.add_pars({F_fit.ave(), dx_1, dx_2});
    fit_dx.add_par_errs({F_fit.err(), dx_1/10, dx_2/10});
    fit_t_res dx_fit_res = fit_dx.fit();
    if(dx_fit_res.chi2/ndof < 1) { cout<<"ch2_dx for OBS: "<<W<<" xg: "<<xg<<" ch2: "<<dx_fit_res.chi2/ndof<<" ITs: "<<its_dx<<endl; a4_guess= dx_fit_res.pars[1]; a5_guess=dx_fit_res.pars[2]; dx_min_found=true;}
      if (its_dx % 100 == 0) cout<<"IT(dx): "<<its_dx<<" OBS: "<<W<<" ,chi2/ndof: "<<dx_fit_res.chi2/ndof<<endl;
    if(dx_fit_res.chi2/ndof < 4 && its_dx > 10000) {
      cout<<"ch2_dx for OBS: "<<W<<" xg: "<<xg<<" ch2: "<<dx_fit_res.chi2/ndof<<" ITs: "<<its_dx<<endl;
      cout<<"do you want to continue?"<<endl;
      cin >> dx_min_found;
      if(dx_min_found) { a1_guess_dx=dx_fit_res.pars[0] ;a4_guess=dx_fit_res.pars[1]; a5_guess=dx_fit_res.pars[2];}
    
    }
    if(its_dx > 20000) {sign_a4*= -1;}
    if(its_dx > 30000) return; //not able to perform the fit 
  }
   //fit_tot
  //############################################################################
  /*
  T_fit fit_dx_sx(T_dx_sx,Y_dx_sx, Y_err_dx_sx);
  fit_dx_sx.ansatz= ansatz_tot;
  fit_dx_sx.add_pars({a1_guess, a2_guess, a3_guess, a4_guess, a5_guess});
  fit_t_res dx_sx_fit_res = fit_dx_sx.fit();
  a1_guess = dx_sx_fit_res.pars[0];
  a2_guess = dx_sx_fit_res.pars[1];
  a3_guess = dx_sx_fit_res.pars[2];
  a4_guess = dx_sx_fit_res.pars[3];
  a5_guess = dx_sx_fit_res.pars[4];
  */
  //############################################################################

 
  //BOOTSTRAP FIT
  //##############################################################################
  bootstrap_fit<ExpFit, Val> fit(Nboots);
  fit.measurement = [=](const ExpFit& par, const Val& val) -> double { return val.meas;};
  fit.error = [=](const ExpFit& par, const Val& val) -> double {return val.err;};
  fit.ansatz =  ansatz;
  fit.Set_number_of_measurements(T_max - T_min +1);
  fit.Set_verbosity(verbosity);
  
  
  cout<<"a1_g: "<<a1_guess<<"  "<<a1_guess_dx<<endl;
  cout<<"a2_g: "<<a2_guess<<endl;
  cout<<"a3_g: "<<a3_guess<<endl;
  cout<<"a4_g: "<<a4_guess<<endl;
  cout<<"a5_g: "<<a5_guess<<endl;
  fit.Add_par("a1", a1_guess, fabs(a1_guess)/30);
  fit.Add_par("a2", a2_guess, fabs(a2_guess)/30);
  fit.Add_par("a3", a3_guess, fabs(a3_guess)/30);
  fit.Add_par("a4", a4_guess, fabs(a4_guess)/30);
  fit.Add_par("a5", a5_guess, fabs(a5_guess)/30);
  

  //prepare gaussian number generator
  GaussianMersenne GM(4355);
  double F_resc=3.0;
 
  //generate bootstrap sample
  for(int iboot=0;iboot<Nboots;iboot++) {
    fit.ib =&iboot;
    for(int i=T_min;i<=T_max;i++) {
      double meas= F.ave()[i];
      double err = F.err()[i];
      Val X;
      X.t= (double)i;
      X.meas=meas+ GM()*err/F_resc;
      X.err=err;
      fit.Append_to_input_par(X);
    }
  }

          

  //Perform bootstrap fit
  boot_fit_data<ExpFit> Fit_output= fit.Perform_bootstrap_fit();



  //retrieve output params
  //##############################################
  
  double ave_ch2=0;
  Vfloat a1_b, a2_b, a3_b, a4_b, a5_b;
  int iboot_min=5;
  for(int iboot=iboot_min;iboot<Nboots;iboot++) {
    a1_b.push_back(Fit_output.par[iboot].a1);
    a2_b.push_back(Fit_output.par[iboot].a2);
    a3_b.push_back(Fit_output.par[iboot].a3);
    a4_b.push_back(Fit_output.par[iboot].a4);
    a5_b.push_back(Fit_output.par[iboot].a5);
    ave_ch2+= Fit_output.chi2[iboot]/(Nboots-iboot_min);
  }


  //print output parameters
  ofstream Print_out("../data/form_factors/"+Meson+"/Contaminations/contaminations_"+W+"_"+Ens_tag+".dat", ofstream::app);
  Print_out<<xg<<setw(20)<<offsh<<setw(20)<<smearing<<setw(20)<<Boot_ave(a1_b)<<setw(20)<<F_resc*Boot_err(a1_b)<<setw(20)<<Boot_ave(a2_b)<<setw(20)<<F_resc*Boot_err(a2_b)<<setw(20)<<Boot_ave(a3_b)<<setw(20)<<F_resc*Boot_err(a3_b)<<setw(20)<<Boot_ave(a5_b)<<setw(20)<<F_resc*Boot_err(a5_b)<<setw(20)<<ave_ch2/ndof<<endl;
  Print_out.close();


  //print fitting function
  boost::filesystem::create_directory("../data/form_factors/"+Meson+"/Contaminations/fit_func");
  ofstream Print_func("../data/form_factors/"+Meson+"/Contaminations/fit_func/"+W+"_"+Ens_tag+"_xg_"+to_string_with_precision(xg,4)+"_offsh_"+to_string_with_precision(offsh,4)+"_smlev_"+to_string(smearing)+"_k0_noise_"+to_string(Include_k0_noise)+".dat");

  int tsteps=500;

  for(int is=0;is<tsteps;is++) {

    Vfloat fitted_t;
    double t= ((double)is)*(Nt/(2.0*(double)tsteps));
    Val X;
    X.t = t;
    for(int iboot=iboot_min;iboot<Nboots;iboot++) fitted_t.push_back( ansatz(Fit_output.par[iboot], X  ));
    Print_func<<t<<setw(20)<<Boot_ave(fitted_t)<<setw(20)<<3*Boot_err(fitted_t)<<endl;
  }


  Print_func.close();
  //###############################################################################

  return;
}

void Get_Tmin_Tmax(string corr_type, string Ens_tag, CorrAnalysis &corr, double xg, string W) {

  if(corr_type=="2pt") {
    corr.Tmax = corr.Nt/2 -4;
    if(Ens_tag.substr(0,1) =="A" || Ens_tag.substr(0,1) == "C" || Ens_tag.substr(0,1) == "B" || Ens_tag.substr(0,1)== "G" || Ens_tag.substr(0,1)=="D" || Ens_tag.substr(0,1) =="E" ) corr.Tmin = 13;
    else if(Ens_tag.substr(0,1) =="B") corr.Tmin = 16;
    else corr.Tmin = 19;
  }

  if(corr_type=="3pt") {
    if(Ens_tag.substr(0,6)=="A40.32") {
      if(W=="A") { 
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if( xg > 0.09 && xg < 0.11)  {corr.Tmax = 17; corr.Tmin= 13;}
	else if( xg > 0.19 && xg < 0.21)  {corr.Tmax = 17; corr.Tmin=13;}
        else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 13;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 16; corr.Tmin= 13;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 15; corr.Tmin= 11;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 15; corr.Tmin= 11;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      else if (W=="V") {
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if (xg >0.09 && xg < 0.11) { corr.Tmax = 29; corr.Tmin= 23;}
	else if( xg > 0.19 && xg < 0.21) { corr.Tmax = 27; corr.Tmin= 22;}
	else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 24; corr.Tmin= 16;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 25; corr.Tmin= 12;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 23; corr.Tmin= 12;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 17; corr.Tmin= 13;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      
    }
    else  if(Ens_tag.substr(0,6)=="B40.32") {
      if(W=="A") { 
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if( xg > 0.09 && xg < 0.11)  {corr.Tmax = 17; corr.Tmin= 13;}
	else if( xg > 0.19 && xg < 0.21)  {corr.Tmax = 16; corr.Tmin=12;}
        else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 13;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 16; corr.Tmin= 13;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 15; corr.Tmin= 11;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 15; corr.Tmin= 11;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      else if (W=="V") {
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if (xg >0.09 && xg < 0.11) { corr.Tmax = 20; corr.Tmin= 15;}
	else if( xg > 0.19 && xg < 0.21) { corr.Tmax = 20; corr.Tmin= 15;}
	else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 20; corr.Tmin= 12;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 16; corr.Tmin= 12;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 16; corr.Tmin= 11;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 16; corr.Tmin= 11;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      
    }
    else if(Ens_tag.substr(0,6)=="A40.24" || Ens_tag.substr(0,6)=="B40.24" || Ens_tag.substr(0,6)=="C40.24") {
      if(W=="A") {
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if( xg > 0.19 && xg < 0.21)  {corr.Tmax = 16; corr.Tmin= 9;}
        else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 9;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 14; corr.Tmin= 8;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 13; corr.Tmin= 8;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 14; corr.Tmin= 8;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      else if (W=="V") {
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 23; corr.Tmin= 18;}
	else if( xg > 0.09 && xg < 0.11) { corr.Tmax = 18; corr.Tmin= 14;}
	else if( xg > 0.19 && xg < 0.21) { corr.Tmax = 18; corr.Tmin= 14;}
	else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 17; corr.Tmin= 12;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 13; corr.Tmin= 12;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 13; corr.Tmin= 10;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 13; corr.Tmin= 10;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      else crash("cannot find the ensemble");
    }
    else if(Ens_tag.substr(0,6)=="A40.40") {
      if(W=="A") {
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 25; corr.Tmin= 18;}
	else if( xg > 0.09 && xg < 0.11)  {corr.Tmax = 25; corr.Tmin= 18;}
	else if( xg > 0.19 && xg < 0.21)  {corr.Tmax = 20; corr.Tmin= 15;}
        else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 19; corr.Tmin= 14;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 19; corr.Tmin= 14;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 17; corr.Tmin= 12;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 16; corr.Tmin= 11;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      else if (W=="V") {
	if( xg > -0.05 && xg < 0.05)  {corr.Tmax = 24; corr.Tmin= 18;}
	else if( xg > 0.09 && xg < 0.11)  {corr.Tmax = 15; corr.Tmin= 9;}
	else if( xg > 0.19 && xg < 0.21) { corr.Tmax = 19; corr.Tmin= 12;}
	else if (xg > 0.39 && xg < 0.41)  {corr.Tmax = 22; corr.Tmin= 16;}
	else if( xg > 0.59 && xg < 0.61)  {corr.Tmax = 20; corr.Tmin= 14;}
	else if( xg > 0.79 && xg < 0.81)  {corr.Tmax = 19; corr.Tmin= 13;}
	else if( xg > 0.9 && xg < 1.1)  {corr.Tmax = 21; corr.Tmin= 12;}
	else {corr.Tmin= 9; corr.Tmax=16;}
      }
      else crash("cannot find the ensemble");
    }
    else if(Ens_tag.substr(0,6)=="A40.48") {
      if(W=="A") { corr.Tmax=23; corr.Tmin=18;}
      else if(W=="V") {corr.Tmax= 24; corr.Tmin=18;}
    }

    
    else { corr.Tmin= 9; corr.Tmax=16;}
  }

  return;

}


void Add_to_mom_list(pt3_momenta_list &M, struct header_virph &header, double& L) {

  vector<pt3_momenta> A;


  for(auto & c: header.comb) A.emplace_back(Vfloat{c.th0[0], c.th0[1], c.th0[2]},Vfloat{c.ths[0], c.ths[1], c.ths[2]},Vfloat{c.tht[0], c.tht[1], c.tht[2] }, c.mu1, c.mu2, c.off,L, header.tmax); 
  M.mom.push_back(A);
  return;
}


void Add_to_mom_list(pt2_momenta_list &M, struct header_virph &header, double& L) {

  vector<pt2_momenta> A;

  for(auto & c: header.comb) A.emplace_back( Vfloat{c.th0[0], c.th0[1], c.th0[2]}, Vfloat{c.ths[0], c.ths[1], c.ths[2]}, c.mu1, c.mu2, L, c.i0, c.is);

  M.mom.push_back(A);

  return;
}


distr_t_list V_ave_unpolarized(vector<vector<distr_t_list>>& distr_mom_k, vector<vector<distr_t_list>>& distr_mom_0, distr_t& Meson_mass, pt3_momenta& Mom,int twall) {

  VVfloat e(4);
  e[0] = {0.0, 0.0, 0.0, 0.0};
  e[3] = {0.0, 0.0, 0.0, 0.0};
  e[1] = {0.0, -1.0/sqrt(2), -1.0/sqrt(2), 0.0};
  e[2] = {0.0, 1.0/sqrt(2), -1.0/sqrt(2), 0.0};

  



  if(distr_mom_k.size() != 4 || distr_mom_0.size() != 4) crash("V_ave_unpolarized called with a v<v<distr_t_list>> with size != alphas=4");
  int size= distr_mom_k[0][0].size();
  if(size==0) crash("In V_ave_unpolarized, distr_mom_k[0][0] has size zero, exiting....");
  int distr_size = distr_mom_k[0][0].distr_list[0].size();
  if(distr_size==0) crash("In V_ave_unpolarized, distr_mom_k[0][0] has zero distr_size, exiting....");


  distr_t_list V_ave(UseJack, size, distr_size);

  Vfloat exp_Nt_t;
  double Egamma=Mom.Egamma();
  double Egamma_T=sinh(Egamma)*(1-exp(-Egamma*Mom.Nt()));
  for(int t=0; t<Mom.Nt();t++) exp_Nt_t.push_back( exp( fabs((twall-t))*Egamma));
  //for(int t=0;t<Mom.Nt(); t++) exp_Nt_t.push_back( exp( ((Mom.Nt()/2.0)-t)*Mom.Egamma()) - exp((t- Mom.Nt()/2.0)*Mom.Egamma()));

  bool is_mom_zero=false;
  Vfloat theta_t_old(Mom.Theta(2));
  if(Mom.k_mod() < eps(5)) { is_mom_zero=true; Vfloat theta_t = {0.0, 0.0, 1.0}; Mom.Set_theta(theta_t, 2);}

  
  
  for(int alpha=1; alpha <= 2; alpha ++) {
   
    if(distr_mom_k[alpha].size() != 4 || distr_mom_0[alpha].size() != 4) crash("V_ave_unpolarized called with a v<v<distr_t_list>> having an element with size != 4");
    for(int mu=0; mu<4;mu++) {

      if(mu != alpha) {
    
	distr_t_list den1(3, Meson_mass);
	distr_t_list den2(3, Meson_mass);
	
	
	den1 = den1*Multiply_vector_by_scalar(external_prod(slicing(e[1],1,3),Mom.k()),-1.0) + Multiply_vector_by_scalar(external_prod(slicing(e[1],1,3), Mom.p()), Mom.Egamma());
    
	den2 = den2*Multiply_vector_by_scalar(external_prod(slicing(e[2],1,3),Mom.k()), -1.0) + Multiply_vector_by_scalar(external_prod(slicing(e[2],1,3), Mom.p()), Mom.Egamma());
    
	V_ave = V_ave+ e[1][mu]*(distr_mom_k[alpha][mu]*exp_Nt_t - ((double)Include_k0_noise)*distr_mom_0[alpha][mu])/den1[alpha-1];

	V_ave = V_ave+ e[2][mu]*(distr_mom_k[alpha][mu]*exp_Nt_t - ((double)Include_k0_noise)*distr_mom_0[alpha][mu])/den2[alpha-1];


      }

     
    }
  }


 
  
  
  distr_t_list DIFF2 = 2.0*( distr_mom_k[2][1]*exp_Nt_t -((double)Include_k0_noise)*distr_mom_0[2][1] - distr_mom_k[1][2]*exp_Nt_t + ((double)Include_k0_noise)*distr_mom_0[1][2]);
  distr_t_list DIFF= DIFF2/(Mom.k()[2]*Meson_mass);
  //for(int t=0; t < V_ave.size(); t++) cout<<V_ave.distr_list[t].ave()<<" "<<V_ave.distr_list[t].err()<<" "<<DIFF.distr_list[t].ave()<<" "<<DIFF.distr_list[t].err()<<endl;
  
  if(is_mom_zero) Mom.Set_theta(theta_t_old, 2);
  
  return V_ave;
  //return DIFF2;


}

distr_t_list A_ave_unpolarized(vector<vector<distr_t_list>>& distr_mom_k) {

  VVfloat e(4);
  e[0] = {0.0, 0.0, 0.0, 0.0};
  e[3] = {0.0, 0.0, 0.0, 0.0};
  e[1] = {0.0, -1.0/sqrt(2), -1.0/sqrt(2), 0.0};
  e[2] = {0.0, 1.0/sqrt(2), -1.0/sqrt(2), 0.0};



  if(distr_mom_k.size() != 4) crash("A_ave_unpolarized called with a v<v<distr_t_list>> with size != alphas= 4");
  for(auto & distr_mom_k_alpha : distr_mom_k) if(distr_mom_k_alpha.size() != 4) crash("A_ave unpolarized called with a with a vv<distr_t_list>> having an element of size != 4");
  int size = distr_mom_k[0][0].size();
  if(size==0) crash(" A_ave_unpolarized called with Nt=0");
  int distr_size = distr_mom_k[0][0].distr_list[0].size();
  if(distr_size ==0) crash(" A_ave unpolarized called with a vv<distr_t_list>> having zero statistics");
  distr_t_list A_ave(UseJack, size, distr_size);

 
  for(int alpha=1; alpha <= 2; alpha ++) {
 
    if(distr_mom_k[alpha].size() != 4) crash(" A_ave_unpolarized called with a v<v<distr_t_list>> having an element with size != 4");
 
    for(int mu=0; mu< 4; mu ++) {

 
      A_ave = A_ave+ distr_mom_k[alpha][mu]*(e[1][mu]/e[1][alpha]);

      A_ave = A_ave+ distr_mom_k[alpha][mu]*(e[2][mu]/e[2][alpha]);
    }
  }


  /*
    cout<<"A##########"<<endl;
    distr_t_list DIFF = 2.0*( distr_mom_k[1][1] +distr_mom_k[2][2]);
    for(int t=0; t < A_ave.size(); t++) cout<<A_ave.distr_list[t].ave()<<" "<<A_ave.distr_list[t].err()<<" "<<DIFF.distr_list[t].ave()<<" "<<DIFF.distr_list[t].err()<<endl;
    cout<<"############"<<endl;
  */

  //return distr_mom_k[1][1];
  return A_ave;
  
}


distr_t_list H_2(const distr_t_list& H30,const distr_t_list& H11,const distr_t_list& H03,const distr_t_list& H30_0,const distr_t_list& H11_0,const distr_t_list& H03_0,pt3_momenta& Mom,const distr_t& m) {//valid for p=0, k= kz

   auto Power = [&](const distr_t& A, int n) -> distr_t{
    distr_t ret_val = A;
    for(int i=1;i<=n;i++) ret_val= ret_val*A;
    return ret_val;
  };

  double Eg= Mom.Egamma();
  double kz= Mom.k()[2];
  double ksq = pow(Mom.virt(),2);
  
  distr_t_list H30_sub = H30 - H11_0*(m-Eg)/( (2.0*m*Eg/kz) -(ksq/kz));
  distr_t_list H11_sub = H11 - H11_0 ;
  distr_t_list H03_sub = H03 - H11_0*(2.0*m -Eg)/(2.0*m*(Eg/kz) -(ksq/kz));

  distr_t_list H_2 = -1.0*H03_sub*Eg/( (pow(Eg,2)-ksq)*kz*Power(m,2));
  H_2 = H_2+ H30_sub*(1.0/(Eg*kz*Power(m,2)));
  H_2 = H_2+ H11_sub*(1.0/( (pow(Eg,2)-ksq)*Power(m,2)));
  return H_2;

}


distr_t_list H_1(const distr_t_list& H30,const distr_t_list& H11,const distr_t_list& H03,const distr_t_list& H30_0,const distr_t_list& H11_0,const distr_t_list& H03_0,pt3_momenta& Mom,const distr_t& m) //valid for p=0, k=kz
{

  auto Power = [&](const distr_t& A, int n) -> distr_t{
    distr_t ret_val = A;
    for(int i=1;i<=n;i++) ret_val= ret_val*A;
    return ret_val;
  };
  
  double Eg= Mom.Egamma();
  double kz= Mom.k()[2];
  double ksq = pow(Mom.virt(),2);
  
  distr_t_list H30_sub = H30 - H11_0*(m-Eg)/( (2.0*m*Eg/kz) -(ksq/kz));
  distr_t_list H11_sub = H11 - H11_0 ;
  distr_t_list H03_sub = H03 - H11_0*(2.0*m -Eg)/(2.0*m*(Eg/kz) -(ksq/kz));


  distr_t_list H_1 = H03_sub*( (Eg-m)*(Eg*m -ksq))/(kz*Power(m,2)*(pow(Eg,2)-ksq));
  H_1 = H_1+ H30_sub*(ksq-Eg*m)/(Eg*kz*Power(m,2));
  H_1 = H_1+ H11_sub*(ksq + m*(m-2.0*Eg))/(Power(m,2)*(pow(Eg,2) -ksq));

  
  return H_1;



}


distr_t_list FA_off(const distr_t_list& H30,const distr_t_list& H11,const distr_t_list& H03,const distr_t_list& H30_0,const distr_t_list& H11_0,const distr_t_list& H03_0,pt3_momenta& Mom,const distr_t &m) //valid for p=0, k=kz
{

  double Eg= Mom.Egamma();
  double kz= Mom.k()[2]; 
  double ksq = pow(Mom.virt(),2);
  
  distr_t_list H30_sub = H30 - H11_0*(m-Eg)/( (2.0*m*Eg/kz) -(ksq/kz));
  distr_t_list H11_sub = H11 - H11_0 ;
  distr_t_list H03_sub = H03 - H11_0*(2.0*m -Eg)/(2.0*m*(Eg/kz) -(ksq/kz));

  distr_t_list FA_off = H03_sub*(ksq*(m-Eg))/(kz*m*(pow(Eg,2)-ksq));
  FA_off = FA_off+ H30_sub*ksq/(Eg*kz*m);
  FA_off = FA_off+ H11_sub*(ksq-Eg*m)/(m*(pow(Eg,2)-ksq));
 
 
  return FA_off;
}







void Compute_form_factors() {




  
  string dir= "../3pt_data/"+Meson;
  int Nens=0;
  vector<string> Ens_tag;
  

  //read number of ensembles and store the paths

  boost::filesystem::directory_iterator end_itr;

  for (boost::filesystem::directory_iterator itr(dir); itr != end_itr; ++itr) {
  

    if(!boost::filesystem::is_regular_file(itr->path())) {

      Nens++;
      Ens_tag.push_back(itr->path().string().substr(dir.length()+1));
    }
  }



  
  pt3_momenta_list mom3_l;
  pt2_momenta_list mom2_l;


  
  
  
  
  

  for (int i_ens=0; i_ens<Nens;i_ens++) {


    FILE *stream_3pt;
    FILE *stream_2pt;

    struct header_virph header_3pt;
    struct header_virph header_2pt;

   

    string Path3pt = dir+"/"+Ens_tag[i_ens]+"/conf.virtualph.dat";
    string Path2pt = dir+"/"+Ens_tag[i_ens]+"/conf.virtualph.dat2";



    stream_3pt= fopen(Path3pt.c_str(), "rb");
    stream_2pt= fopen(Path2pt.c_str(), "rb");



    read_header_bin(stream_3pt, header_3pt);
    read_header_bin(stream_2pt, header_2pt);

    //Forward lattice  informations to 3pt_momenta_list
    double ml, l, Nt;
    Read_pars_from_ensemble_tag(Ens_tag[i_ens], ml, l , Nt);
    Add_to_mom_list(mom3_l, header_3pt, l);
    Add_to_mom_list(mom2_l, header_2pt, l);
    mom2_l.Nt.push_back(Nt);
    mom3_l.Nt.push_back(Nt);


    //print size of header file
    cout<<"Ens: "<<Ens_tag[i_ens]<<", header size 2pt (Bytes): "<<header_2pt.header_size<<endl;
    cout<<"Number of configs (from 2pt): "<<Get_number_of_configs_2pt(stream_2pt, header_2pt)<<endl;
    cout<<"Ens: "<<Ens_tag[i_ens]<<", header size 3pt (Bytes): "<<header_3pt.header_size<<endl;
    cout<<"Number of configs (from 3pt): "<<Get_number_of_configs_3pt(stream_3pt, header_3pt)<<endl;
    cout<<"Number of smearing levels present: "<<header_3pt.nqsml<<endl;


    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/form_factors");
    boost::filesystem::create_directory("../data/form_factors/"+Meson);
    boost::filesystem::create_directory("../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]);
    boost::filesystem::create_directory("../data/form_factors/"+Meson+"/C_"+Ens_tag[i_ens]);
    boost::filesystem::create_directory("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST");
    boost::filesystem::create_directory("../data/form_factors/"+Meson+"/Contaminations");
    
    //PRINT HEADERS
    //#########################################################
    ofstream P_head;
    if(Determine_contaminations) {
    P_head.open("../data/form_factors/"+Meson+"/Contaminations/contaminations_RV_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"a1"<<setw(20)<<"a1_err"<<setw(20)<<"a2"<<setw(20)<<"a2_err"<<setw(20)<<"a3"<<setw(20)<<"a3_err"<<setw(20)<<"a5"<<setw(20)<<"a5_err"<<setw(20)<<"ch2/ndof"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/Contaminations/contaminations_RA_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"a1"<<setw(20)<<"a1_err"<<setw(20)<<"a2"<<setw(20)<<"a2_err"<<setw(20)<<"a3"<<setw(20)<<"a3_err"<<setw(20)<<"a5"<<setw(20)<<"a5_err"<<setw(20)<<"ch2/ndof"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/Contaminations/contaminations_FA_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"a1"<<setw(20)<<"a1_err"<<setw(20)<<"a2"<<setw(20)<<"a2_err"<<setw(20)<<"a3"<<setw(20)<<"a3_err"<<setw(20)<<"a5"<<setw(20)<<"a5_err"<<setw(20)<<"ch2/ndof"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/Contaminations/contaminations_FV_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"a1"<<setw(20)<<"a1_err"<<setw(20)<<"a2"<<setw(20)<<"a2_err"<<setw(20)<<"a3"<<setw(20)<<"a3_err"<<setw(20)<<"a5"<<setw(20)<<"a5_err"<<setw(20)<<"ch2/ndof"<<endl;
    P_head.close();
    }
    P_head.open("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/vector_form_factors_list_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"F"<<setw(20)<<"F_err"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/axial_form_factors_list_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"F"<<setw(20)<<"F_err"<<endl;
    P_head.close();
    P_head.open("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/virtual_form_factors_list_"+Ens_tag[i_ens]+".dat");
    P_head<<"#xg"<<setw(20)<<"offsh"<<setw(20)<<"sm_lev"<<setw(20)<<"H1"<<setw(20)<<"H1_err"<<setw(20)<<"H2"<<setw(20)<<"H2_err"<<setw(20)<<"FA_off"<<setw(20)<<"FA_off_err"<<setw(20)<<"FV_off"<<setw(20)<<"FV_off_err"<<endl;
    P_head.close();
    //##########################################################


    
    for(int nqsm=0; nqsm<header_3pt.nqsml;nqsm++) {

      int Analyze_sm_lev=1;
      if(header_3pt.nqsml == 1) sm_lev=0;
      else {
	cout<<"Should I analyze smearing level: "<<nqsm<<" ?"<<endl;
	cout<<"1->yes, 0->no"<<endl;
	cin>>Analyze_sm_lev;
      }

      
      if(Analyze_sm_lev) {
	sm_lev=nqsm;


	//STORAGE FOR 2-pt FUNC OBS
	vector<distr_t_list> fp_distr_list, m_distr_list;
	vector<distr_t> fp_fit_distr, m_fit_distr, sqrt_overlap_fit_distr;

	



	//init CorrAnalysis class
	//CorrAnalysis corr(UseJack, Get_number_of_configs_3pt(stream_3pt, header_3pt), nboots);
	CorrAnalysis corr(UseJack, Njacks, nboots); //in case you want to use a different number of NJacks
	corr.Nt=Nt;


	//DEFINE LAMBDA TO REMOVE EXPONENTIAL DEPENDENCE ON MESON ENERGY FROM THE CORRELATORS
	//###############################################################################à
	auto Exp_lat = [&] (double a, double b, double c) { return (b<c/2)?1.0/(exp(-a*b)):1.0/(exp(-a*(c-b)));};
	//#############################################################################

    

	//2pts
	Get_Tmin_Tmax("2pt", Ens_tag[i_ens], corr, 0, "V");
	//#########################################################
	

	int ncomb2pt= header_2pt.ncomb;
	int ncomb3pt= header_3pt.ncomb;

	cout<<"#########################"<<endl;
	cout<<"number of combs in 2pt: "<<ncomb2pt<<endl;
	cout<<"number of combs in 3pt: "<<ncomb3pt<<endl;
	cout<<"#########################"<<endl;

      

	//THIS PART OF THE CODE HANDLES THE TWO-PT FUNCTION. IT COMPUTES EFFECTIVE MASSES AND MESON'S DECAY CONSTANT, FOR ALL KINEMATICS.
	//#############################################################################
	//#############################################################################
	//#############################################################################

	for(int icomb2pt=0; icomb2pt<ncomb2pt;icomb2pt++) {

	  auto c = header_2pt.comb[icomb2pt];
	  string tag=to_string(icomb2pt);
	  distr_t_list overlap_distr_list = corr.residue_t(Get_obs_2pt(stream_2pt,header_2pt,0,icomb2pt,0,0), "");
	  distr_t_list overlap_smeared_distr_list = corr.residue_t(Get_obs_2pt(stream_2pt,header_2pt,0,icomb2pt,0,sm_lev), "");
	  distr_t overlap_fit_distr =corr.Fit_distr(overlap_distr_list);
	  distr_t overlap_smeared_fit_distr=corr.Fit_distr( overlap_smeared_distr_list);
	  fp_distr_list.push_back( (c.mu1+c.mu2)*corr.decay_constant_t( Get_obs_2pt(stream_2pt, header_2pt, 0, icomb2pt, 0, 0), "../data/form_factors/"+Meson+"/fp_"+Ens_tag[i_ens]+"_sm_lev"+to_string(sm_lev)+"_"+tag));
	  m_distr_list.push_back( corr.effective_mass_t(Get_obs_2pt(stream_2pt, header_2pt, 0, icomb2pt, 0, sm_lev), "../data/form_factors/"+Meson+"/meson_mass_"+Ens_tag[i_ens]+"_sm_lev"+to_string(sm_lev)+"_"+tag));
	  fp_fit_distr.push_back(corr.Fit_distr(fp_distr_list[icomb2pt]));
	  m_fit_distr.push_back(corr.Fit_distr(m_distr_list[icomb2pt]));
	  auto sq= [](double x)->double {return sqrt(x);};
	  sqrt_overlap_fit_distr.push_back( overlap_smeared_fit_distr/(distr_t::f_of_distr(sq,overlap_fit_distr)));
	}

	//#############################################################################
	//#############################################################################
	//#############################################################################
	//end anaysis 2pts


    
	//BEGINNING ANALYSIS OF 3PT-FUNCTIONS
	//#############################################################################
	//#############################################################################
	//#############################################################################
      
	//loop over kinematics
	for(int icomb3pt=0; icomb3pt<ncomb3pt;icomb3pt++) {

	  auto c = header_3pt.comb[icomb3pt];
	  
	  
	  vector<vector<distr_t_list>> distr_V(4), distr_V_k0(4), distr_A(4), distr_A_k0(4);
	  vector<vector<distr_t_list>> distr_no_symm_V(4), distr_no_symm_V_k0(4), distr_no_symm_A(4), distr_no_symm_A_k0(4);
	  for(int i=0;i<4;i++) { distr_V[i].resize(4); distr_V_k0[i].resize(4); distr_A[i].resize(4); distr_A_k0[i].resize(4);}
	  for(int i=0;i<4;i++) { distr_no_symm_V[i].resize(4); distr_no_symm_V_k0[i].resize(4); distr_no_symm_A[i].resize(4); distr_no_symm_A_k0[i].resize(4);}
	  
	 

	  //infos about the kinematic
	  //####################################################

      
	  int icomb_k0= Get_comb_k0(header_3pt, icomb3pt);
	  int symmetric_comb=Get_symmetric_comb(header_3pt, icomb3pt);
	  int symmetric_comb_k0 = Get_symmetric_comb(header_3pt, icomb_k0);
	  int pt2_k0p0 = Get_2pt_k0p0(header_2pt, c.mu1, c.mu2);
	  int pt2_p = Get_2pt_p(header_2pt, c.i0, c.is);
	  double Egamma = mom3_l.mom[i_ens][icomb3pt].Egamma();
	  double xg = mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[pt2_k0p0]).ave();
	  double offsh= mom3_l.mom[i_ens][icomb3pt].virt();
	  double Egamma_T= sinh(Egamma)*(1.0-exp(-Egamma*Nt));
	  
	 
	  	 
     
     
	  //PRINT INFOS
	  cout<<"###BEG###"<<endl;
	  cout<<"icomb: "<<icomb3pt<<" th_t: "<<c.tht[2]<<endl;
	  cout<<"nconfs: "<<Get_number_of_configs_3pt(stream_3pt, header_3pt)<<endl;
	  cout<<"icomb_k0: "<<icomb_k0<<" th_t: "<<header_3pt.comb[icomb_k0].tht[2]<<endl;
	  cout<<"icomb_symm: "<<symmetric_comb<<" th_t: "<<header_3pt.comb[symmetric_comb].tht[2]<<endl;
	  cout<<"icomb_symm_k0: "<<symmetric_comb_k0<<" th_t: "<<header_3pt.comb[symmetric_comb_k0].tht[2]<<endl;
	  cout<<"fp: "<<fp_fit_distr[pt2_k0p0].ave()<<"("<<fp_fit_distr[pt2_k0p0].err()<<")"<<endl;
	  cout<<"k_z: "<<mom3_l.mom[i_ens][icomb3pt].k()[2]<<endl;
	  cout<<"Egamma["<<icomb3pt<<"] :"<<Egamma<<endl;
	  cout<<"EgammaT["<<icomb3pt<<"] :"<<Egamma_T<<endl;
	  cout<<"Mass: "<<m_fit_distr[pt2_k0p0].ave()<<"("<<m_fit_distr[pt2_k0p0].err()<<")"<<endl;
	  cout<<"x_gamma: "<<xg<<endl;
	  cout<<"offshellness: "<<offsh<<endl;
	  cout<<"##############"<<endl;
	


	  //####################################################

	  //TO REMOVE EXPONENTIAL IN THE PHOTON ENERGY
	  //####################################################
	  Vfloat exp_Nt_t;
	  for(int t=0; t<Nt;t++) exp_Nt_t.push_back( exp( fabs((Nt/2-t))*Egamma));
	  
	  //####################################################


	  //loop over alpha and mu
        
	  for(int alpha=0; alpha<4;alpha++) {
	    for(int mu=0; mu<4;mu++) {

	  
	      // compute_vector_correlators C_V_alpha^mu for each alpha and mu
	      //######################################################################

	      corr.Reflection_sign = -1;
	      corr.Perform_Nt_t_average=1;
   
	      distr_no_symm_V[alpha][mu] = e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, icomb3pt, alpha, mu, "V", sm_lev), "");
	      distr_no_symm_V_k0[alpha][mu] =e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1,icomb_k0, alpha, mu, "V", sm_lev), "");
	      distr_V[alpha][mu] = distr_no_symm_V[alpha][mu] + e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, symmetric_comb, alpha, mu, "V", sm_lev), "");
	      distr_V_k0[alpha][mu] = distr_no_symm_V_k0[alpha][mu]  + e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, 1, symmetric_comb_k0, alpha, mu, "V", sm_lev), "");

	      //######################################################################


	      //compute_axial_correlators C_A_alpha^\mu for each alpha and mu
	      //######################################################################
	      corr.Reflection_sign = 1;
	      int Im_Re= 0; double parity=1.0; double sign=1;
	      if( (alpha==0 || mu==0) && (alpha != 0 || mu != 0)  ) {Im_Re=1; corr.Reflection_sign=-1;sign=1; parity=-1;}
	      distr_no_symm_A[alpha][mu] = parity*e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, Im_Re, icomb3pt, alpha, mu, "A", sm_lev), "");
	      distr_no_symm_A_k0[alpha][mu] =parity*e_f1*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt,Im_Re, icomb_k0, alpha, mu, "A", sm_lev), "");
	      distr_A[alpha][mu] = distr_no_symm_A[alpha][mu]  - parity*sign*e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, Im_Re, symmetric_comb, alpha, mu, "A", sm_lev), "");
	      distr_A_k0[alpha][mu] = distr_no_symm_A_k0[alpha][mu]   - parity*sign*e_f2*corr.corr_t(Get_obs_3pt(stream_3pt, header_3pt, Im_Re, symmetric_comb_k0, alpha, mu, "A", sm_lev), "");

	      //######################################################################
	      corr.Reflection_sign= 1;
	      sign=1;
	      parity=1.0;
	      

	      //PRINT TO FILE THE TENSOR C_W_alpha^\mu
	      Print_To_File({}, {distr_no_symm_A[alpha][mu].ave(), distr_no_symm_A[alpha][mu].err(), distr_no_symm_V[alpha][mu].ave(), distr_no_symm_V[alpha][mu].err()}, "../data/form_factors/"+Meson+"/C_"+Ens_tag[i_ens]+"/no_symm_icomb"+to_string(icomb3pt)+"_alpha"+to_string(alpha)+"_mu"+to_string(mu)+"_sm_lev_"+to_string(sm_lev)+".dat","", "#");

	      Print_To_File({}, {distr_A[alpha][mu].ave(), distr_A[alpha][mu].err(), distr_V[alpha][mu].ave(), distr_V[alpha][mu].err()}, "../data/form_factors/"+Meson+"/C_"+Ens_tag[i_ens]+"/icomb_"+to_string(icomb3pt)+"_alpha"+to_string(alpha)+"_mu"+to_string(mu)+"_sm_lev_"+to_string(sm_lev)+".dat","", "#");




	      //COMPUTE THE TENSOR H_W_alpha^\mu
	      //########################################################################
	 
	      distr_t_list distr_no_symm_A_exp= (m_fit_distr[pt2_p]/sqrt_overlap_fit_distr[pt2_p])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p], Nt)*(distr_no_symm_A[alpha][mu]*exp_Nt_t); 
	      distr_t_list distr_no_symm_V_exp=  (m_fit_distr[pt2_p]/sqrt_overlap_fit_distr[pt2_p])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p], Nt)*(distr_no_symm_V[alpha][mu]*exp_Nt_t);
	      distr_t_list distr_A_exp= (m_fit_distr[pt2_p]/sqrt_overlap_fit_distr[pt2_p])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p], Nt)*(distr_A[alpha][mu]*exp_Nt_t); 
	      distr_t_list distr_V_exp=  (m_fit_distr[pt2_p]/sqrt_overlap_fit_distr[pt2_p])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p], Nt)*(distr_V[alpha][mu]*exp_Nt_t);

	      //PRINT TO FILE THE TENSOR H_alpha^mu
	      Print_To_File({}, {distr_no_symm_A_exp.ave(), distr_no_symm_A_exp.err(), distr_no_symm_V_exp.ave(), distr_no_symm_V_exp.err()}, "../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]+"/no_symm_icomb"+to_string(icomb3pt)+"_alpha"+to_string(alpha)+"_mu"+to_string(mu)+"_sm_lev_"+to_string(sm_lev)+".dat","", "#");
	      Print_To_File({}, {distr_A_exp.ave(), distr_A_exp.err(), distr_V_exp.ave(), distr_V_exp.err()}, "../data/form_factors/"+Meson+"/H_"+Ens_tag[i_ens]+"/icomb_"+to_string(icomb3pt)+"_alpha"+to_string(alpha)+"_mu"+to_string(mu)+"_sm_lev_"+to_string(sm_lev)+".dat","", "#");


	      //########################################################################

	 
	    }
	  }

	 
           
	       
	  //COMPUTE R and bar_R and form factors
	  distr_t_list bar_R_A_distr, bar_R_V_distr, R_A_distr, R_V_distr;

	  if(icomb3pt < ncomb3pt/2.0) {
      
	  //####################################################################################################################################
	    if(offsh == 0) {
	      distr_t_list V_ave = V_ave_unpolarized(distr_V,distr_V_k0, m_fit_distr[pt2_p], mom3_l.mom[i_ens][icomb3pt],Nt/2);
	      distr_t_list A_ave = A_ave_unpolarized(distr_A);
	      distr_t_list A_ave_zero_mom = A_ave_unpolarized(distr_A_k0);  
	      distr_t F_A_distr_factor= (2.0*fp_fit_distr[pt2_k0p0]/(m_fit_distr[pt2_k0p0]*mom3_l.mom[i_ens][icomb3pt].x_gamma(m_fit_distr[pt2_k0p0])));
	      auto Exp_photon= [&](double a, int t) { return exp(fabs(corr.Nt/2-t)*a);};

	      distr_t_list eff_photon_energy= corr.effective_mass_t( A_ave_zero_mom/A_ave,"");
	      
	      //bar_R_A_distr=(distr_t_list::f_of_distr_list(Exp_photon, eff_photon_energy)*A_ave/A_ave_zero_mom);
	      bar_R_A_distr=(exp_Nt_t*A_ave/A_ave_zero_mom);
	      cout<<"Eff photon energy: "<<corr.Fit_distr(eff_photon_energy).ave()<<" "<<corr.Fit_distr(eff_photon_energy).err()<<endl;
	      bar_R_V_distr= (V_ave/A_ave_zero_mom)*m_fit_distr[pt2_k0p0]*fp_fit_distr[pt2_k0p0];
	      R_A_distr= (-1.0*m_fit_distr[pt2_p]/(2.0*m_fit_distr[pt2_k0p0]))*(A_ave/sqrt_overlap_fit_distr[icomb3pt])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p], Nt)*exp_Nt_t;
	      R_V_distr= (m_fit_distr[pt2_p]*m_fit_distr[pt2_k0p0]/4.0)*(V_ave/sqrt_overlap_fit_distr[icomb3pt])*2.0*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p], Nt);
	      distr_t_list F_A_distr = (bar_R_A_distr-1.0)*F_A_distr_factor;
	      
	      //#####################################################################################################################################
	      
	      //FIT form factors and R estimators
	      // #########################################
	      Get_Tmin_Tmax("3pt", Ens_tag[i_ens], corr,xg , "V");
	      int T_min_V= corr.Tmin;
	      int T_max_V= corr.Tmax;
	      distr_t fit_result_V = corr.Fit_distr(bar_R_V_distr);
	      distr_t fit_result_R_V = corr.Fit_distr(R_V_distr);
	      Get_Tmin_Tmax("3pt", Ens_tag[i_ens], corr,xg , "A");
	      int T_min_A= corr.Tmin;
	      int T_max_A=corr.Tmax;
	      distr_t fit_result_A = corr.Fit_distr(F_A_distr);
	      distr_t fit_result_R_A = corr.Fit_distr(R_A_distr);
	      // #####################################
      

      
	      //PRINT R-estimators and form factors to file
	      //###############################################################################
	      VVfloat To_print_bar_R({bar_R_A_distr.ave(), bar_R_A_distr.err(), bar_R_V_distr.ave(), bar_R_V_distr.err(), R_A_distr.ave(), R_A_distr.err(), R_V_distr.ave(), R_V_distr.err()});
	       
	     
	      Print_To_File({}, To_print_bar_R, "../data/form_factors/"+Meson+"/R_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+"_smlev_"+to_string(sm_lev)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "#t       bar_RA    bar_RA_err     bar_RV      bar_RV_err    RA      RA_err        RV      RV_err");
	      Print_To_File({}, {F_A_distr.ave(), F_A_distr.err(), bar_R_V_distr.ave(), bar_R_V_distr.err() } , "../data/form_factors/"+Meson+"/F_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+"_smlev_"+to_string(sm_lev)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "#t       FA    FA_err     FV      FV_err");
	      if(xg > 0) { //print form factors only for xg > 0
		ofstream Print_Fitted_form_factors("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/axial_form_factors_list_"+Ens_tag[i_ens]+".dat",ofstream::app);
		Print_Fitted_form_factors<<xg<<setw(20)<<offsh<<setw(20)<<sm_lev<<setw(20)<<fit_result_A.ave()<<setw(20)<<fit_result_A.err()<<endl;
		Print_Fitted_form_factors.close();
		Print_Fitted_form_factors.open("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/vector_form_factors_list_"+Ens_tag[i_ens]+".dat",ofstream::app);
		Print_Fitted_form_factors<<xg<<setw(20)<<offsh<<setw(20)<<sm_lev<<setw(20)<<fit_result_V.ave()<<setw(20)<<fit_result_V.err()<<endl;
		Print_Fitted_form_factors.close();
	      }
	      //###############################################################################
      
	      
	      //plot form factors
	      //##################################################################################################
	      if(xg > 1.0e-8) { //plot only if xg > 0
		Plot_form_factors("V", bar_R_V_distr, fit_result_V, T_min_V, T_max_V, Nt, Ens_tag[i_ens], xg, offsh, sm_lev);
		Plot_form_factors("A", F_A_distr, fit_result_A, T_min_A, T_max_A, Nt, Ens_tag[i_ens], xg, offsh, sm_lev);
	      }
	      //##################################################################################################
      

	      //FIT EXPONENTIAL CONTAMINATIONS FROM R and bar_R 
	      if(Determine_contaminations && xg > 1.0e-8) { //fit contaminations only if xg > 0
	   
		Fit_contaminations("FA", F_A_distr, fit_result_A, Ens_tag[i_ens], Nt, xg, offsh, sm_lev);
		Fit_contaminations("FV", bar_R_V_distr, fit_result_V, Ens_tag[i_ens], Nt, xg, offsh, sm_lev);
		Fit_contaminations("RA", R_A_distr, fit_result_R_A, Ens_tag[i_ens], Nt, xg, offsh, sm_lev);
		Fit_contaminations("RV", R_V_distr, fit_result_R_V, Ens_tag[i_ens],  Nt, xg, offsh, sm_lev); 
	      }
	    }

	    //determine form_factors H1 H2 e FA_offsh (only valid if meson at rest and k=kx)
	    if( xg > 1.0e-8 && mom3_l.mom[i_ens][icomb3pt].k()[2] !=  0) { //only if xg, kz != 0
	      //distr_t_list H_resc= (m_fit_distr[pt2_p]/sqrt_overlap_fit_distr[pt2_p])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p], Nt)*exp_Nt_t;
	      //distr_t_list H_resc_0 = (m_fit_distr[pt2_k0p0]/sqrt_overlap_fit_distr[pt2_k0p0])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_k0p0], Nt);
	      distr_t_list H1 = H_1(distr_A[0][3]*exp_Nt_t, 0.5*(distr_A[1][1]+distr_A[2][2])*exp_Nt_t, distr_A[3][0]*exp_Nt_t, distr_A_k0[0][3], 0.5*(distr_A_k0[1][1]+distr_A_k0[2][2]), distr_A_k0[3][0], mom3_l.mom[i_ens][icomb3pt], m_fit_distr[pt2_k0p0])*(-2.0*fp_fit_distr[pt2_k0p0]/(distr_A_k0[1][1]+distr_A_k0[2][2]));
	      distr_t_list H2 = H_2(distr_A[0][3]*exp_Nt_t, 0.5*(distr_A[1][1]+distr_A[2][2])*exp_Nt_t, distr_A[3][0]*exp_Nt_t, distr_A_k0[0][3], 0.5*(distr_A_k0[1][1]+distr_A_k0[2][2]), distr_A_k0[3][0], mom3_l.mom[i_ens][icomb3pt], m_fit_distr[pt2_k0p0])*(-2.0*fp_fit_distr[pt2_k0p0]/(distr_A_k0[1][1]+distr_A_k0[2][2]));
	      distr_t_list FAoff = FA_off(distr_A[0][3]*exp_Nt_t, 0.5*(distr_A[1][1]+distr_A[2][2])*exp_Nt_t, distr_A[3][0]*exp_Nt_t, distr_A_k0[0][3], 0.5*(distr_A_k0[1][1]+distr_A_k0[2][2]), distr_A_k0[3][0], mom3_l.mom[i_ens][icomb3pt], m_fit_distr[pt2_k0p0])*(-2.0*fp_fit_distr[pt2_k0p0]/(distr_A_k0[1][1]+distr_A_k0[2][2]));
	      //print to file
	      distr_t_list FVoff = (V_ave_unpolarized(distr_V,distr_V_k0, m_fit_distr[pt2_p], mom3_l.mom[i_ens][icomb3pt],Nt/2)/A_ave_unpolarized(distr_A_k0))*m_fit_distr[pt2_k0p0]*fp_fit_distr[pt2_k0p0];
	      Print_To_File({}, {H1.ave(), H1.err(), H2.ave(), H2.err(), FAoff.ave(), FAoff.err(), FVoff.ave(), FVoff.err()}, "../data/form_factors/"+Meson+"/FH_off_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+"_smlev_"+to_string(sm_lev)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "#t       H1    H1_err     H2      H2_err    FA_off      FA_off_err     FV_off     FV_off_err ");
	      //fit in time interval
	      corr.Tmin=12;
	      corr.Tmax=18;
	      distr_t H1_fit = corr.Fit_distr(H1);
	      corr.Tmin=10;
	      corr.Tmax=15;
	      distr_t H2_fit = corr.Fit_distr(H2);
	      corr.Tmax=15;
	      distr_t FAoff_fit = corr.Fit_distr(FAoff);
	      corr.Tmin=8;
	      corr.Tmax=14;
	      distr_t FVoff_fit = corr.Fit_distr(FVoff);
	      //save in file
	       if(xg > 1.0e-8) { //print form factors only for xg > 0
		ofstream Print_virtual_form_factors("../data/form_factors/"+Meson+"/FORM_FACTOR_LIST/virtual_form_factors_list_"+Ens_tag[i_ens]+".dat",ofstream::app);
		Print_virtual_form_factors<<xg<<setw(20)<<offsh<<setw(20)<<sm_lev<<setw(20)<<H1_fit.ave()<<setw(20)<<H1_fit.err()<<setw(20)<<H2_fit.ave()<<setw(20)<<H2_fit.err()<<setw(20)<<FAoff_fit.ave()<<setw(20)<<FAoff_fit.err()<<setw(20)<<FVoff_fit.ave()<<setw(20)<<FVoff_fit.err()<<endl;
		Print_virtual_form_factors.close();


		//print subtracted hadronic tensor
		auto Power = [&](const distr_t& A, int n) -> distr_t{
		  distr_t ret_val = A;
		  for(int i=1;i<=n;i++) ret_val= ret_val*A;
		  return ret_val;
		};
		
		double Eg= mom3_l.mom[i_ens][icomb3pt].Egamma();
		double kz= mom3_l.mom[i_ens][icomb3pt].k()[2];
		double ksq = pow(mom3_l.mom[i_ens][icomb3pt].virt(),2);
		distr_t m=m_fit_distr[pt2_k0p0];
		distr_t_list H_resc_0 = (m_fit_distr[pt2_k0p0]/sqrt_overlap_fit_distr[pt2_k0p0])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_k0p0], Nt);
		distr_t_list H_resc= (m_fit_distr[pt2_p]/sqrt_overlap_fit_distr[pt2_p])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_p], Nt)*exp_Nt_t;

		distr_t_list H30_sub = distr_A[0][3]*H_resc - distr_A_k0[1][1]*H_resc_0*(m-Eg)/( (2.0*m*Eg/kz) -(ksq/kz));
		distr_t_list H11_sub = distr_A[1][1]*H_resc - distr_A_k0[1][1]*H_resc_0;
		distr_t_list H33_sub = distr_A[3][3]*H_resc -distr_A_k0[3][3]*H_resc_0*(2.0*m*Eg -ksq -pow(kz,2))/(2.0*m*Eg - ksq);
		distr_t_list H03_sub = distr_A[3][0]*H_resc - distr_A_k0[1][1]*H_resc_0*(2.0*m -Eg)/(2.0*m*(Eg/kz) -(ksq/kz));
		distr_t_list H00_sub = distr_A[0][0]*H_resc - distr_A_k0[0][0]*H_resc_0*(2.0*m -Eg + pow(Eg,2)/m -ksq/m)/(2.0*m -ksq/Eg);

		Print_To_File({}, {H30_sub.ave(), H30_sub.err(), H11_sub.ave(), H11_sub.err(), H33_sub.ave(), H33_sub.err(), H03_sub.ave(), H03_sub.err(), H00_sub.ave(), H00_sub.err()},  "../data/form_factors/"+Meson+"/Hadr_tens_subtracted_"+Ens_tag[i_ens]+"_"+mom3_l.mom[i_ens][icomb3pt].name()+"_smlev_"+to_string(sm_lev)+"_k0_noise_"+to_string(Include_k0_noise)+".dat", "", "#t     H30         H30_err        H11        H11_err     H33       H33_err       H03         H03_err            H00               H00_err");
        
		
	       }

	      
	    
	    }

	    if(icomb3pt==0) {
	      //print fp determined from various matrix elements of the hadronic tensor
	      distr_t_list H_resc_0 = (m_fit_distr[pt2_k0p0]/sqrt_overlap_fit_distr[pt2_k0p0])*distr_t_list::f_of_distr(Exp_lat, m_fit_distr[pt2_k0p0], Nt);
	      distr_t_list fp1= -1.0*distr_A[3][0]*H_resc_0;
	      distr_t_list fp2= -2.0*distr_A[0][3]*H_resc_0;
	      distr_t_list fp3= -1.0*distr_A[1][1]*H_resc_0;
	      distr_t_list fp4= -1.0*distr_A[3][3]*H_resc_0;

	     
	      Print_To_File({}, {fp1.ave(), fp1.err(), fp2.ave(),fp2.err(),fp3.ave(),fp3.err(),fp4.ave(),fp4.err()},"../data/form_factors/"+Meson+"/fp_3pt_"+Ens_tag[i_ens]+"_sm_"+to_string(sm_lev)+".dat", "", "#t     fp1       fp1_err     fp2       fp2_err         fp3        fp3_err            fp4              fp4_err");
	      
	      
	    }
	   
	  }
	}
      }
    }
    fclose(stream_3pt);
    fclose(stream_2pt);
  }
  
  return;
}

