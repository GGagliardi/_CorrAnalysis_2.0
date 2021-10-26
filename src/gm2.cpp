#include "../include/gm2.h"


using namespace std;

//set constants
namespace plt = matplotlibcpp;

constexpr double kappa=2.837297;
const double MPiPhys=0.135;
const double alpha = 1.0/137.04;
const double e2 = alpha*4.0*M_PI;
const double r0 = pow(0.672/0.197,2);
const int Nbranches = 8;
const int nboots= 200;
const bool UseJack=1;
const int Njacks=30;
const int Nboots=200;
const double qu= 2.0/3.0;
const double qd= -1.0/3.0;
const double qs= qd;
const double qc= qu;
const int Upper_Limit_Time_Integral_strange= 300;
const int Upper_Limit_Time_Integral_charm= 300;
const int Upper_Limit_Time_Integral_light=300;
const double fm_to_inv_Gev= 1.0/0.1975;
const bool verbosity=1;
const double Nresonances=12;
const int Luscher_num_zeroes= 30;
const int npts_spline= 1000;
bool Use_OS_vector_current=false;
bool Use_Mpi_OS=false;
bool Include_light_disco= false;
bool Combined_fit_disco=false;
double Qfact= Combined_fit_disco?1.0:10.0/9.0;
const double m_phi= 3.0969;
const double m_Jpsi= 1.0195;
const double m_pi = 0.135;
const double m_rho= 0.775;
const double t0 = 0.4*fm_to_inv_Gev;
const double t1 = 1.0*fm_to_inv_Gev;
const double Delta= 0.15*fm_to_inv_Gev;

class ipar {

public:
  ipar() : V_light(0.0), V_light_err(0.0) {}

  
  double Mp, Mp_OS;
  double t,  L;
  double V_light, V_light_err;
  double fp;
   


};

class fit_par {

public:
  fit_par() {}
  fit_par(const Vfloat &par) {
    if((signed)par.size() != 6) crash("In class fit_par in fitting analytic representation of V(t)_light, class constructor Vfloat par has size != 5");
    Rd=par[0];
    Ed=par[1];
    Mrho=par[2];
    gpi=par[3];
    kappa= par[4];
    Z=par[5];
    
  }

  double Rd,Ed, Mrho, gpi,kappa, Z;



};



void Gm2() {

  omp_set_num_threads(1);

  //Init LL_functions;
  //find first  zeros of the Lusher functions
  Vfloat Luscher_zeroes;
  Zeta_function_zeroes(Luscher_num_zeroes, Luscher_zeroes);

  //############################################INTERPOLATE PHI FUNCTION AND DERIVATIVES#############################

  VVfloat phi_data, phi_der_data;
  Vfloat sx_int;
  Vfloat sx_der, dx_der;
  Vfloat Dz;

  for(int L_zero=0;L_zero<Nresonances+1;L_zero++) {
    double sx, dx;
    //interpolating between the Luscher_zero[L_zero-1] and Luscher_zero[L_zero];
    if(L_zero==0) { sx_int.push_back(0.0); sx=0.0;}
    else {sx=Luscher_zeroes[L_zero-1];  sx_int.push_back(sx);}
    dx= Luscher_zeroes[L_zero];
    phi_data.resize(L_zero+1);
    phi_der_data.resize(L_zero+1);
    phi_data[L_zero].push_back(L_zero==0?0.0:-M_PI/2.0);
    //divide interval into thousand points;
    double dz = (dx-sx)/npts_spline;
    Dz.push_back(dz);


    for(int istep=1;istep<=npts_spline-1;istep++) { double pt= sx+dz*istep; phi_data[L_zero].push_back( phi(sqrt(pt)));}

    phi_data[L_zero].push_back(M_PI/2.0);
    double sx_der_loc =  phi_der_for_back(sqrt(sx)+1e-14, 1);
    double dx_der_loc =  phi_der_for_back(sqrt(dx)-1e-14, -1);
    sx_der.push_back(sx_der_loc);
    dx_der.push_back(dx_der_loc);

    phi_der_data[L_zero].push_back(sx_der_loc);
    for(int istep=1;istep<=npts_spline-1;istep++) { double pt= sx+dz*istep; phi_der_data[L_zero].push_back( phi_der(sqrt(pt)));}
    phi_der_data[L_zero].push_back(dx_der_loc);
    
  }
    
   

  LL_functions LL(phi_data,phi_der_data,sx_der, dx_der, sx_int, Dz, Nresonances, Luscher_zeroes);
    
  //###########################################END INTERPOLATION PHI FUNCTION AND DERIVATIVES################################
  cout<<"####Spline for phi(x) and phi'(x) successfully generated!"<<endl;


  

  ///////   TEST //////////////

  double MR= 7.75260e-01;
  double gr= 5.90;
  double M = 0.134980;
  double ka = -3;
  double domega = 0.1;
  double omega_max = 10.0;
  double omega_start = 0.3;
  int nsteps = (int)omega_max/domega;

  for(int istep=0; istep<nsteps;istep++) {
    double omega_n = omega_start + istep*domega;
    cout<<"omega (GeV): "<<omega_n<<" "<<"  F_pi: "<<LL.F_pi_GS_mod(omega_n, MR, gr , M, ka)<<" V_pipi(L = \infty, t=1/omega):  "<< Qfact*LL.V_pipi_infL(1.0/omega_n, MR, gr, M, ka) <<endl;
  }
  exit(-1);
  
  
  
 
  

  //init Gaussian number generator
  GaussianMersenne GM(981832);
  string channel="";

  
  data_t  V_light_1, V_light_2, V_light_3, pt2_pion;
  data_t pt2_pion_strange, pt2_pion_charm;
  data_t  V_strange_1, V_strange_2, V_strange_3;
  data_t  V_charm_1, V_charm_2, V_charm_3;
  //to compute ZV
  data_t corr_A0P5;
  //to compute ZA
  data_t corr_A0P5_OS, corr_P5P5_OS;
  data_t disco_light;

  //Read data

  //Custom sorting for V_light to account for the two replica r0 and r1
  auto Sort_light_confs = [](string A, string B) {
		      string rA = A.substr(A.length()-2);
		      string rB = B.substr(B.length()-2);
		      int n1 = stoi(A.substr(A.length()-1));
		      int n2 = stoi(B.substr(B.length() -1));
		      if(rA==rB) return A<B;
		      else return n1<n2;
		      return A<B;
		    };
  string R_MODE="1";
  if(Use_OS_vector_current) R_MODE="2";
  //#################################END CUSTOM SORTING#################
  V_light_1.Read("../gm2_data/light", "mes_contr_2pts_ll_"+R_MODE, "V1V1", Sort_light_confs);
  V_light_2.Read("../gm2_data/light", "mes_contr_2pts_ll_"+R_MODE, "V2V2", Sort_light_confs);
  V_light_3.Read("../gm2_data/light", "mes_contr_2pts_ll_"+R_MODE, "V3V3", Sort_light_confs);
  pt2_pion.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  corr_A0P5.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "P5A0", Sort_light_confs);
  corr_A0P5_OS.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "P5A0", Sort_light_confs);
  corr_P5P5_OS.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs);
  if(Include_light_disco) disco_light.Read("../gm2_data/light", "disco", "", Sort_light_confs);
  
  
  V_strange_1.Read("../gm2_data/strange", "mes_contr_2pts_sl_ave", "V1V1");
  V_strange_2.Read("../gm2_data/strange", "mes_contr_2pts_sl_ave", "V2V2");
  V_strange_3.Read("../gm2_data/strange", "mes_contr_2pts_sl_ave", "V3V3");
  pt2_pion_strange.Read("../gm2_data/strange", "mes_contr_2pts_sl_ave", "P5P5");
  V_charm_1.Read("../gm2_data/charm", "mes_contr_2pts_cl_ave", "V1V1");
  V_charm_2.Read("../gm2_data/charm", "mes_contr_2pts_cl_ave", "V2V2");
  V_charm_3.Read("../gm2_data/charm", "mes_contr_2pts_cl_ave", "V3V3");
  pt2_pion_charm.Read("../gm2_data/charm", "mes_contr_2pts_cl_ave", "P5P5");
  

  int Nens_light= V_light_1.size;
  int Nens_strange = V_strange_1.size;
  int Nens_charm = V_charm_1.size;
  cout<<"N_ens light: "<<Nens_light<<endl;
  cout<<"N_ens strange: "<<Nens_strange<<endl;
  cout<<"N_ens charm: "<<Nens_charm<<endl;


  //create directories
  boost::filesystem::create_directory("../data/gm2");
  boost::filesystem::create_directory("../data/gm2/light");
  boost::filesystem::create_directory("../data/gm2/light/OS");
  boost::filesystem::create_directory("../data/gm2/light/tm");
  boost::filesystem::create_directory("../data/gm2/light/disco");
  boost::filesystem::create_directory("../data/gm2/strange");
  boost::filesystem::create_directory("../data/gm2/charm");

  //define distr_t_list to be used in chiral+continuum analysis

  distr_t_list agm2_strange(UseJack), agm2_charm(UseJack), agm2_light(UseJack), agm2_light_fit(UseJack), agm2_light_2L_fit(UseJack), agm2_light_Lprime_fit(UseJack)  , agm2_light_infL_fit(UseJack), agm2_light_W(UseJack);
  vector<vector<fit_par>> par_list_anal_repr;//only used for light quark contribution

  distr_t_list MV_fit_light(UseJack), MV_fit_strange(UseJack), MV_fit_charm(UseJack);
  distr_t_list ZV_fit_light(UseJack), ZV_fit_strange(UseJack), ZV_fit_charm(UseJack);
  distr_t_list Mpi_fit(UseJack), Mpi_OS_fit(UseJack),  fp_fit(UseJack);
  distr_t_list Zv_fit(UseJack), Za_fit(UseJack), Zp_ov_Zs_fit(UseJack);
  Vfloat L_list, a_list, ml_list;
  distr_t_list a_distr_list(UseJack);
  Vfloat L_strange_list, a_strange_list, ml_strange_list;
  Vfloat L_charm_list, a_charm_list, ml_charm_list;

  //define lambda for convolution with kernel
  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
  auto exp_MV = [&](double Mv, double t, double size) -> double { return exp(-Mv*t);};
  
  
  //strange
  channel="s";
  for(int i_ens=0;i_ens<Nens_strange;i_ens++) {
    
  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = V_strange_1.nrows[i_ens];


  //resample lattice spacing
  distr_t a_distr(UseJack), ZA_distr(UseJack);
  LatticeInfo L_info;
  L_info.LatInfo_new_ens(V_strange_1.Tag[i_ens]);
  //generate jackknife sample of input parameters
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      a_distr.distr.push_back( fm_to_inv_Gev*( L_info.a + GM()*L_info.a_err*(1.0/sqrt(Njacks-1.0))));
      ZA_distr.distr.push_back(  L_info.Za + GM()*L_info.Za_err*(1.0/sqrt(Njacks-1.0)));
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      a_distr.distr.push_back( fm_to_inv_Gev*( L_info.a + GM()*L_info.a_err));
      ZA_distr.distr.push_back(  L_info.Za + GM()*L_info.Za_err);
    }
  }
  
  //push_back lattice info
  L_strange_list.push_back(L_info.L);
  a_strange_list.push_back(L_info.a);
  ml_strange_list.push_back(L_info.ml);
  
 
  distr_t_list  V_strange_1_distr, V_strange_2_distr, V_strange_3_distr;
  distr_t_list  V_strange_distr, MV_strange_distr, ZV_strange_distr;
  distr_t MV_strange, ZV_strange;

  //Analyze correlators
  if(Corr.Nt > 48) {
  Corr.Tmin=15;
  Corr.Tmax=30;
  }
  else { Corr.Tmin= 12; Corr.Tmax=20;}
  
 
  //strange
  V_strange_1_distr = Corr.corr_t(V_strange_1.col(0)[i_ens], "../data/gm2/strange/corr_1_"+V_strange_1.Tag[i_ens]+".dat");
  V_strange_2_distr = Corr.corr_t(V_strange_2.col(0)[i_ens], "../data/gm2/strange/corr_2_"+V_strange_2.Tag[i_ens]+".dat");
  V_strange_3_distr = Corr.corr_t(V_strange_3.col(0)[i_ens], "../data/gm2/strange/corr_3_"+V_strange_3.Tag[i_ens]+".dat");
 

 
  //sum over the Lorenz indices of the e.m. current
  V_strange_distr= (pow(qs,2)/3.0)*(V_strange_1_distr+ V_strange_2_distr + V_strange_3_distr);
 
  // print summed correlators to file
  Print_To_File({}, {V_strange_distr.ave(), V_strange_distr.err()}, "../data/gm2/strange/corr_sum_"+V_strange_1.Tag[i_ens]+".dat.t", "", "");
  
  //extract effective masses, overlap from V and fit

 
  if(V_strange_1.Tag[i_ens].substr(1,1) == "C") { Corr.Tmin=15; Corr.Tmax=32;}
  else if(V_strange_1.Tag[i_ens].substr(1,1) == "B") { Corr.Tmin=15; Corr.Tmax=20;}
  else if(V_strange_1.Tag[i_ens].substr(1,1) == "A") {
    if( L_info.L == 24) { Corr.Tmin= 13; Corr.Tmax = 17;}
    else { Corr.Tmin= 15; Corr.Tmax=20;}
  }
  else crash("Ensemble tag not valid");
 
  //strange
  MV_strange_distr= Corr.effective_mass_t(V_strange_distr, "../data/gm2/strange/MV_mass_"+V_strange_1.Tag[i_ens]+".dat");
  ZV_strange_distr= Corr.residue_t(V_strange_distr, "../data/gm2/strange/ZV_overlap_"+V_strange_1.Tag[i_ens]+".dat");
  MV_strange = Corr.Fit_distr(MV_strange_distr);
  ZV_strange = Corr.Fit_distr(ZV_strange_distr);

  //push_back MV_strange and ZV_strange
  MV_fit_strange.distr_list.push_back(MV_strange);
  ZV_fit_strange.distr_list.push_back(ZV_strange);




  int Tdata_min= 8;
  int Tdata_max = Corr.Nt/2.0 -4;
  int Tdata_fit = (Corr.Tmax+Corr.Tmin)/2;
  
  //compute kernel distribution
  distr_t_list Kernel_distr = distr_t_list::f_of_distr(K,MV_strange/m_phi, Upper_Limit_Time_Integral_strange+1);
  //compute exp(-Mv*t) distribution
  distr_t_list exp_MVs = distr_t_list::f_of_distr(exp_MV, MV_strange, Upper_Limit_Time_Integral_strange+1);

  //Print single-exponential prediction to file
  Print_To_File({}, {(exp_MVs*(ZV_strange/(2.0*MV_strange))).ave(), (exp_MVs*(ZV_strange/(2.0*MV_strange))).err()}, "../data/gm2/strange/corr_gsd_sum_"+V_strange_1.Tag[i_ens]+".dat.t", "", "");

  distr_t_list agm2_distr_Tdata(UseJack);
  Vfloat Tdata_vec;
  
  for(int Tdata=Tdata_min;Tdata<Tdata_max;Tdata++) {
    //compute 4\pia^2 using lattice data up to Tdata (included)
   
    distr_t agm2(UseJack, UseJack?Njacks:Nboots);
    
    for(int t=1;t<=Upper_Limit_Time_Integral_strange;t++) {
      if(t<=Tdata) agm2 = agm2 + 4.0*pow(alpha,2)*V_strange_distr.distr_list[t]*Kernel_distr.distr_list[t];
      else agm2= agm2 + 4.0*pow(alpha,2)*(ZV_strange/(2.0*MV_strange))*exp_MVs.distr_list[t]*Kernel_distr.distr_list[t];
    }
    Tdata_vec.push_back((double)Tdata);
    agm2_distr_Tdata.distr_list.push_back(ZA_distr*ZA_distr*agm2);
    if(Tdata==Tdata_fit) agm2_strange.distr_list.push_back(ZA_distr*ZA_distr*agm2);
  }
  //print to file
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata.ave(), agm2_distr_Tdata.err()}, "../data/gm2/strange/agm2_Tdata_"+V_strange_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m agm2_err");
  }

  cout<<"strange quark correlator analyzed!"<<endl;

  //charm
  channel="c";
  for(int i_ens=0;i_ens<Nens_charm;i_ens++) {
    
  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = V_charm_1.nrows[i_ens];


  //resample lattice spacing
  distr_t a_distr(UseJack), ZA_distr(UseJack);
  LatticeInfo L_info;
  L_info.LatInfo_new_ens(V_charm_1.Tag[i_ens]);
  //generate jackknife sample of input parameters
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      a_distr.distr.push_back( fm_to_inv_Gev*( L_info.a + GM()*L_info.a_err*(1.0/sqrt(Njacks-1.0))));
      ZA_distr.distr.push_back(  L_info.Za + GM()*L_info.Za_err*(1.0/sqrt(Njacks-1.0)));
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      a_distr.distr.push_back( fm_to_inv_Gev*( L_info.a + GM()*L_info.a_err));
      ZA_distr.distr.push_back(  L_info.Za + GM()*L_info.Za_err);
    }
  }


  //push_back lattice info
  L_charm_list.push_back(L_info.L);
  a_charm_list.push_back(L_info.a);
  ml_charm_list.push_back(L_info.ml);
  
  distr_t_list  V_charm_1_distr, V_charm_2_distr, V_charm_3_distr;
  distr_t_list  V_charm_distr, MV_charm_distr, ZV_charm_distr;
  distr_t MV_charm , ZV_charm;

  //Analyze correlators

 
    

  if(Corr.Nt > 48) {
  Corr.Tmin=15;
  Corr.Tmax=30;
  }
  else { Corr.Tmin= 12; Corr.Tmax=20;}
  
  //charm
  V_charm_1_distr = Corr.corr_t(V_charm_1.col(0)[i_ens], "../data/gm2/charm/corr_1_"+V_charm_1.Tag[i_ens]+".dat");
  V_charm_2_distr = Corr.corr_t(V_charm_2.col(0)[i_ens], "../data/gm2/charm/corr_2_"+V_charm_2.Tag[i_ens]+".dat");
  V_charm_3_distr = Corr.corr_t(V_charm_3.col(0)[i_ens], "../data/gm2/charm/corr_3_"+V_charm_3.Tag[i_ens]+".dat");

  //sum over the Lorenz indices of the e.m. current
  V_charm_distr= (pow(qc,2)/3.0)*(V_charm_1_distr+ V_charm_2_distr + V_charm_3_distr);


  // print summed correlators to file
  Print_To_File({}, {V_charm_distr.ave(), V_charm_distr.err()}, "../data/gm2/charm/corr_sum_"+V_charm_1.Tag[i_ens]+".dat.t", "", "");


 //extract effective masses, overlap from V and fit

  //charm

  if(V_charm_1.Tag[i_ens].substr(1,1) == "C") { Corr.Tmin=33; Corr.Tmax=50;}
  else if(V_charm_1.Tag[i_ens].substr(1,1) == "B") { Corr.Tmin=25; Corr.Tmax=40;}
  else if(V_charm_1.Tag[i_ens].substr(1,1) == "A") {
    if( L_info.L == 24) { Corr.Tmin= 17; Corr.Tmax = 22;}
    else { Corr.Tmin= 22; Corr.Tmax=30;}
  }
  else crash("Ensemble tag not valid");
  MV_charm_distr= Corr.effective_mass_t(V_charm_distr, "../data/gm2/charm/MV_mass_"+V_charm_1.Tag[i_ens]+".dat");
  ZV_charm_distr= Corr.residue_t(V_charm_distr, "../data/gm2/charm/ZV_overlap_"+V_charm_1.Tag[i_ens]+".dat");
  MV_charm = Corr.Fit_distr(MV_charm_distr);
  ZV_charm = Corr.Fit_distr(ZV_charm_distr);
  //push_back MV_charm and ZV_charm
  MV_fit_charm.distr_list.push_back(MV_charm);
  ZV_fit_charm.distr_list.push_back(ZV_charm);
 


  int Tdata_min= 8;
  int Tdata_max = Corr.Nt/2.0 -4;
  int Tdata_fit = (Corr.Tmax+Corr.Tmin)/2;
  
  //compute kernel distribution
  distr_t_list Kernel_distr = distr_t_list::f_of_distr(K,MV_charm/m_Jpsi, Upper_Limit_Time_Integral_strange+1);
  //compute exp(-Mv*t) distribution
  distr_t_list exp_MVc = distr_t_list::f_of_distr(exp_MV, MV_charm, Upper_Limit_Time_Integral_charm+1);

  //Print single-exponential prediction to file
  Print_To_File({}, {(exp_MVc*(ZV_charm/(2.0*MV_charm))).ave(), (exp_MVc*(ZV_charm/(2.0*MV_charm))).err()}, "../data/gm2/charm/corr_gsd_sum_"+V_charm_1.Tag[i_ens]+".dat.t", "", "");

  distr_t_list agm2_distr_Tdata(UseJack);
  Vfloat Tdata_vec;
  
  for(int Tdata=Tdata_min;Tdata<Tdata_max;Tdata++) {
    //compute 4\pia^2 using lattice data up to Tdata (included)
   
    distr_t agm2(UseJack, UseJack?Njacks:Nboots);
    
    for(int t=1;t<=Upper_Limit_Time_Integral_charm;t++) {
      if(t<=Tdata) agm2 = agm2 + 4.0*pow(alpha,2)*V_charm_distr.distr_list[t]*Kernel_distr.distr_list[t];
      else agm2= agm2 + 4.0*pow(alpha,2)*(ZV_charm/(2.0*MV_charm))*exp_MVc.distr_list[t]*Kernel_distr.distr_list[t];
    }
    Tdata_vec.push_back((double)Tdata);
    agm2_distr_Tdata.distr_list.push_back(ZA_distr*ZA_distr*agm2);
    if(Tdata==Tdata_fit) agm2_charm.distr_list.push_back(ZA_distr*ZA_distr*agm2);
  }
  //print to file
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata.ave(), agm2_distr_Tdata.err()}, "../data/gm2/charm/agm2_Tdata_"+V_charm_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m agm2_err");
  }


  cout<<"charm quark correlator analyzed!"<<endl;

  //light
  channel="l";
  string light_suff=Use_OS_vector_current?"OS":"tm";
  for(int i_ens=0;i_ens<Nens_light;i_ens++) {
    
  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = V_light_1.nrows[i_ens];

  //resample lattice spacing
  distr_t a_distr(UseJack), ZA_distr(UseJack), ZV_distr(UseJack);
  LatticeInfo L_info;
  L_info.LatInfo_new_ens(V_light_1.Tag[i_ens]);
  //generate jackknife sample of input parameters
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      a_distr.distr.push_back( fm_to_inv_Gev*( L_info.a + GM()*L_info.a_err*(1.0/sqrt(Njacks-1.0))));
      ZA_distr.distr.push_back(  L_info.Za + GM()*L_info.Za_err*(1.0/sqrt(Njacks-1.0)));
      ZV_distr.distr.push_back(  L_info.Zv + GM()*L_info.Zv_err*(1.0/sqrt(Njacks-1.0)));
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      a_distr.distr.push_back( fm_to_inv_Gev*( L_info.a + GM()*L_info.a_err));
      ZA_distr.distr.push_back(  L_info.Za + GM()*L_info.Za_err);
      ZV_distr.distr.push_back( L_info.Zv + GM()*L_info.Zv_err);
    }
  }



  
  distr_t_list  V_light_1_distr, V_light_2_distr, V_light_3_distr, disco_distr, conn_disco_distr;
  distr_t_list  V_light_distr, MV_light_distr, ZV_light_distr, MV_light_disco_distr;
  distr_t_list Mpi_distr, Mpi_OS_distr, fp_distr, overlap_P5P5_distr, overlap_P5P5_OS_distr, pion_corr;
  distr_t_list Omega_distr;
  distr_t_list A0P5_distr, A0P5_OS_distr;
  distr_t MV_light, ZV_light, Mpi, Mpi_OS, fp;
  distr_t_list P5P5_OS_distr,  RV, RA,RA0, ratio_P5P5_overlap_OS_tm, Zp_ov_Zs_distr;
  distr_t Zv, Za, Zp_ov_Zs, fp_ov_Mpi;
  

  //Analyze correlators
   if(V_light_1.Tag[i_ens].substr(1,1)=="A") {Corr.Tmin=14;Corr.Tmax=19;}
   else if(V_light_1.Tag[i_ens].substr(1,1)=="B") {Corr.Tmin=23; Corr.Tmax=33;}
   else if(V_light_1.Tag[i_ens].substr(1,1)=="C") {Corr.Tmin=33; Corr.Tmax=43;}
   else if(V_light_1.Tag[i_ens].substr(1,1)=="D") {Corr.Tmin=30; Corr.Tmax=50;}
   else crash("Cannot find ensemble: "+V_light_1.Tag[i_ens]+" while modifying fit interval for MV_light");

  //define lambda function to compute function of distributions
  auto sqr= [=](double a, double b) {return sqrt(a);};
  auto SINH= [](double m) -> double  {return sinh(m);};
  auto SINH2 = [](double m, double t) -> double {return sinh(m);};

  
  //vector light sector
  V_light_1_distr = Corr.corr_t(V_light_1.col(0)[i_ens], "../data/gm2/light/"+light_suff+"/corr_1_"+V_light_1.Tag[i_ens]+".dat");
  V_light_2_distr = Corr.corr_t(V_light_2.col(0)[i_ens], "../data/gm2/light/"+light_suff+"/corr_2_"+V_light_2.Tag[i_ens]+".dat");
  V_light_3_distr = Corr.corr_t(V_light_3.col(0)[i_ens], "../data/gm2/light/"+light_suff+"/corr_3_"+V_light_3.Tag[i_ens]+".dat");
  //sum over the Lorenz indices of the e.m. current
  V_light_distr= ((pow(qu,2)+pow(qd,2))/3.0)*(V_light_1_distr+ V_light_2_distr + V_light_3_distr);
  if(Include_light_disco) {
    disco_distr = Corr.corr_t(disco_light.col(0)[i_ens], "");
  
    MV_light_disco_distr= Corr.effective_mass_t(conn_disco_distr, "../data/gm2/light/disco/MV_"+V_light_1.Tag[i_ens]+".dat");
    Omega_distr = Corr.effective_mass_t( 2.0*V_light_1_distr + disco_distr, "../data/gm2/light/disco/Omega_"+V_light_1.Tag[i_ens]+".dat");
    disco_distr = disco_distr*(pow(qu+qd,2)/4.0);
    conn_disco_distr = V_light_distr + disco_distr;
   
  }
  MV_light_distr= Corr.effective_mass_t(V_light_distr, "../data/gm2/light/"+light_suff+"/MV_"+V_light_1.Tag[i_ens]+".dat");
  ZV_light_distr= Corr.residue_t(V_light_distr, "../data/gm2/light/"+light_suff+"/ZV_overlap_"+V_light_1.Tag[i_ens]+"dat");


  //tm pion sector
  pion_corr= Corr.corr_t(pt2_pion.col(0)[i_ens], "");
  Mpi_distr= Corr.effective_mass_t(pion_corr, "../data/gm2/light/Mpi_"+V_light_1.Tag[i_ens]+".dat");
  overlap_P5P5_distr = Corr.residue_t(pion_corr, "");
  fp_distr = Corr.decay_constant_t(pow(2.0*L_info.ml,2)*pion_corr, "../data/gm2/light/"+light_suff+"/fp_"+V_light_1.Tag[i_ens]+".dat");


  //OS sector
  P5P5_OS_distr = Corr.corr_t(corr_P5P5_OS.col(0)[i_ens], "");
  Mpi_OS_distr= Corr.effective_mass_t(P5P5_OS_distr, "../data/gm2/light/Mpi_OS_"+V_light_1.Tag[i_ens]+".dat");
  overlap_P5P5_OS_distr= Corr.residue_t(P5P5_OS_distr, "");



  //take ratio between OS and tm pion amplitude to compute Zp/Zs RC.
  ratio_P5P5_overlap_OS_tm= overlap_P5P5_OS_distr/overlap_P5P5_distr;
  Zp_ov_Zs_distr = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm);



  //antysymmetrize w.r.t. t -> T-t for A0P5 correlators
  //Corr.Reflection_sign = -1;
  A0P5_distr= Corr.corr_t(corr_A0P5.col(0)[i_ens], "");
  A0P5_OS_distr = Corr.corr_t(corr_A0P5_OS.col(0)[i_ens], "");

  
  //restore symmetrization
  // Corr.Reflection_sign = 1;

  //compute RV (estimator for Zv)


  RV= 2.0*L_info.ml*pion_corr/distr_t_list::derivative(A0P5_distr, 0); //central derivative

 
  

  //extract effective masses, overlap from V and fit
  int Tmin_old=Corr.Tmin;
  int Tmax_old=Corr.Tmax;
  //define restricted time intervals since Rho plateaux are tricky
  if(V_light_1.Tag[i_ens].substr(1,1)=="C" || V_light_1.Tag[i_ens].substr(1,1) == "D") {Corr.Tmin=14;Corr.Tmax=19;}
  else if(V_light_1.Tag[i_ens].substr(1,1)=="B" || V_light_1.Tag[i_ens].substr(1,1) =="A") {Corr.Tmin=16; Corr.Tmax=25;}
  else crash("Cannot find ensemble: "+V_light_1.Tag[i_ens]+" while modifying fit interval for MV_light");
  //################################################

  
  MV_light = Corr.Fit_distr(MV_light_distr);
  ZV_light = Corr.Fit_distr(ZV_light_distr);
  //push_back MV_light and ZV_light
  MV_fit_light.distr_list.push_back(MV_light);
  ZV_fit_light.distr_list.push_back(ZV_light);


  //restore old time intervals
  Corr.Tmin=Tmin_old;
  Corr.Tmax=Tmax_old;
  //########################

  
  Mpi=Corr.Fit_distr(Mpi_distr);
  Mpi_OS= Corr.Fit_distr(Mpi_OS_distr);
  fp= Corr.Fit_distr(fp_distr);
  fp_ov_Mpi= fp/Mpi;
  Zp_ov_Zs = Corr.Fit_distr(Zp_ov_Zs_distr);
  Zv= Corr.Fit_distr(RV);
  RA0 = (Mpi_OS_distr/Mpi_distr)*(distr_t_list::f_of_distr_list(SINH2, Mpi_OS_distr)/distr_t_list::f_of_distr_list(SINH2, Mpi_distr));
  RA = 2.0*L_info.ml*(P5P5_OS_distr/distr_t_list::derivative(A0P5_OS_distr, 0))*(Mpi_OS/Mpi)*(distr_t::f_of_distr(SINH, Mpi_OS)/distr_t::f_of_distr(SINH, Mpi))*(1.0/Zp_ov_Zs);
  Za = Corr.Fit_distr(RA);

  //#################################################################


  
  //Print To File additional observables
  // print summed connected correlators to file
  Print_To_File({}, {V_light_distr.ave(), V_light_distr.err(), Use_OS_vector_current?(Zv*Zv*V_light_distr).ave():(Za*Za*V_light_distr).ave(),  Use_OS_vector_current?(Zv*Zv*V_light_distr).err():(Za*Za*V_light_distr).err()}, "../data/gm2/light/"+light_suff+"/corr_sum_"+V_light_1.Tag[i_ens]+".dat.t", "", "");
  //Print disco and iscoscalar conn disco combination
  if(Include_light_disco) Print_To_File({}, {disco_distr.ave(), disco_distr.err(), conn_disco_distr.ave(), conn_disco_distr.err()}, "../data/gm2/light/disco/disc_"+V_light_1.Tag[i_ens]+".dat.t","",""); 
  //print RV
  Print_To_File({}, {RV.ave(), RV.err()}, "../data/gm2/light/RV_"+V_light_1.Tag[i_ens]+".dat.t", "", "");
  //print RA
  Print_To_File({}, {RA.ave(), RA.err(), RA0.ave(), RA0.err()}, "../data/gm2/light/RA_"+V_light_1.Tag[i_ens]+".dat.t", "", "");
  //print Zp_ov_Zs
  Print_To_File({}, {Zp_ov_Zs_distr.ave(), Zp_ov_Zs_distr.err()}, "../data/gm2/light/Zp_ov_Zs_"+V_light_1.Tag[i_ens]+".dat.t", "", "");


  cout<<"//#########################ENSEMBLE INFO: "<<endl;
  cout<<"Ensemble Tag: "<<V_light_1.Tag[i_ens]<<endl; 
  cout<<"L: "<<L_info.L<<" T: "<<L_info.T<<endl;
  cout<<"lat spacing: "<<L_info.a<<"+- "<<L_info.a_err<<" (fm)"<<endl;
  cout<<"ml: "<<L_info.ml<<endl;
  cout<<"Mpi: "<<Mpi.ave()<<" +- "<<Mpi.err()<<endl;
  cout<<"Mpi*L: "<<Mpi.ave()*L_info.L<<" +- "<<Mpi.err()*L_info.L<<endl;
  cout<<"Mpi OS: "<<Mpi_OS.ave()<<" +- "<<Mpi_OS.err()<<endl;
  cout<<"fp: "<<fp.ave()<<" +- "<<fp.err()<<endl;
  cout<<"fp/Mpi: "<<fp_ov_Mpi.ave()<<" +- "<<fp_ov_Mpi.err()<<" Phys point: "<<0.966<<endl;
  cout<<"Zp/Zs: "<<Zp_ov_Zs.ave()<<" +- "<<Zp_ov_Zs.err()<<endl;
  cout<<"Za: "<<Za.ave()<<" +- "<<Za.err()<<endl;
  cout<<"Zv: "<<Zv.ave()<<" +- "<<Zv.err()<<endl;
  cout<<"MV: "<<MV_light.ave()<<" +- "<<MV_light.err()<<endl;
 

  
  //push_back Mpi Mpi_OS and fp, Rcs a and L
  Mpi_fit.distr_list.push_back(Mpi);
  Mpi_OS_fit.distr_list.push_back(Mpi_OS);
  fp_fit.distr_list.push_back(fp);
  Za_fit.distr_list.push_back(Za);
  Zv_fit.distr_list.push_back(Zv);
  Zp_ov_Zs_fit.distr_list.push_back(Zp_ov_Zs);
  ml_list.push_back(L_info.ml);
  a_list.push_back(L_info.a);
  L_list.push_back((double)L_info.L);
  

  


  int Tdata_min= 8;
  int Tdata_max = Corr.Nt/2.0 -4;
  int Tdata_fit = (Corr.Tmax+Corr.Tmin)/2;
  
  //compute kernel distribution
  distr_t_list Kernel_distr = distr_t_list::f_of_distr(K,Mpi/m_pi, Upper_Limit_Time_Integral_light+1);
  //compute exp(-Mv*t) distribution
  distr_t_list exp_MVl = distr_t_list::f_of_distr(exp_MV, MV_light, Upper_Limit_Time_Integral_light+1);


  //Print single-exponential prediction to file
  Print_To_File({}, {(exp_MVl*(ZV_light/(2.0*MV_light))).ave(), (exp_MVl*(ZV_light/(2.0*MV_light))).err()}, "../data/gm2/light/"+light_suff+"/corr_gsd_sum_"+V_light_1.Tag[i_ens]+".dat.t", "", "");

  distr_t_list agm2_distr_Tdata(UseJack);
  Vfloat Tdata_vec;
  
  for(int Tdata=Tdata_min;Tdata<Tdata_max;Tdata++) {
    //compute 4\pia^2 using lattice data up to Tdata (included)
   
    distr_t agm2(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2 to zero by default
    
    
    for(int t=1;t<=Upper_Limit_Time_Integral_light;t++) {
      if(t<=Tdata) agm2 = agm2 + 4.0*pow(alpha,2)*V_light_distr.distr_list[t]*Kernel_distr.distr_list[t];
      else agm2= agm2 + 4.0*pow(alpha,2)*(ZV_light/(2.0*MV_light))*exp_MVl.distr_list[t]*Kernel_distr.distr_list[t];
    }
    Tdata_vec.push_back((double)Tdata);
    if(!Use_OS_vector_current) agm2 = agm2*Za*Za;
    else agm2 = agm2*Zv*Zv;
    agm2_distr_Tdata.distr_list.push_back(agm2);
    if(Tdata==Tdata_fit) agm2_light.distr_list.push_back(agm2);
  }
  //print to file
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata.ave(), agm2_distr_Tdata.err()}, "../data/gm2/light/"+light_suff+"/agm2_Tdata_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m agm2_err");


  /*

  //#######################  INTERMEDIATE WINDOW ###################################

  distr_t agm2_W(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default

  for(int t=1; t< Corr.Nt; t++) {
    agm2_W = agm2_W + + 4.0*pow(alpha,2)*V_light_distr.distr_list[t]*Kernel_distr.distr_list[t];



  }


  */


  //################################################################################


  
  //#################################################################################
  //#################################################################################
  //#################################################################################
  //#################################################################################
  //fit lattice data using analytic representation for V(t)_light



  bootstrap_fit<fit_par,ipar> bf(Njacks);
  bf.set_warmup_lev(3); //sets warmup

  //######################FIT INTERVALS####################
  //##############INCLUDE only times t for which dC(t)/C(t) < 0.1#####################

  //Tfit_min always starts at 0.2fm
  //double t02fm= 0.2/L_info.a +1.0;
  int Tfit_min= 5;
  if(Use_OS_vector_current) Tfit_min=6;
  if(V_light_1.Tag[i_ens].substr(1,1) == "D") Tfit_min++;
  int Tfit_max= Tfit_min;
  bool Found_error_less_10_percent=false;
  while(!Found_error_less_10_percent && Tfit_max < Corr.Nt/2) {
    distr_t Z= Use_OS_vector_current?Zv:Za;
    if( (Z*V_light_distr.distr_list[Tfit_max]).err()/fabs( (Z*V_light_distr.distr_list[Tfit_max]).ave()) <  0.1) Tfit_max++;
    else Found_error_less_10_percent=true;
  }
  int Tfit_points = Tfit_max +1 - Tfit_min;
  //###################################################################################

  
  bf.Set_number_of_measurements(Tfit_points);
  bf.Set_verbosity(verbosity);

  cout<<"FITTING WITH PIPI+DUAL REPRESENTATION ENSEMBLE: "<<V_light_1.Tag[i_ens]<<endl;
  cout<<"Tmin: "<<Tfit_min<<", Tmax: "<<Tfit_max<<endl;

  //initial guesses are based on the results on the old ETMC confs
  bf.Add_par("Rd", 1.1, 0.1);
  bf.Set_limits("Rd", 0.7, 1.8);
  if(!Use_Mpi_OS) {
  bf.Add_par("Ed", 2.5, 0.2);
  bf.Set_limits("Ed", 0.7, 5.0);
  bf.Add_par("Mrho",5.5, 0.1);
  bf.Set_limits("Mrho", 3.0, 7.0);
  }
  else {
  bf.Add_par("Ed", 1.3, 0.2);
  bf.Set_limits("Ed", 0.5, 3.0);
  bf.Add_par("Mrho",3.0, 0.1);
  bf.Set_limits("Mrho", 2.2, 5);
  }
  bf.Add_prior_par("gpi", 1.0, 0.01);
  bf.Add_prior_par("kappa", -3.0, 0.1);
  //bf.Set_limits("gpi",0.75, 1.5);
  //bf.Set_limits("kappa", -5.0, 1.0);
  if(!Use_OS_vector_current)  bf.Add_prior_par("Z", Za.ave(), Za.err());
  else bf.Add_prior_par("Z", Zv.ave(), Zv.err());
  bf.Fix_n_release("Z");

  bf.Fix_par("gpi",1.0);
  bf.Fix_par("kappa", -3.0);

  map<pair<pair<double,double>, pair<double,double>>,Vfloat> Energy_lev_list;
  

  bf.ansatz =  [&](const fit_par &p, const ipar &ip) -> double {


		 double Pi_M = Use_Mpi_OS?ip.Mp_OS:ip.Mp;

		 double GPI = p.gpi*p.Mrho*Pi_M/ip.fp;
		 Vfloat Knpp;
		 pair<double,double> Mass_par= make_pair(p.Mrho*Pi_M, Pi_M);
		 pair<double,double> Couplings = make_pair(GPI, p.kappa);
		 pair< pair<double,double>,pair<double, double>> input_pars = make_pair( Mass_par, Couplings);
		 map<pair<pair<double,double>, pair<double, double>>,Vfloat>::iterator it;
		 it= Energy_lev_list.find(input_pars);
		 if(it != Energy_lev_list.end()) Knpp= it->second;
		 else {
		   LL.Find_pipi_energy_lev(ip.L,p.Mrho*Pi_M, GPI, Pi_M, p.kappa, Knpp);
		   //add to Energy_lev_list
		   Energy_lev_list.insert( make_pair(input_pars, Knpp));
		 }
		
		 return Qfact*LL.V_pipi(ip.t, ip.L, p.Mrho*Pi_M, GPI, Pi_M, p.kappa, Knpp) + LL.Vdual(ip.t, p.Mrho*Pi_M, p.Ed*Pi_M, p.Rd);
  };
  bf.measurement = [&](const fit_par& p,const ipar& ip) -> double {
		     return pow(p.Z,2)*ip.V_light;
  };
  bf.error =  [&](const fit_par& p,const ipar &ip) -> double {
		return pow(p.Z,2)*ip.V_light_err;
  };

  //fill the data
  vector<vector<ipar>> data(Njacks);
  //allocate space for output result
  boot_fit_data<fit_par> Bt_fit;

  

  for(auto &data_iboot: data) data_iboot.resize(Tfit_points);
  
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int t=Tfit_min;t<Tfit_points+Tfit_min;t++) {
	int tt=t-Tfit_min;
	data[ijack][tt].V_light = Combined_fit_disco?conn_disco_distr.distr_list[t].distr[ijack]:V_light_distr.distr_list[t].distr[ijack];
	data[ijack][tt].V_light_err = Combined_fit_disco?conn_disco_distr.distr_list[t].err():V_light_distr.distr_list[t].err();
	data[ijack][tt].L = L_info.L;
	data[ijack][tt].Mp = Mpi.distr[ijack];
	data[ijack][tt].Mp_OS = Mpi_OS.distr[ijack];
	data[ijack][tt].t = t;
	data[ijack][tt].fp = fp.distr[ijack];
      } 
    }

    
    //add prior values
    for(int ijack=0;ijack<Njacks;ijack++) {
     bf.ib= &ijack;
     bf.Append_to_prior("gpi", 1.0, 0.2);
     bf.Append_to_prior("kappa", -3.0, 3.0);
     bf.Append_to_prior("Z", Use_OS_vector_current?Zv.distr[ijack]:Za.distr[ijack], Use_OS_vector_current?Zv.err():Za.err());
    }

    
    //append
    bf.Append_to_input_par(data);

    
    
    //fit
    Bt_fit= bf.Perform_bootstrap_fit();
    
    Bt_fit.ch2_ave();

    //###################################################
    //print fitted func
    
    double step_size_time= 0.1;
    int nsteps= (int)((1.0/step_size_time)*Corr.Nt); //500 points in time
    cout<<"Printing fit function with: istep: "<<step_size_time<<" and nstep: "<<nsteps<<endl;
    distr_t_list func(UseJack, nsteps);
    Vfloat times;

    //generate data
    for(int istep=0; istep<nsteps;istep++) {
      ipar my_ipar;
      my_ipar.t = step_size_time*(1.0 + istep);
      my_ipar.L = (double)L_info.L;
      times.push_back(my_ipar.t);
      for(int ijack=0;ijack<Njacks;ijack++) {
	my_ipar.Mp=Mpi.distr[ijack];
	my_ipar.Mp_OS=Mpi_OS.distr[ijack];
	my_ipar.fp = fp.distr[ijack];
	func.distr_list[istep].distr.push_back( bf.ansatz( Bt_fit.par[ijack], my_ipar )/pow(Bt_fit.par[ijack].Z,2));
      }
    }

    //print to file
    Print_To_File({}, {times, func.ave(), func.err()}, "../data/gm2/light/"+light_suff+"/Corr_anal_rep_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id_col  time  val   err  (fit interval: ["+to_string(Tfit_min)+","+to_string(Tfit_points+Tfit_min-1)+"])  a = "+to_string_with_precision(L_info.a,6)+" +- "+to_string_with_precision(L_info.a_err, 6)+" fm");
    //####################################################
    
    
    //push_back jack distribution of fitted params
    par_list_anal_repr.push_back(Bt_fit.par);
    cout<<"SIZE OF ENERGY LEVEL MAP IS: "<<Energy_lev_list.size()<<endl;
    cout<<"#########################END ENSEMBLE FIT"<<endl;

  }

  cout<<"light quark correlator analyzed!"<<endl;






  


  cout<<"####### Reconstructing agm2_light using fit paramters Ed, Rd, Mrho, gpi....."<<endl;
  distr_t_list Edual(UseJack), Rdual(UseJack), Mrho(UseJack), gpi(UseJack), Kappa(UseJack);
  for(int i_ens=0; i_ens<Nens_light;i_ens++) {
    cout<<"######### ENSEMBLE: "<<V_light_1.Tag[i_ens]<<endl;
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_light_1.Tag[i_ens]);
  

    
    distr_t Ed(UseJack), Rd(UseJack), Mr(UseJack), g(UseJack), kap(UseJack), agm2_dual(UseJack), mp(UseJack), agm2_2L_dual(UseJack), agm2_infL_dual(UseJack), agm2_Lprime_dual(UseJack) ,  a_distr(UseJack);
    for(int ijack=0;ijack<Njacks;ijack++) {
    //retrieve fit_parameters
    fit_par my_fit_par = par_list_anal_repr[i_ens][ijack];
    Ed.distr.push_back(my_fit_par.Ed);
    Rd.distr.push_back(my_fit_par.Rd);
    Mr.distr.push_back(my_fit_par.Mrho);
    g.distr.push_back(my_fit_par.gpi);
    kap.distr.push_back(my_fit_par.kappa);
    a_distr.distr.push_back( fm_to_inv_Gev*(L_info.a + GM()*L_info.a_err*(1.0/sqrt(Njacks-1.0))));

   

    double Mp,fp, csi;
    double L= 1.0*L_list[i_ens];
    fp= fp_fit[i_ens].distr[ijack];
    //double csi = pow(Mp,2)/pow(4.0*M_PI*fp,2);
    if(!Use_Mpi_OS) Mp=Mpi_fit.distr_list[i_ens].distr[ijack];
    else Mp=Mpi_OS_fit.distr_list[i_ens].distr[ijack];
    mp.distr.push_back(Mp);

    
    Vfloat Elev;
    
    LL.Find_pipi_energy_lev(L,my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa, Elev);

    //Find energy levs corresponding to 2L;

    Vfloat Elev_2L;

    LL.Find_pipi_energy_lev(1.5*L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa, Elev_2L);


    Vfloat Elev_MpiL_4dot2;

    double Lprime= 4.2/(!Use_Mpi_OS?Mpi_fit.ave(i_ens):Mpi_OS_fit.ave(i_ens));

    LL.Find_pipi_energy_lev(Lprime, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa, Elev_MpiL_4dot2);

    



    auto F_int= [&](double t) {


		  
		  //double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
		  double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
		  //double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  double func_val =  Qfact*LL.V_pipi(t, L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa, Elev)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

		  double F_int_val = kern_val*func_val;
        

		  return F_int_val;
		};


    auto F_int_2L = [&](double t) {

		      //double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
		      double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
		      //double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		      double func_val =  Qfact*LL.V_pipi(t, 1.5*L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa, Elev_2L)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

		      double F_int_val = kern_val*func_val;
        

		      return F_int_val;
		};

    auto F_int_Lprime = [&](double t) {

			  //  double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
		      double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
		      // double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		      double func_val =  Qfact*LL.V_pipi(t, Lprime, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa,  Elev_MpiL_4dot2)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

		      double F_int_val = kern_val*func_val;
        

		      return F_int_val;
		};

    auto F_infL = [&](double t) {

		    //double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
		    double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
		    //double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);

		    double func_val =  Qfact*LL.V_pipi_infL(t, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

		    double F_int_val = kern_val*func_val;
        

		    return F_int_val;

		  };

   

    
    
    //Vfloat Sum_Series_Terms;
    double agm2_summ=0.0;
    double agm2_2L_summ=0.0;
    double agm2_infL_summ=0.0;
    double agm2_Lprime_summ=0.0;
    double Nterms = 2000;
    double Nterms_2L= 2000;
    
    //actual volume
    for(int iterm=1;iterm<Nterms;iterm++) {
      double F_n = F_int((double)iterm);
      //Sum_Series 
      agm2_summ += F_n;}
    agm2_dual.distr.push_back(agm2_summ);
   

    //1.5L volume
    for(int iterm=1;iterm<Nterms_2L;iterm++) {
      double F_n_2L = F_int_2L((double)iterm);
      //Sum Series 
      agm2_2L_summ += F_n_2L;}
    agm2_2L_dual.distr.push_back(agm2_2L_summ);


     //Mpi*L=4.2 volume
    for(int iterm=1;iterm<Nterms_2L;iterm++) {
      double F_n_Lprime = F_int_Lprime((double)iterm);
      //Sum Series 
      agm2_Lprime_summ += F_n_Lprime;}
    agm2_Lprime_dual.distr.push_back(agm2_Lprime_summ);
    

    //inf volume
    for(int iterm=1; iterm<Nterms_2L;iterm++) {
      double F_n_infL = F_infL((double)iterm);
      //Sum series
      agm2_infL_summ += F_n_infL;
    }
    agm2_infL_dual.distr.push_back(agm2_infL_summ);

      
    

    }

    a_distr_list.distr_list.push_back(a_distr);
    

    //push_back
    Edual.distr_list.push_back(Ed);
    Rdual.distr_list.push_back(Rd);
    Mrho.distr_list.push_back(Mr);
    gpi.distr_list.push_back(g);
    Kappa.distr_list.push_back(kap);
    agm2_light_fit.distr_list.push_back(agm2_dual);
    agm2_light_2L_fit.distr_list.push_back(agm2_2L_dual);
    agm2_light_infL_fit.distr_list.push_back(agm2_infL_dual);
    agm2_light_Lprime_fit.distr_list.push_back(agm2_Lprime_dual);
    cout<<"Edual/Mpi: "<<(Ed).ave()<<" +- "<<(Ed).err()<<endl;
    cout<<"Edual: "<<(Ed*mp).ave()<<" +- "<<(Ed*mp).err()<<endl;
    cout<<"Rdual: "<<Rd.ave()<<" +- "<<Rd.err()<<endl;
    cout<<"Mrho/Mpi: "<<(Mr).ave()<<" +- "<<(Mr).err()<<endl;
    cout<<"Mrho (lattice units): "<<(Mr*mp).ave()<<" +- "<<(Mr*mp).err()<<endl;
    cout<<"Mrho (GeV): "<<(Mr*mp/a_distr).ave()<<" +- "<<(Mr*mp/a_distr).err()<<endl;
    cout<<"Mrho/fp: "<<(Mr*mp/fp_fit[i_ens]).ave()<<" +- "<<(Mr*mp/fp_fit[i_ens]).err()<<endl;
    cout<<"gpi: "<<g.ave()<<" +- "<<g.err()<<endl;
    cout<<"kappa: "<<kap.ave()<<" +- "<<kap.err()<<endl;
    cout<<"Mpi: "<<mp.ave()<<" +- "<<mp.err()<<endl;
    cout<<"Mpi*L : "<<mp.ave()*L_list[i_ens]<<" +- "<<mp.err()*L_list[i_ens]<<endl;
    cout<<"L: "<<L_list[i_ens]<<endl;
    cout<<"agm2 pp+dual (L): "<<agm2_dual.ave()<<" +- "<<agm2_dual.err()<<endl;
    cout<<"agm2 pp+dual (1.5*L): "<<agm2_2L_dual.ave()<<" +- "<<agm2_2L_dual.err()<<endl;
    cout<<"agm2 pp+dual (Mpi*L=4.2): "<<agm2_Lprime_dual.ave()<<" +- "<<agm2_Lprime_dual.err()<<endl;
    cout<<"agm2 pp+dual (L -> infty): "<<agm2_infL_dual.ave()<<" +- "<<agm2_infL_dual.err()<<endl; 
    cout<<"#######################################"<<endl;
    
    



    }


  
  cout<<"########### DONE #############"<<endl;
  




  //strange
  Print_To_File(V_strange_1.Tag, {L_strange_list, a_strange_list, ml_strange_list, ZV_fit_strange.ave(), ZV_fit_strange.err(), MV_fit_strange.ave(), MV_fit_strange.err(), agm2_strange.ave(), agm2_strange.err()}, "../data/gm2/strange/agm2_fit.list", "", "# Ens L a ml  ZV   MV    agm2");

     

  //charm
  Print_To_File(V_charm_1.Tag, {L_charm_list, a_charm_list, ml_charm_list, ZV_fit_charm.ave(), ZV_fit_charm.err(), MV_fit_charm.ave(), MV_fit_charm.err(), agm2_charm.ave(), agm2_charm.err()}, "../data/gm2/charm/agm2_fit.list", "", "# Ens  L a ml ZV MV  agm2");

 
  //light
  Print_To_File(V_light_1.Tag, {ZV_fit_light.ave(), ZV_fit_light.err(), MV_fit_light.ave(), MV_fit_light.err()}, "../data/gm2/light/"+light_suff+"/ZV_MV_fitted.list", "", "#Ens ZV MV");

 
  //print fitted pars for all ensembles
  Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), Edual.ave(), Edual.err(), Rdual.ave(), Rdual.err(), Mrho.ave(), Mrho.err(), (Mrho*Mpi_fit/a_distr_list).ave(), (Mrho*Mpi_fit/a_distr_list).err(), gpi.ave(), gpi.err()}, "../data/gm2/light/"+light_suff+"/fit_pars.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs  Edual Rdual Mrho g");

 


  
  //tm and OS pion mass, decay constant, RCs and agm2
  Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), agm2_light_fit.ave(), agm2_light_fit.err(), agm2_light_2L_fit.ave(), agm2_light_2L_fit.err(), agm2_light_Lprime_fit.ave(), agm2_light_Lprime_fit.err()   , agm2_light_infL_fit.ave(), agm2_light_infL_fit.err()}, "../data/gm2/light/"+light_suff+"/agm2_fit.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs    agm2(L)    agm2(1.5L) agm2(Mpi*L=4.2)   agm2(infL)");

 
   



  return; 

}
