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
const bool Use_JB_distribution= false;
const bool UseJack=1;
const int Njacks=15;
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

class ipar {

public:
  ipar() : V_light(0.0), V_light_err(0.0) {}

  
  double fp, Mp,  Mv;
  double t,  L;
  double V_light, V_light_err;
   


};

class fit_par {

public:
  fit_par() {}
  fit_par(const Vfloat &par) {
    if((signed)par.size() != 5) crash("In class fit_par in fitting analytic representation of V(t)_light, class constructor Vfloat par has size != 5");
    Rd=par[0];
    Ed=par[1];
    Mrho=par[2];
    gpi=par[3];
    Za=par[4];
  }

  double Rd,Ed, Mrho, gpi,Za;



};



void Gm2() {


  //init Gaussian number generator
  GaussianMersenne GM(981832);
  string channel="";

  
  data_t  V_light_1, V_light_2, V_light_3, pt2_pion;
  data_t pt2_pion_strange, pt2_pion_charm;
  data_t  V_strange_1, V_strange_2, V_strange_3;
  data_t  V_charm_1, V_charm_2, V_charm_3;

  //Read data
 
  V_light_1.Read("../gm2_data/light", "mes_contr_2pts_ll_ave", "V1V1");
  V_light_2.Read("../gm2_data/light", "mes_contr_2pts_ll_ave", "V2V2");
  V_light_3.Read("../gm2_data/light", "mes_contr_2pts_ll_ave", "V3V3");
  pt2_pion.Read("../gm2_data/light", "mes_contr_2pts_ll_ave", "P5P5");
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
  boost::filesystem::create_directory("../data/gm2/strange");
  boost::filesystem::create_directory("../data/gm2/charm");

  //define distr_t_list to be used in chiral+continuum analysis

  distr_t_list agm2_strange, agm2_charm, agm2_light;
  vector<vector<fit_par>> par_list_anal_repr;//only used for light quark contribution

  distr_t_list MV_fit_light(UseJack), MV_fit_strange(UseJack), MV_fit_charm(UseJack);
  distr_t_list ZV_fit_light(UseJack), ZV_fit_strange(UseJack), ZV_fit_charm(UseJack);
  distr_t_list Mpi_fit(UseJack), fp_fit(UseJack);
  distr_t_list Mpi_dim;

  //define lambda for convolution with kernel
  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv, channel);};
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
  Print_To_File({}, {V_strange_distr.ave(), V_strange_distr.err()}, "../data/gm2/strange/corr_sum_"+V_strange_1.Tag[i_ens]+".dat", "", "");
  
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
  distr_t_list Kernel_distr = distr_t_list::f_of_distr(K,MV_strange, Upper_Limit_Time_Integral_strange+1);
  //compute exp(-Mv*t) distribution
  distr_t_list exp_MVs = distr_t_list::f_of_distr(exp_MV, MV_strange, Upper_Limit_Time_Integral_strange+1);

  //Print single-exponential prediction to file
  Print_To_File({}, {(exp_MVs*(ZV_strange/(2.0*MV_strange))).ave(), (exp_MVs*(ZV_strange/(2.0*MV_strange))).err()}, "../data/gm2/strange/corr_gsd_sum_"+V_strange_1.Tag[i_ens]+".dat", "", "");

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
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata.ave(), agm2_distr_Tdata.err()}, "../data/gm2/strange/agm2_Tdata_"+V_strange_1.Tag[i_ens]+".dat", "", "#id  Tdata   ag2m agm2_err");
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
  Print_To_File({}, {V_charm_distr.ave(), V_charm_distr.err()}, "../data/gm2/charm/corr_sum_"+V_charm_1.Tag[i_ens]+".dat", "", "");


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
  distr_t_list Kernel_distr = distr_t_list::f_of_distr(K,MV_charm, Upper_Limit_Time_Integral_strange+1);
  //compute exp(-Mv*t) distribution
  distr_t_list exp_MVc = distr_t_list::f_of_distr(exp_MV, MV_charm, Upper_Limit_Time_Integral_charm+1);

  //Print single-exponential prediction to file
  Print_To_File({}, {(exp_MVc*(ZV_charm/(2.0*MV_charm))).ave(), (exp_MVc*(ZV_charm/(2.0*MV_charm))).err()}, "../data/gm2/charm/corr_gsd_sum_"+V_charm_1.Tag[i_ens]+".dat", "", "");

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
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata.ave(), agm2_distr_Tdata.err()}, "../data/gm2/charm/agm2_Tdata_"+V_charm_1.Tag[i_ens]+".dat", "", "#id  Tdata   ag2m agm2_err");
  }


  cout<<"charm quark correlator analyzed!"<<endl;

  //light
  channel="l";
  for(int i_ens=0;i_ens<Nens_light;i_ens++) {
    
  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = V_light_1.nrows[i_ens];

  //resample lattice spacing
  distr_t a_distr(UseJack), ZA_distr(UseJack);
  LatticeInfo L_info;
  L_info.LatInfo_new_ens(V_light_1.Tag[i_ens]);
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



  
  distr_t_list  V_light_1_distr, V_light_2_distr, V_light_3_distr;
  distr_t_list  V_light_distr, MV_light_distr, ZV_light_distr;
  distr_t_list Mpi_distr, fp_distr;
  distr_t MV_light, ZV_light, Mpi, fp;
  

  //Analyze correlators
  if(Corr.Nt > 48) {
  Corr.Tmin=15;
  Corr.Tmax=30;
  }
  else { Corr.Tmin= 12; Corr.Tmax=20;}
  
  //light
  V_light_1_distr = Corr.corr_t(V_light_1.col(0)[i_ens], "../data/gm2/light/corr_1_"+V_light_1.Tag[i_ens]+".dat");
  V_light_2_distr = Corr.corr_t(V_light_2.col(0)[i_ens], "../data/gm2/light/corr_2_"+V_light_2.Tag[i_ens]+".dat");
  V_light_3_distr = Corr.corr_t(V_light_3.col(0)[i_ens], "../data/gm2/light/corr_3_"+V_light_3.Tag[i_ens]+".dat");
  Mpi_distr = Corr.effective_mass_t( pt2_pion.col(0)[i_ens], "../data/gm2/light/Mpi_"+pt2_pion.Tag[i_ens]+".dat");
  fp_distr = Corr.decay_constant_t(pt2_pion.col(0)[i_ens], "../data/gm2/light/fp_"+pt2_pion.Tag[i_ens]+".dat");


 
  //sum over the Lorenz indices of the e.m. current
  V_light_distr= ((pow(qu,2)+pow(qd,2))/3.0)*(V_light_1_distr+ V_light_2_distr + V_light_3_distr);
 
  // print summed correlators to file
  Print_To_File({}, {V_light_distr.ave(), V_light_distr.err()}, "../data/gm2/light/corr_sum_"+V_light_1.Tag[i_ens]+".dat", "", "");

  //extract effective masses, overlap from V and fit

  //light
  MV_light_distr= Corr.effective_mass_t(V_light_distr, "../data/gm2/light/MV_mass_"+V_light_1.Tag[i_ens]+".dat");
  ZV_light_distr= Corr.residue_t(V_light_distr, "../data/gm2/light/ZV_overlap_"+V_light_1.Tag[i_ens]+".dat");
  int Tmin_old=Corr.Tmin;
  int Tmax_old=Corr.Tmax;
  if(V_light_1.Tag[i_ens].substr(1,1)=="C") {Corr.Tmin=14;Corr.Tmax=19;}
  else if(V_light_1.Tag[i_ens].substr(1,1)=="B" || V_light_1.Tag[i_ens].substr(1,1) =="A") {Corr.Tmin=12; Corr.Tmax=18;}
  else crash("Cannot find ensemble: "+V_light_1.Tag[i_ens]+" while modifying fit interval for MV_light");
  MV_light = Corr.Fit_distr(MV_light_distr);
  ZV_light = Corr.Fit_distr(ZV_light_distr);
  //push_back MV_light and ZV_light
  MV_fit_light.distr_list.push_back(MV_light);
  ZV_fit_light.distr_list.push_back(ZV_light);
  Corr.Tmin=Tmin_old;
  Corr.Tmax=Tmax_old;
  Mpi=Corr.Fit_distr(Mpi_distr);
  fp= Corr.Fit_distr(fp_distr);
  //push_back Mpi and fp
  Mpi_fit.distr_list.push_back(Mpi);
  Mpi_dim.distr_list.push_back(Mpi/a_distr);
  fp_fit.distr_list.push_back(fp);

  int Tdata_min= 8;
  int Tdata_max = Corr.Nt/2.0 -4;
  int Tdata_fit = (Corr.Tmax+Corr.Tmin)/2;
  
  //compute kernel distribution
  distr_t_list Kernel_distr = distr_t_list::f_of_distr(K,a_distr, Upper_Limit_Time_Integral_light+1);
  //compute exp(-Mv*t) distribution
  distr_t_list exp_MVl = distr_t_list::f_of_distr(exp_MV, MV_light, Upper_Limit_Time_Integral_light+1);

  //Print single-exponential prediction to file
  Print_To_File({}, {(exp_MVl*(ZV_light/(2.0*MV_light))).ave(), (exp_MVl*(ZV_light/(2.0*MV_light))).err()}, "../data/gm2/light/corr_gsd_sum_"+V_light_1.Tag[i_ens]+".dat", "", "");

  distr_t_list agm2_distr_Tdata(UseJack);
  Vfloat Tdata_vec;
  
  for(int Tdata=Tdata_min;Tdata<Tdata_max;Tdata++) {
    //compute 4\pia^2 using lattice data up to Tdata (included)
   
    distr_t agm2(UseJack, UseJack?Njacks:Nboots);
    
    for(int t=1;t<=Upper_Limit_Time_Integral_light;t++) {
      if(t<=Tdata) agm2 = agm2 + 4.0*pow(alpha,2)*V_light_distr.distr_list[t]*Kernel_distr.distr_list[t];
      else agm2= agm2 + 4.0*pow(alpha,2)*(ZV_light/(2.0*MV_light))*exp_MVl.distr_list[t]*Kernel_distr.distr_list[t];
    }
    Tdata_vec.push_back((double)Tdata);
    agm2_distr_Tdata.distr_list.push_back(ZA_distr*ZA_distr*agm2);
    if(Tdata==Tdata_fit) agm2_light.distr_list.push_back(ZA_distr*ZA_distr*agm2);
  }
  //print to file
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata.ave(), agm2_distr_Tdata.err()}, "../data/gm2/light/agm2_Tdata_"+V_light_1.Tag[i_ens]+".dat", "", "#id  Tdata   ag2m agm2_err");



  //fit lattice data using analytic representation for V(t)_light


  bootstrap_fit<fit_par,ipar> bf(Njacks);

  int Tfit_points=Corr.Nt/2 -5;
  int Tfit_min=3;

  bf.Set_number_of_measurements(Tfit_points);
  bf.Set_verbosity(verbosity);

  //initial guesses are based on the results on the old ETMC confs
  bf.Add_par("Rd", 1.3, 0.03);
  bf.Add_par("Ed", 2.0, 0.03);
  bf.Add_par("Mrho",2.6, 0.04);
  bf.Add_par("gpi", 5.0, 0.05);
  bf.Add_prior_par("Za", ZA_distr.ave(), ZA_distr.err());
  bf.Fix_n_release("Za");

  bf.ansatz =  [=](const fit_par &p, const ipar &ip) -> double {
		 return V_pipi(ip.t, ip.L, p.Mrho*ip.Mp, p.gpi, ip.Mp) + Vdual(ip.t, p.Mrho*ip.Mp, p.Ed*ip.Mp, p.Rd);
  };
  bf.measurement = [=](const fit_par& p,const ipar& ip) -> double {
		     return pow(p.Za,2)*ip.V_light;
  };
  bf.error =  [=](const fit_par& p,const ipar &ip) -> double {
		return pow(p.Za,2)*ip.V_light_err;
  };

  //fill the data
  vector<vector<ipar>> data(Njacks);
  //allocate space for output result
  boot_fit_data<fit_par> Bt_fit;

  

  for(auto &data_iboot: data) data_iboot.resize(Tfit_points);
  
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int t=Tfit_min;t<Tfit_points+Tfit_min;t++) {
	int tt=t-Tfit_min;
	data[ijack][tt].V_light = V_light_distr.distr_list[t].distr[ijack];
	data[ijack][tt].V_light_err = V_light_distr.distr_list[t].err();
	data[ijack][tt].L = L_info.L;
	data[ijack][tt].Mp = Mpi.distr[ijack];
	data[ijack][tt].fp = fp.distr[ijack];
	data[ijack][tt].Mv = MV_light.distr[ijack];
	
      } 
    }

    

    //add prior values
    for(int ijack=0;ijack<Njacks;ijack++) {
     bf.ib= &ijack;
     bf.Append_to_prior("Za", ZA_distr.distr[ijack], ZA_distr.err());
    }

    
    //append
    bf.Append_to_input_par(data);

    
    
    //fit
    Bt_fit= bf.Perform_bootstrap_fit();

    //###################################################
    //print fitted func

    int nsteps=1000; //500 points in time
    int step_size_time= 0.2;
    distr_t_list func(UseJack, nsteps);
    Vfloat times;

    //generate data
    for(int istep=0; istep<nsteps;istep++) {
      for(int ijack=0;ijack<Njacks;ijack++) {
	ipar my_ipar;
	my_ipar.L= L_info.L;
	my_ipar.Mp=Mpi.distr[ijack];
	my_ipar.fp=fp.distr[ijack];
	my_ipar.Mv = MV_light.distr[ijack];
	my_ipar.t = 0.2 +istep*step_size_time;
	times.push_back(my_ipar.t);
	func.distr_list[istep].distr.push_back( bf.ansatz( Bt_fit.par[ijack], my_ipar ));
      }
    }

    //print to file
    Print_To_File({}, {times, func.ave(), func.err()}, "../data/gm2/light/Corr_anal_rep_"+V_light_1.Tag[i_ens]+".dat", "", "#id_col  time  val   err  (fit interval: ["+to_string(Tfit_min)+","+to_string(Tfit_points+Tfit_min-1)+"])");
    //####################################################
    
    
    //push_back jack distribution of fitted params
    par_list_anal_repr.push_back(Bt_fit.par);

   
  }

   cout<<"light quark correlator analyzed!"<<endl;


   distr_t_list Mpi_dim_strange(UseJack), Mpi_dim_charm(UseJack);

   for(int iens=0;iens<Nens_strange;iens++) Mpi_dim_strange.distr_list.push_back( Mpi_dim.distr_list[V_light_1.Get_iens_from_tag(V_strange_1.Tag[iens])]);
   for(int iens=0;iens<Nens_charm;iens++) Mpi_dim_charm.distr_list.push_back( Mpi_dim.distr_list[V_light_1.Get_iens_from_tag(V_charm_1.Tag[iens])]);

   //print ZV_fit, MV_fit , Mpi_fit, fp_fit, agm2 (l,s,c) 

   //strange
   Print_To_File(V_strange_1.Tag, {ZV_fit_strange.ave(), ZV_fit_strange.err()}, "../data/gm2/strange/ZV_fitted.list", "", "");
   Print_To_File(V_strange_1.Tag, {MV_fit_strange.ave(), MV_fit_strange.err()}, "../data/gm2/strange/MV_fitted.list", "", "");
   Print_To_File(V_strange_1.Tag, {Mpi_dim_strange.ave(), Mpi_dim_strange.err(), agm2_strange.ave(), agm2_strange.err()}, "../data/gm2/strange/agm2.list", "", "# Ens   Mpi   Mpi_err    agm2    agm2_err");
   
   

   //charm
   Print_To_File(V_charm_1.Tag, {ZV_fit_charm.ave(), ZV_fit_charm.err()}, "../data/gm2/charm/ZV_fitted.list", "", "");
   Print_To_File(V_charm_1.Tag, {MV_fit_charm.ave(), MV_fit_charm.err()}, "../data/gm2/charm/MV_fitted.list", "", "");
   Print_To_File(V_charm_1.Tag, {Mpi_dim_charm.ave(), Mpi_dim_charm.err(), agm2_charm.ave(), agm2_charm.err()}, "../data/gm2/charm/agm2.list", "", "# Ens   Mpi   Mpi_err    agm2    agm2_err");

   //light
   Print_To_File(V_light_1.Tag, {ZV_fit_light.ave(), ZV_fit_light.err()}, "../data/gm2/light/ZV_fitted.list", "", "");
   Print_To_File(V_light_1.Tag, {MV_fit_light.ave(), MV_fit_light.err()}, "../data/gm2/light/MV_fitted.list", "", "");
   Print_To_File(V_light_1.Tag, {Mpi_dim.ave(), Mpi_dim.err(), agm2_light.ave(), agm2_light.err()}, "../data/gm2/light/agm2.list", "", "# Ens   Mpi   Mpi_err    agm2    agm2_err");
   Print_To_File(V_light_1.Tag, {Mpi_fit.ave(), Mpi_fit.err()}, "../data/gm2/light/Mpi_fitted.list", "", "");
   Print_To_File(V_light_1.Tag, {fp_fit.ave(), fp_fit.err()}, "../data/gm2/light/fp_fitted.list", "", "");
   



  return; 

}
