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





void Gm2() {


  data_t  V_light_1, V_light_2, V_light_3;
  data_t  V_strange_1, V_strange_2, V_strange_3;
  data_t   V_charm_1, V_charm_2, V_charm_3;

  //Read data
 
  V_light_1.Read("../gm2_data/light", "mes_contr_2pts_ll", "V1V1");
  V_light_2.Read("../gm2_data/light", "mes_contr_2pts_ll", "V2V2");
  V_light_3.Read("../gm2_data/light", "mes_contr_2pts_ll", "V3V3");
  V_strange_1.Read("../gm2_data/strange", "mes_contr_2pts_sl", "V1V1");
  V_strange_2.Read("../gm2_data/strange", "mes_contr_2pts_sl", "V2V2");
  V_strange_3.Read("../gm2_data/strange", "mes_contr_2pts_sl", "V3V3");
  V_charm_1.Read("../gm2_data/charm", "mes_contr_2pts_cl", "V1V1");
  V_charm_2.Read("../gm2_data/charm", "mes_contr_2pts_cl", "V2V2");
  V_charm_3.Read("../gm2_data/charm", "mes_contr_2pts_cl", "V3V3");
  

  int Nens_light= V_light_1.size;
  int Nens_strange = V_strange_1.size;
  int Nens_charm = V_charm_1.size;
  cout<<"N_ens light: "<<Nens_light<<endl;
  cout<<"N_ens strange: "<<Nens_strange<<endl;
  cout<<"N_ens charm: "<<Nens_charm<<endl;


  //create directories
  boost::filesystem::create_directory("../gm2");
  boost::filesystem::create_directory("../gm2/light");
  boost::filesystem::create_directory("../gm2/strange");
  boost::filesystem::create_directory("../gm2/charm");

  //define distr_t_list to be used in chiral+continuum analysis

  distr_t_list agm2_strange, agm2_charm;
  
  
  //strange 
  for(int i_ens=0;i_ens<Nens_strange;i_ens++) {
    
  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = V_strange_1.nrows[i_ens];
 
  distr_t_list  V_strange_1_distr, V_strange_2_distr, V_strange_3_distr;
  distr_t_list  V_strange_distr, MV_strange_distr, ZV_strange_distr;
  distr_t MV_strange, ZV_strange;

   //Analyze correlators

  Corr.Tmin=12;
  Corr.Tmax=20;
  
 
  //strange
  V_strange_1_distr = Corr.corr_t(V_strange_1.col(0)[i_ens], "../data/gm2/strange/corr_1_"+V_strange_1.Tag[i_ens]+".dat");
  V_strange_2_distr = Corr.corr_t(V_strange_2.col(0)[i_ens], "../data/gm2/strange/corr_2_"+V_strange_2.Tag[i_ens]+".dat");
  V_strange_3_distr = Corr.corr_t(V_strange_3.col(0)[i_ens], "../data/gm2/strange/corr_3_"+V_strange_3.Tag[i_ens]+".dat");

 
  //sum over the Lorenz indices of the e.m. current
  V_strange_distr= (pow(qs,2)/3.0)*(V_strange_1_distr+ V_strange_2_distr + V_strange_3_distr);
 
  // print summed correlators to file
  Print_To_File({}, {V_strange_distr.ave(), V_strange_distr.err()}, "../data/gm2/strange/corr_sum_"+V_strange_1.Tag[i_ens]+".dat", "", "");
  
  //extract effective masses, overlap from V and fit
 
  //strange
  MV_strange_distr= Corr.effective_mass_t(V_strange_distr, "../data/gm2/strange/MV_mass_"+V_strange_1.Tag[i_ens]+".dat");
  ZV_strange_distr= Corr.residue_t(V_strange_distr, "../data/gm2/strange/ZV_overlap_"+V_strange_1.Tag[i_ens]+".dat");
  MV_strange = Corr.Fit_distr(MV_strange_distr);
  ZV_strange = Corr.Fit_distr(ZV_strange_distr);
  }


  //charm 
  for(int i_ens=0;i_ens<Nens_charm;i_ens++) {
    
  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = V_strange_1.nrows[i_ens];
  distr_t_list  V_charm_1_distr, V_charm_2_distr, V_charm_3_distr;
  distr_t_list  V_charm_distr, MV_charm_distr, ZV_charm_distr;
  distr_t MV_charm , ZV_charm;

  //Analyze correlators

  Corr.Tmin=12;
  Corr.Tmax=20;
  
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
  MV_charm_distr= Corr.effective_mass_t(V_charm_distr, "../data/gm2/charm/MV_mass_"+V_charm_1.Tag[i_ens]+".dat");
  ZV_charm_distr= Corr.residue_t(V_charm_distr, "../data/gm2/charm/ZV_overlap_"+V_charm_1.Tag[i_ens]+".dat");
  MV_charm = Corr.Fit_distr(MV_charm_distr);
  ZV_charm = Corr.Fit_distr(ZV_charm_distr);
  }



  //light
  for(int i_ens=0;i_ens<Nens_light;i_ens++) {
    
  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = V_strange_1.nrows[i_ens];
  distr_t_list  V_light_1_distr, V_light_2_distr, V_light_3_distr;
  distr_t_list  V_light_distr, MV_light_distr, ZV_light_distr;
  distr_t MV_light, ZV_light;

   //Analyze correlators

  Corr.Tmin=12;
  Corr.Tmax=20;
  
  //light
  V_light_1_distr = Corr.corr_t(V_light_1.col(0)[i_ens], "../data/gm2/light/corr_1_"+V_light_1.Tag[i_ens]+".dat");
  V_light_2_distr = Corr.corr_t(V_light_2.col(0)[i_ens], "../data/gm2/light/corr_2_"+V_light_2.Tag[i_ens]+".dat");
  V_light_3_distr = Corr.corr_t(V_light_3.col(0)[i_ens], "../data/gm2/light/corr_3_"+V_light_3.Tag[i_ens]+".dat");

 
  //sum over the Lorenz indices of the e.m. current
  V_light_distr= ((pow(qu,2)+pow(qd,2))/3.0)*(V_light_1_distr+ V_light_2_distr + V_light_3_distr);
 
  // print summed correlators to file
  Print_To_File({}, {V_light_distr.ave(), V_light_distr.err()}, "../data/gm2/light/corr_sum_"+V_light_1.Tag[i_ens]+".dat", "", "");

  //extract effective masses, overlap from V and fit

  //light
  MV_light_distr= Corr.effective_mass_t(V_light_distr, "../data/gm2/light/MV_mass_"+V_light_1.Tag[i_ens]+".dat");
  ZV_light_distr= Corr.residue_t(V_light_distr, "../data/gm2/light/ZV_overlap_"+V_light_1.Tag[i_ens]+".dat");
  MV_light = Corr.Fit_distr(MV_light_distr);
  ZV_light = Corr.Fit_distr(ZV_light_distr);

  }







  




  return; 

}
