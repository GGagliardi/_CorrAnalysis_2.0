#include "../include/gm2.h"


using namespace std;

//set constants
namespace plt = matplotlibcpp;

constexpr double kappa=2.837297;
const double MPiPhys=0.135;
const double alpha = 1.0/137.04;
const double e2 = alpha*4.0*M_PI;
const int nboots= 200;
const bool UseJack=1;
const int Njacks=50;
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
const int Luscher_num_zeroes= 20;
const int npts_spline= 1000;
bool Use_Mpi_OS=false;
bool Include_light_disco= true;
double Qfact= 10.0/9.0;
const double m_Jpsi= 3.0969;
const double m_phi= 1.019461;
const double m_phi_err= 0.000016;
const double m_pi = 0.135;
const double m_rho= 0.775;
const double m_k= 0.497611;
const double m_d= 1.86484;
const double m_etac= 2.9839;
const double m_etas = 0.6859;
const double m_etas_err= 0.0040;
const double fp_phys= 0.1304;
const double fk_phys= 0.1557; //FLAG REV ("2021")
const double fd_phys= 0.212; //FLAG REV ("2021")
const double csi_phys= pow(m_pi/(4.0*M_PI*fp_phys),2);
const double t0 = 0.4*fm_to_inv_Gev;
const double t1 = 1.0*fm_to_inv_Gev;
const double Delta= 0.15*fm_to_inv_Gev;
const double l1ph= -0.4; //-0.4
const double l2ph= 4.3; //4.3
const double l3ph= 3.2; //3.2
const double l4ph= 4.4; //4.4
const double s0= 2.0-M_PI/2.0;
const double s1 = M_PI/4.0 - 0.5;
const double s2 = 0.5 - M_PI/8.0;
const double s3 = 3.0*M_PI/16.0 - 0.5;
const double Csi_K_phys= m_k*m_k/(fk_phys*fk_phys) ;
const double Csi_D_phys= m_d*m_d/(fd_phys*fd_phys) ;
bool Skip_total_light_calc= true;
string Mod="FLEE_THEORY";
bool Use_Za_Zv_from_strange_run=true;
bool Use_Za_Zv_from_charm_run = true;
bool Use_phi_etas = true;

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
    if((signed)par.size() != 5) crash("In class fit_par in fitting analytic representation of V(t)_light, class constructor Vfloat par has size != 5");
    Rd=par[0];
    Ed=par[1];
    Mrho=par[2];
    gpi=par[3];
    kappa= par[4];
    
  }

  double Rd,Ed, Mrho, gpi,kappa;



};




void Gm2() {

  omp_set_num_threads(1);

  if(Mod=="FREE_THEORY") { Compute_SD_window_Free(); exit(-1);}


  //Compute_clogSD_free(20);

  //exit(-1);

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

   
 
 
  //init Gaussian number generator
  GaussianMersenne GM(981832);
  string channel="";

  
  data_t  V_light_1, V_light_2, V_light_3, pt2_pion;
  data_t  V_light_OS_1, V_light_OS_2, V_light_OS_3;
  data_t corr_P5P5_strange, corr_P5P5_strange_heavy,  pt2_pion_charm, corr_P5P5_OS_strange, corr_P5P5_OS_strange_heavy,  pt2_pion_OS_charm;
 

  //L
  data_t  V_strange_1_L, V_strange_2_L, V_strange_3_L, V_strange_OS_1_L, V_strange_OS_2_L, V_strange_OS_3_L;
  data_t  V_charm_1_L, V_charm_2_L, V_charm_3_L, V_charm_OS_1_L, V_charm_OS_2_L, V_charm_OS_3_L;
  data_t pt2_K_L, pt2_D_L, pt2_etaC_L, pt2_etaC_OS_L;
  //M
  data_t  V_strange_1_M, V_strange_2_M, V_strange_3_M, V_strange_OS_1_M, V_strange_OS_2_M, V_strange_OS_3_M;
  data_t  V_charm_1_M, V_charm_2_M, V_charm_3_M, V_charm_OS_1_M, V_charm_OS_2_M, V_charm_OS_3_M;
  data_t pt2_K_M, pt2_D_M, pt2_etaC_M, pt2_etaC_OS_M;
  //H
  data_t  V_strange_1_H, V_strange_2_H, V_strange_3_H, V_strange_OS_1_H, V_strange_OS_2_H, V_strange_OS_3_H;
  data_t  V_charm_1_H, V_charm_2_H, V_charm_3_H, V_charm_OS_1_H, V_charm_OS_2_H, V_charm_OS_3_H;
  data_t pt2_K_H, pt2_D_H, pt2_etaC_H, pt2_etaC_OS_H;
  //to compute ZV 
  data_t corr_A0P5;
  //to compute ZA
  data_t corr_A0P5_OS, corr_P5P5_OS;
  //disco
  data_t disco_light;

  //to compute ZV (strange)
  data_t corr_A0P5_strange, corr_A0P5_strange_heavy;
  //to compute ZA (strange)
  data_t corr_A0P5_OS_strange, corr_A0P5_OS_strange_heavy;


  
  //to compute ZV and ZA (charm)
  data_t corr_A0P5_charm_L , corr_A0P5_OS_charm_L;
  data_t corr_A0P5_charm_M , corr_A0P5_OS_charm_M;
  data_t corr_A0P5_charm_H , corr_A0P5_OS_charm_H;

  
  //Read data

  //Custom sorting for V_light to account for the two replica r0 and r1
  auto Sort_light_confs = [](string A, string B) {
		      
		      string rA = A.substr(A.length()-2);
		      string rB = B.substr(B.length()-2);
		      if(rA.substr(0,1) == "r") { 
		      int n1 = stoi(A.substr(A.length()-1));
		      int n2 = stoi(B.substr(B.length() -1));
		      if(rA==rB) return A<B;
		      else return n1<n2;
		      }
		      
		      return A<B;
		    };

  //#################################END CUSTOM SORTING#################
  V_light_1.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs);
  V_light_2.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "V2V2", Sort_light_confs);
  V_light_3.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "V3V3", Sort_light_confs);

  V_light_OS_1.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs);
  V_light_OS_2.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "V2V2", Sort_light_confs);
  V_light_OS_3.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "V3V3", Sort_light_confs);
  
  pt2_pion.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  corr_A0P5.Read("../gm2_data/light", "mes_contr_2pts_ll_1", "P5A0", Sort_light_confs);
  corr_A0P5_OS.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "P5A0", Sort_light_confs);
  corr_P5P5_OS.Read("../gm2_data/light", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs);
  if(Include_light_disco) {
    disco_light.Read("../gm2_data/disco_light/data", "disco", "", Sort_light_confs);
  }
  
 

  //strange
  //L
  V_strange_1_L.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs); //17-18 for lighter ms  //21-22 for medium  ms,     //25-26 for heavier ms
  V_strange_2_L.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_1", "V2V2", Sort_light_confs);
  V_strange_3_L.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_1", "V3V3", Sort_light_confs);
  V_strange_OS_1_L.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs); 
  V_strange_OS_2_L.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_2", "V2V2", Sort_light_confs);
  V_strange_OS_3_L.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_2", "V3V3", Sort_light_confs);
  //M
  V_strange_1_M.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs); 
  V_strange_2_M.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_1", "V2V2", Sort_light_confs);
  V_strange_3_M.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_1", "V3V3", Sort_light_confs);
  V_strange_OS_1_M.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_2", "V1V1", Sort_light_confs); 
  V_strange_OS_2_M.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_2", "V2V2", Sort_light_confs);
  V_strange_OS_3_M.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_2", "V3V3", Sort_light_confs);

  //H
  V_strange_1_H.Read("../gm2_data/strange", "mes_contr_2pts_sl_25", "V1V1");
  V_strange_2_H.Read("../gm2_data/strange", "mes_contr_2pts_sl_25", "V2V2");
  V_strange_3_H.Read("../gm2_data/strange", "mes_contr_2pts_sl_25", "V3V3");
  V_strange_OS_1_H.Read("../gm2_data/strange", "mes_contr_2pts_sl_26", "V1V1"); 
  V_strange_OS_2_H.Read("../gm2_data/strange", "mes_contr_2pts_sl_26", "V2V2");
  V_strange_OS_3_H.Read("../gm2_data/strange", "mes_contr_2pts_sl_26", "V3V3");
  //P5P5
  corr_P5P5_strange.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  corr_P5P5_OS_strange.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs);
  corr_P5P5_strange_heavy.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  corr_P5P5_OS_strange_heavy.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs); 
  //A0P5
  corr_A0P5_strange.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_1", "P5A0", Sort_light_confs);
  corr_A0P5_OS_strange.Read("../gm2_data/strange_Nhits64/light", "mes_contr_2pts_ll_2", "P5A0", Sort_light_confs);
  corr_A0P5_strange_heavy.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_1", "P5A0", Sort_light_confs);
  corr_A0P5_OS_strange_heavy.Read("../gm2_data/strange_Nhits64/heavy", "mes_contr_2pts_ll_2", "P5A0", Sort_light_confs);
  //L
  pt2_K_L.Read("../gm2_data/strange", "mes_contr_2pts_sl_5", "P5P5");  //5 for lighter ms //9 for medium ms, 13 for heavier ms
  //M
  pt2_K_M.Read("../gm2_data/strange", "mes_contr_2pts_sl_9", "P5P5");
  //H
  pt2_K_H.Read("../gm2_data/strange", "mes_contr_2pts_sl_13", "P5P5");
  

  //charm
  //L
  V_charm_1_L.Read("../gm2_data/charm", "mes_contr_2pts_cl_17", "V1V1"); //17-18 for lighter mc  //21-22 for medium  mc,     //25-26 for heavier mc
  V_charm_2_L.Read("../gm2_data/charm", "mes_contr_2pts_cl_17", "V2V2");
  V_charm_3_L.Read("../gm2_data/charm", "mes_contr_2pts_cl_17", "V3V3");
  V_charm_OS_1_L.Read("../gm2_data/charm", "mes_contr_2pts_cl_18", "V1V1"); 
  V_charm_OS_2_L.Read("../gm2_data/charm", "mes_contr_2pts_cl_18", "V2V2");
  V_charm_OS_3_L.Read("../gm2_data/charm", "mes_contr_2pts_cl_18", "V3V3");
  //M
  V_charm_1_M.Read("../gm2_data/charm", "mes_contr_2pts_cl_21", "V1V1"); 
  V_charm_2_M.Read("../gm2_data/charm", "mes_contr_2pts_cl_21", "V2V2");
  V_charm_3_M.Read("../gm2_data/charm", "mes_contr_2pts_cl_21", "V3V3");
  V_charm_OS_1_M.Read("../gm2_data/charm", "mes_contr_2pts_cl_22", "V1V1"); 
  V_charm_OS_2_M.Read("../gm2_data/charm", "mes_contr_2pts_cl_22", "V2V2");
  V_charm_OS_3_M.Read("../gm2_data/charm", "mes_contr_2pts_cl_22", "V3V3");
  //H
  V_charm_1_H.Read("../gm2_data/charm", "mes_contr_2pts_cl_25", "V1V1");
  V_charm_2_H.Read("../gm2_data/charm", "mes_contr_2pts_cl_25", "V2V2");
  V_charm_3_H.Read("../gm2_data/charm", "mes_contr_2pts_cl_25", "V3V3");
  V_charm_OS_1_H.Read("../gm2_data/charm", "mes_contr_2pts_cl_26", "V1V1"); 
  V_charm_OS_2_H.Read("../gm2_data/charm", "mes_contr_2pts_cl_26", "V2V2");
  V_charm_OS_3_H.Read("../gm2_data/charm", "mes_contr_2pts_cl_26", "V3V3");

  

  
  //P5P5
  pt2_pion_charm.Read("../gm2_data/charm", "mes_contr_2pts_cl_1", "P5P5");
  pt2_pion_OS_charm.Read("../gm2_data/charm", "mes_contr_2pts_cl_2", "P5P5");
  //L
  pt2_D_L.Read("../gm2_data/charm", "mes_contr_2pts_cl_5", "P5P5");  //5 for lighter ms //9 for medium ms, 13 for heavier ms
  //M
  pt2_D_M.Read("../gm2_data/charm", "mes_contr_2pts_cl_9", "P5P5");
  //H
  pt2_D_H.Read("../gm2_data/charm", "mes_contr_2pts_cl_13", "P5P5");
  //eta_c
  //L
  pt2_etaC_L.Read("../gm2_data/charm", "mes_contr_2pts_cl_17", "P5P5");
  pt2_etaC_OS_L.Read("../gm2_data/charm", "mes_contr_2pts_cl_18", "P5P5");
  //M
  pt2_etaC_M.Read("../gm2_data/charm", "mes_contr_2pts_cl_21", "P5P5");
  pt2_etaC_OS_M.Read("../gm2_data/charm", "mes_contr_2pts_cl_22", "P5P5");
  //H
  pt2_etaC_H.Read("../gm2_data/charm", "mes_contr_2pts_cl_25", "P5P5");
  pt2_etaC_OS_H.Read("../gm2_data/charm", "mes_contr_2pts_cl_26", "P5P5");

  //A0P5 L
  corr_A0P5_charm_L.Read("../gm2_data/charm", "mes_contr_2pts_cl_17", "A0P5");
  corr_A0P5_OS_charm_L.Read("../gm2_data/charm", "mes_contr_2pts_cl_18", "A0P5");

  //A0P5 M
  corr_A0P5_charm_M.Read("../gm2_data/charm", "mes_contr_2pts_cl_21", "A0P5");
  corr_A0P5_OS_charm_M.Read("../gm2_data/charm", "mes_contr_2pts_cl_22", "A0P5");

  //A0P5 H
  corr_A0P5_charm_H.Read("../gm2_data/charm", "mes_contr_2pts_cl_25", "A0P5");
  corr_A0P5_OS_charm_H.Read("../gm2_data/charm", "mes_contr_2pts_cl_26", "A0P5");




  
  //create directories
  boost::filesystem::create_directory("../data/gm2");
  boost::filesystem::create_directory("../data/gm2/light");
  boost::filesystem::create_directory("../data/gm2/light/OS");
  boost::filesystem::create_directory("../data/gm2/light/tm");
  boost::filesystem::create_directory("../data/gm2/light/windows_fit_func");
  boost::filesystem::create_directory("../data/gm2/light/disco");
  boost::filesystem::create_directory("../data/gm2/strange");
  boost::filesystem::create_directory("../data/gm2/strange/tm");
  boost::filesystem::create_directory("../data/gm2/strange/OS");
  boost::filesystem::create_directory("../data/gm2/charm");
  boost::filesystem::create_directory("../data/gm2/charm/tm");
  boost::filesystem::create_directory("../data/gm2/charm/OS");


  /*

  //charm Nhits 20
  data_t V_charm_20_1_L, V_charm_20_1_M, V_charm_20_2_L, V_charm_20_2_M, V_charm_20_3_L, V_charm_20_3_M;
  data_t P5P5_charm_20_L, P5P5_charm_20_M;
  V_charm_20_1_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs);
  V_charm_20_2_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_1", "V2V2", Sort_light_confs);
  V_charm_20_3_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_1", "V3V3", Sort_light_confs);
  P5P5_charm_20_L.Read("../gm2_data/charm_Nhits20/light", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
  //M
  V_charm_20_1_M.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_1", "V1V1", Sort_light_confs); 
  V_charm_20_2_M.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_1", "V2V2", Sort_light_confs);
  V_charm_20_3_M.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_1", "V3V3", Sort_light_confs);
  P5P5_charm_20_M.Read("../gm2_data/charm_Nhits20/heavy", "mes_contr_2pts_ll_1", "P5P5", Sort_light_confs);
 
  int Nens_charm_20 = V_charm_20_1_L.size;

  for(int i_ens=0; i_ens < Nens_charm_20;i_ens++) {

    distr_t_list Vcharm1L, Vcharm2L, Vcharm3L, Vcharm1M, Vcharm2M, Vcharm3M;

    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    Corr.Nt = V_charm_20_1_L.nrows[i_ens];
    cout<<"Nt: "<<Corr.Nt<<endl;
    cout<<"Ens: "<<V_charm_20_1_L.Tag[i_ens]<<endl;


  //resample lattice spacing
  distr_t a_distr(UseJack), Za(UseJack), Zv(UseJack), Za_RIM(UseJack), Zv_RIM(UseJack);
  LatticeInfo L_info;
  L_info.LatInfo_new_ens(V_charm_1_L.Tag[i_ens]);
  //generate jackknife sample of input parameters
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      a_distr.distr.push_back( fm_to_inv_Gev*( L_info.a + GM()*L_info.a_err*(1.0/sqrt(Njacks-1.0))));
      Za_RIM.distr.push_back(  L_info.Za + GM()*L_info.Za_err*(1.0/sqrt(Njacks-1.0)));
      Zv_RIM.distr.push_back(  L_info.Zv + GM()*L_info.Zv_err*(1.0/sqrt(Njacks-1.0)));
      Za.distr.push_back(  L_info.Za_WI + GM()*L_info.Za_WI_err*(1.0/sqrt(Njacks-1.0)));
      Zv.distr.push_back(  L_info.Zv_WI + GM()*L_info.Zv_WI_err*(1.0/sqrt(Njacks-1.0)));
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      a_distr.distr.push_back( fm_to_inv_Gev*( L_info.a + GM()*L_info.a_err));
      Za_RIM.distr.push_back(  L_info.Za + GM()*L_info.Za_err);
      Zv_RIM.distr.push_back(  L_info.Zv + GM()*L_info.Zv_err);
      Za.distr.push_back(  L_info.Za_WI + GM()*L_info.Za_WI_err);
      Zv.distr.push_back(  L_info.Zv_WI + GM()*L_info.Zv_WI_err);
    }
  }



  Vcharm1L = Corr.corr_t(V_charm_20_1_L.col(0)[i_ens], "");
  Vcharm2L = Corr.corr_t(V_charm_20_2_L.col(0)[i_ens], "");
  Vcharm3L = Corr.corr_t(V_charm_20_3_L.col(0)[i_ens], "");
  Vcharm1M = Corr.corr_t(V_charm_20_1_M.col(0)[i_ens], "");
  Vcharm2M = Corr.corr_t(V_charm_20_2_M.col(0)[i_ens], "");
  Vcharm3M = Corr.corr_t(V_charm_20_3_M.col(0)[i_ens], "");
  distr_t_list VL = (qc*qc/3.0)*(Vcharm1L + Vcharm2L + Vcharm3L);
  distr_t_list VM = (qc*qc/3.0)*(Vcharm1M + Vcharm2M + Vcharm3M);

  distr_t_list etac_L, etac_M, Jpsi_L, Jpsi_M;

  distr_t_list P5P5_L, P5P5_M;
  P5P5_L= Corr.corr_t(P5P5_charm_20_L.col(0)[i_ens], "../data/gm2/charm/tm/P5P5_corr_Nhits20_L_"+V_charm_20_1_L.Tag[i_ens]);
  P5P5_M= Corr.corr_t(P5P5_charm_20_M.col(0)[i_ens], "../data/gm2/charm/tm/P5P5_corr_Nhits20_M_"+V_charm_20_1_L.Tag[i_ens]);
  
  etac_L = Corr.effective_mass_t(P5P5_charm_20_L.col(0)[i_ens], "../data/gm2/charm/tm/etac_Nhits20_L_"+V_charm_20_1_L.Tag[i_ens]);
  etac_M = Corr.effective_mass_t(P5P5_charm_20_M.col(0)[i_ens], "../data/gm2/charm/tm/etac_Nhits20_M_"+V_charm_20_1_L.Tag[i_ens]);
  Jpsi_L = Corr.effective_mass_t(VL, "../data/gm2/charm/tm/Jpsi_Nhits20_L_"+V_charm_20_1_L.Tag[i_ens]);
  Jpsi_M = Corr.effective_mass_t(VM, "../data/gm2/charm/tm/Jpsi_Nhits20_M_"+V_charm_20_1_L.Tag[i_ens]);

  //print vector correlator

  Print_To_File({}, {VL.ave(), VL.err(), (VL*Za*Za).ave(), (VL*Zv*Zv).err()}, "../data/gm2/charm/tm/corr_Nhits20_L_"+V_charm_20_1_L.Tag[i_ens], "", "");
  Print_To_File({}, {VM.ave(), VM.err(), (VM*Za*Za).ave(), (VM*Zv*Zv).err()}, "../data/gm2/charm/tm/corr_Nhits20_M_"+V_charm_20_1_L.Tag[i_ens], "", "");
		
		
  }

  //exit(-1);
  */




  //###########################################################################################
  //#######################  generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D ###########################
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err;
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a;
  a_A_err= a_info.a_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a;
  a_B_err= a_info.a_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a;
  a_C_err= a_info.a_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a;
  a_D_err= a_info.a_err;
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err));
    }
  }
  //#####################################################################################################################
 
  

  int Nens_light= V_light_1.size;
  int Nens_strange = V_strange_1_L.size;
  int Nens_charm = V_charm_1_L.size;
  int Nens_disco_light = disco_light.size;
  vector<string> disco_Tags;
  cout<<"N_ens light: "<<Nens_light<<endl;
  cout<<"N_ens disco_light: "<<Nens_disco_light<<endl;
  cout<<"N_ens strange: "<<Nens_strange<<endl;
  cout<<"N_ens charm: "<<Nens_charm<<endl;



  //define distr_t_list to be used in chiral+continuum analysis

  //strange and charm
  //L
  distr_t_list agm2_strange_L(UseJack), agm2_charm_L(UseJack), agm2_strange_OS_L(UseJack), agm2_charm_OS_L(UseJack);
  //M
  distr_t_list agm2_strange_M(UseJack), agm2_charm_M(UseJack), agm2_strange_OS_M(UseJack), agm2_charm_OS_M(UseJack);
  //H
  distr_t_list agm2_strange_H(UseJack), agm2_charm_H(UseJack), agm2_strange_OS_H(UseJack), agm2_charm_OS_H(UseJack);
  //Extr
  distr_t_list agm2_strange_Extr(UseJack), agm2_charm_Extr(UseJack), agm2_strange_OS_Extr(UseJack), agm2_charm_OS_Extr(UseJack);

  //strange and charm NON_ELM
  //L
  distr_t_list agm2_strange_No_ELM_L(UseJack), agm2_charm_No_ELM_L(UseJack), agm2_strange_OS_No_ELM_L(UseJack), agm2_charm_OS_No_ELM_L(UseJack);
  //M
  distr_t_list agm2_strange_No_ELM_M(UseJack), agm2_charm_No_ELM_M(UseJack), agm2_strange_OS_No_ELM_M(UseJack), agm2_charm_OS_No_ELM_M(UseJack);
  //H
  distr_t_list agm2_strange_No_ELM_H(UseJack), agm2_charm_No_ELM_H(UseJack), agm2_strange_OS_No_ELM_H(UseJack), agm2_charm_OS_No_ELM_H(UseJack);
  //Extr
  distr_t_list agm2_strange_No_ELM_Extr(UseJack), agm2_charm_No_ELM_Extr(UseJack), agm2_strange_OS_No_ELM_Extr(UseJack), agm2_charm_OS_No_ELM_Extr(UseJack);
  

  //light
  distr_t_list agm2_light(UseJack), agm2_light_fit(UseJack), agm2_light_2L_fit(UseJack), agm2_light_Lprime_fit(UseJack)  , agm2_light_infL_fit(UseJack);
  distr_t_list agm2_light_OS(UseJack), agm2_light_fit_OS(UseJack), agm2_light_2L_fit_OS(UseJack), agm2_light_Lprime_fit_OS(UseJack), agm2_light_infL_fit_OS(UseJack);

  //#############################################################################################################
  //window contributions
  //light
  distr_t_list agm2_light_W(UseJack), agm2_light_SD(UseJack), agm2_light_W_ELM(UseJack), agm2_light_SD_ELM(UseJack);
  distr_t_list agm2_light_W_OS(UseJack), agm2_light_SD_OS(UseJack), agm2_light_W_ELM_OS(UseJack), agm2_light_SD_ELM_OS(UseJack);
  //disco light
  distr_t_list agm2_disco_light_W(UseJack), agm2_disco_light_SD(UseJack), agm2_disco_light_W_ELM(UseJack), agm2_disco_light_SD_ELM(UseJack);

  //strange
  //L
  distr_t_list agm2_strange_W_L(UseJack), agm2_strange_SD_L(UseJack), agm2_strange_W_ELM_L(UseJack), agm2_strange_SD_ELM_L(UseJack);
  distr_t_list agm2_strange_W_OS_L(UseJack), agm2_strange_SD_OS_L(UseJack), agm2_strange_W_ELM_OS_L(UseJack), agm2_strange_SD_ELM_OS_L(UseJack);
  //M
  distr_t_list agm2_strange_W_M(UseJack), agm2_strange_SD_M(UseJack), agm2_strange_W_ELM_M(UseJack), agm2_strange_SD_ELM_M(UseJack);
  distr_t_list agm2_strange_W_OS_M(UseJack), agm2_strange_SD_OS_M(UseJack), agm2_strange_W_ELM_OS_M(UseJack), agm2_strange_SD_ELM_OS_M(UseJack);
  //H
  distr_t_list agm2_strange_W_H(UseJack), agm2_strange_SD_H(UseJack), agm2_strange_W_ELM_H(UseJack), agm2_strange_SD_ELM_H(UseJack);
  distr_t_list agm2_strange_W_OS_H(UseJack), agm2_strange_SD_OS_H(UseJack), agm2_strange_W_ELM_OS_H(UseJack), agm2_strange_SD_ELM_OS_H(UseJack);
  //Extr
  distr_t_list agm2_strange_W_Extr(UseJack), agm2_strange_SD_Extr(UseJack), agm2_strange_W_ELM_Extr(UseJack), agm2_strange_SD_ELM_Extr(UseJack);
  distr_t_list agm2_strange_W_OS_Extr(UseJack), agm2_strange_SD_OS_Extr(UseJack), agm2_strange_W_ELM_OS_Extr(UseJack), agm2_strange_SD_ELM_OS_Extr(UseJack);
  
  //charm
  //L
  distr_t_list agm2_charm_W_L(UseJack), agm2_charm_SD_L(UseJack), agm2_charm_W_ELM_L(UseJack), agm2_charm_SD_ELM_L(UseJack);
  distr_t_list agm2_charm_W_OS_L(UseJack), agm2_charm_SD_OS_L(UseJack), agm2_charm_W_ELM_OS_L(UseJack), agm2_charm_SD_ELM_OS_L(UseJack);
  //M
  distr_t_list agm2_charm_W_M(UseJack), agm2_charm_SD_M(UseJack), agm2_charm_W_ELM_M(UseJack), agm2_charm_SD_ELM_M(UseJack);
  distr_t_list agm2_charm_W_OS_M(UseJack), agm2_charm_SD_OS_M(UseJack), agm2_charm_W_ELM_OS_M(UseJack), agm2_charm_SD_ELM_OS_M(UseJack);
  //H
  distr_t_list agm2_charm_W_H(UseJack), agm2_charm_SD_H(UseJack), agm2_charm_W_ELM_H(UseJack), agm2_charm_SD_ELM_H(UseJack);
  distr_t_list agm2_charm_W_OS_H(UseJack), agm2_charm_SD_OS_H(UseJack), agm2_charm_W_ELM_OS_H(UseJack), agm2_charm_SD_ELM_OS_H(UseJack);
  //Extr
  distr_t_list agm2_charm_W_Extr(UseJack), agm2_charm_SD_Extr(UseJack), agm2_charm_W_ELM_Extr(UseJack), agm2_charm_SD_ELM_Extr(UseJack);
  distr_t_list agm2_charm_W_OS_Extr(UseJack), agm2_charm_SD_OS_Extr(UseJack), agm2_charm_W_ELM_OS_Extr(UseJack), agm2_charm_SD_ELM_OS_Extr(UseJack);
  //#############################################################################################################

 
  
  vector<vector<fit_par>> par_list_anal_repr;//only used for light quark contribution with opposite r 
  vector<vector<fit_par>> par_list_anal_repr_OS; //only used for light quark contribution with same r

  //light MV and ZV
  distr_t_list MV_fit_light(UseJack), MV_fit_light_OS(UseJack), ZV_fit_light(UseJack), ZV_fit_light_OS(UseJack);

  //strange and charm MV and ZV
  //L
  distr_t_list MV_fit_strange_L(UseJack), MV_fit_charm_L(UseJack), MV_fit_strange_OS_L(UseJack), MV_fit_charm_OS_L(UseJack);
  distr_t_list ZV_fit_strange_L(UseJack), ZV_fit_charm_L(UseJack), ZV_fit_strange_OS_L(UseJack), ZV_fit_charm_OS_L(UseJack);
  //M
  distr_t_list MV_fit_strange_M(UseJack), MV_fit_charm_M(UseJack), MV_fit_strange_OS_M(UseJack), MV_fit_charm_OS_M(UseJack);
  distr_t_list ZV_fit_strange_M(UseJack), ZV_fit_charm_M(UseJack), ZV_fit_strange_OS_M(UseJack), ZV_fit_charm_OS_M(UseJack);
  //H
  distr_t_list MV_fit_strange_H(UseJack), MV_fit_charm_H(UseJack), MV_fit_strange_OS_H(UseJack), MV_fit_charm_OS_H(UseJack);
  distr_t_list ZV_fit_strange_H(UseJack), ZV_fit_charm_H(UseJack), ZV_fit_strange_OS_H(UseJack), ZV_fit_charm_OS_H(UseJack);

  //kaon mass fit
  distr_t_list MK_fit_L(UseJack), MK_fit_M(UseJack), MK_fit_H(UseJack);
  //D mass fit
  distr_t_list MD_fit_L(UseJack), MD_fit_M(UseJack), MD_fit_H(UseJack);
  
  distr_t_list Mpi_fit_strange(UseJack), Mpi_OS_fit_strange(UseJack), Mpi_fit_charm(UseJack), Mpi_OS_fit_charm(UseJack);
  distr_t_list Mpi_fit(UseJack), Mpi_OS_fit(UseJack),  fp_fit(UseJack), Mpi_fit_disco(UseJack), Mpi_OS_fit_disco(UseJack), fp_fit_disco(UseJack);
  distr_t_list Zv_fit(UseJack), Za_fit(UseJack), Zp_ov_Zs_fit(UseJack), Zv_fit_disco(UseJack);
  distr_t_list Zv_fit_strange(UseJack), Za_fit_strange(UseJack);
  distr_t_list Zv_fit_strange_heavy(UseJack), Za_fit_strange_heavy(UseJack);
  distr_t_list Zv_fit_charm_L(UseJack), Za_fit_charm_L(UseJack);
  distr_t_list Zv_fit_charm_M(UseJack), Za_fit_charm_M(UseJack);
  distr_t_list Zv_fit_charm_H(UseJack), Za_fit_charm_H(UseJack);
  distr_t_list Zv_RIMOM(UseJack), Za_RIMOM(UseJack), Zv_WI(UseJack), Za_WI(UseJack);
  Vfloat L_list, a_list, ml_list, L_list_disco, a_list_disco, ml_list_disco;
  distr_t_list a_distr_list(UseJack), a_distr_list_strange(UseJack), a_distr_list_charm(UseJack);
  Vfloat L_strange_list, a_strange_list, ml_strange_list;
  Vfloat L_charm_list, a_charm_list, ml_charm_list;

  //define lambda for convolution with kernel
  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
  auto exp_MV = [&](double Mv, double t, double size) -> double { return exp(-Mv*t);};
  
  
  //strange
  channel="s";
  for(int i_ens=0;i_ens<Nens_strange;i_ens++) {
    
  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = V_strange_1_L.nrows[i_ens];


  //resample lattice spacing
  distr_t a_distr(UseJack), Za(UseJack), Zv(UseJack), Za_RIM(UseJack), Zv_RIM(UseJack);
  distr_t Zv_heavy(UseJack), Za_heavy(UseJack);
  LatticeInfo L_info;
  L_info.LatInfo_new_ens(V_strange_1_L.Tag[i_ens]);
  //generate jackknife sample of input parameters
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      Za_RIM.distr.push_back(  L_info.Za + GM()*L_info.Za_err*(1.0/sqrt(Njacks-1.0)));
      Zv_RIM.distr.push_back(  L_info.Zv + GM()*L_info.Zv_err*(1.0/sqrt(Njacks-1.0)));
      Za.distr.push_back(  L_info.Za_WI + GM()*L_info.Za_WI_err*(1.0/sqrt(Njacks-1.0)));
      Zv.distr.push_back(  L_info.Zv_WI + GM()*L_info.Zv_WI_err*(1.0/sqrt(Njacks-1.0)));
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      Za_RIM.distr.push_back(  L_info.Za + GM()*L_info.Za_err);
      Zv_RIM.distr.push_back(  L_info.Zv + GM()*L_info.Zv_err);
      Za.distr.push_back(  L_info.Za_WI + GM()*L_info.Za_WI_err);
      Zv.distr.push_back(  L_info.Zv_WI + GM()*L_info.Zv_WI_err);
    }
  }
  
  //push_back lattice info
  L_strange_list.push_back(L_info.L);
  a_strange_list.push_back(L_info.a);
  ml_strange_list.push_back(L_info.ml);

  //push_back Zv and Za RI-MOM
  Zv_RIMOM.distr_list.push_back(Zv_RIM);
  Za_RIMOM.distr_list.push_back(Za_RIM);
  Zv_WI.distr_list.push_back(Zv);
  Za_WI.distr_list.push_back(Za);
  Zv_heavy = Zv;
  Za_heavy = Za;

  int Tmin_P5P5;
  int Tmax_P5P5;
  int Tmin_VV;
  int Tmax_VV;
  int Tmin_VV_OS;
  int Tmax_VV_OS;


  //set time intervals for pseudoscalar obs
  if(V_strange_1_L.Tag[i_ens].substr(1,1) == "C") {
    a_distr = a_C;
    if(V_strange_1_L.Tag[i_ens]=="cC211a.06.80") { Corr.Tmin=40; Corr.Tmax=60;}
    else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
  }
  else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "B") {
    a_distr = a_B;
    if(V_strange_1_L.Tag[i_ens]== "cB211a.14.64") {Corr.Tmin=31; Corr.Tmax=58;}
    else if(V_strange_1_L.Tag[i_ens] == "cB211a.25.48") {Corr.Tmin=23;Corr.Tmax=44;}
    else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.64") {Corr.Tmin=35; Corr.Tmax= 55;}
    else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.96") {Corr.Tmin=35; Corr.Tmax= 60;}
    else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
  }
  else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "A") {
    a_distr = a_A;
    if(V_strange_1_L.Tag[i_ens] == "cA211a.12.48") {Corr.Tmin=19; Corr.Tmax=33;}
    else if(V_strange_1_L.Tag[i_ens] == "cA211a.40.24") {Corr.Tmin=16; Corr.Tmax=22;}
    else if(V_strange_1_L.Tag[i_ens] == "cA211a.53.24") {Corr.Tmin=16; Corr.Tmax=22;}
    else if(V_strange_1_L.Tag[i_ens] == "cA211ab.30.32") {Corr.Tmin=21; Corr.Tmax=30;}
    else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
  }
  else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "D") {
    a_distr = a_D;
    if(V_strange_1_L.Tag[i_ens] == "cD211a.054.96") {Corr.Tmin=50; Corr.Tmax=75;}
    else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
  }
  else crash("Ensemble tag not valid");


  a_distr_list_strange.distr_list.push_back(a_distr);

  //set Tmin_P5P5 and Tmax_P5P5 to the values Corr.Tmin and Corr.Tmax
  Tmin_P5P5 = Corr.Tmin;
  Tmax_P5P5 = Corr.Tmax;
  
  //set time intervals for vector obs
  if(V_strange_1_L.Tag[i_ens].substr(1,1) == "C") {
    if(V_strange_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_VV=25; Tmax_VV=35; Tmin_VV_OS=24; Tmax_VV_OS= 33;}
    else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
  }
  else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "B") {
    if(V_strange_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_VV=18; Tmax_VV=28; Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
    else if(V_strange_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_VV=19;Tmax_VV=25;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
    else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_VV=25; Tmax_VV= 32; Tmin_VV_OS= 25; Tmax_VV_OS=32;}
    else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.96") {Tmin_VV=25; Tmax_VV= 32; Tmin_VV_OS=25; Tmax_VV_OS=32;}
    else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
  }
  else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "A") {
    if(V_strange_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_VV=18; Tmax_VV=28;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
    else if(V_strange_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_VV=14; Tmax_VV=20;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
    else if(V_strange_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_VV=12; Tmax_VV=19;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
    else if(V_strange_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_VV=17;Tmax_VV=23;  Tmin_VV_OS=Tmin_VV; Tmax_VV_OS=Tmax_VV;}
    else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
  }
  else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "D") {
    if(V_strange_1_L.Tag[i_ens]=="cD211a.054.96") { Tmin_VV=32; Tmax_VV=40;  Tmin_VV_OS=36; Tmax_VV_OS=45;}
    else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
  }
  else crash("Ensemble tag not valid");

  

  
  
  //tm
  //L
  distr_t_list  V_strange_1_distr_L, V_strange_2_distr_L, V_strange_3_distr_L;
  distr_t_list  V_strange_distr_L, MV_strange_distr_L, ZV_strange_distr_L;
  distr_t_list M_etas_distr, M_etas_distr_heavy, MK_distr_L, FK_distr_L;
  distr_t MV_strange_L, ZV_strange_L, M_etas, M_etas_heavy, MK_L, FK_L;
  //M
  distr_t_list  V_strange_1_distr_M, V_strange_2_distr_M, V_strange_3_distr_M;
  distr_t_list  V_strange_distr_M, MV_strange_distr_M, ZV_strange_distr_M;
  distr_t_list  MK_distr_M, FK_distr_M;
  distr_t MV_strange_M, ZV_strange_M, MK_M, FK_M;
  //H
  distr_t_list  V_strange_1_distr_H, V_strange_2_distr_H, V_strange_3_distr_H;
  distr_t_list  V_strange_distr_H, MV_strange_distr_H, ZV_strange_distr_H;
  distr_t_list  MK_distr_H, FK_distr_H;
  distr_t MV_strange_H, ZV_strange_H, MK_H, FK_H;
  


  //OS
  //L
  distr_t_list  V_strange_OS_1_distr_L, V_strange_OS_2_distr_L, V_strange_OS_3_distr_L;
  distr_t_list  V_strange_OS_distr_L, MV_strange_OS_distr_L, ZV_strange_OS_distr_L;
  distr_t_list M_etas_OS_distr, M_etas_OS_distr_heavy;
  distr_t MV_strange_OS_L, ZV_strange_OS_L, M_etas_OS, M_etas_OS_heavy;
  //M
  distr_t_list  V_strange_OS_1_distr_M, V_strange_OS_2_distr_M, V_strange_OS_3_distr_M;
  distr_t_list  V_strange_OS_distr_M, MV_strange_OS_distr_M, ZV_strange_OS_distr_M;
  distr_t MV_strange_OS_M, ZV_strange_OS_M;
  //H
  distr_t_list  V_strange_OS_1_distr_H, V_strange_OS_2_distr_H, V_strange_OS_3_distr_H;
  distr_t_list  V_strange_OS_distr_H, MV_strange_OS_distr_H, ZV_strange_OS_distr_H;
  distr_t MV_strange_OS_H, ZV_strange_OS_H;


  //To compute Za and Zv (hadronic method)
  distr_t_list overlap_P5P5_distr, fetas_distr, overlap_P5P5_distr_heavy, fetas_distr_heavy;
  distr_t_list etas_corr, P5P5_OS_distr;
  distr_t_list etas_corr_heavy, P5P5_OS_distr_heavy;
  distr_t_list overlap_P5P5_OS_distr, overlap_P5P5_OS_distr_heavy;
  distr_t_list ratio_P5P5_overlap_OS_tm,  Zp_ov_Zs_distr, ratio_P5P5_overlap_OS_tm_heavy, Zp_ov_Zs_distr_heavy;
  distr_t_list A0P5_distr, A0P5_OS_distr, A0P5_distr_heavy, A0P5_OS_distr_heavy;
  distr_t_list RV, RA, RV_heavy, RA_heavy;
  distr_t Zv_hadr, Za_hadr, Zp_ov_Zs;
  distr_t Zv_hadr_heavy, Za_hadr_heavy, Zp_ov_Zs_heavy;


   
  //tm
  
  etas_corr= Corr.corr_t(corr_P5P5_strange.col(0)[i_ens], "");
  M_etas_distr= Corr.effective_mass_t(corr_P5P5_strange.col(0)[i_ens], "../data/gm2/strange/tm/P5P5_ss_mass"+V_strange_1_L.Tag[i_ens]+".dat");
  etas_corr_heavy = Corr.corr_t(corr_P5P5_strange_heavy.col(0)[i_ens], "");
  M_etas_distr_heavy = Corr.effective_mass_t(corr_P5P5_strange_heavy.col(0)[i_ens], "../data/gm2/strange/tm/P5P5_ss_mass_heavy"+V_strange_1_L.Tag[i_ens]+".dat");
  
  
  overlap_P5P5_distr = Corr.residue_t(corr_P5P5_strange.col(0)[i_ens], "");
  overlap_P5P5_distr_heavy = Corr.residue_t(corr_P5P5_strange_heavy.col(0)[i_ens], "");
  fetas_distr = 2.0*L_info.ms_L*Corr.decay_constant_t(corr_P5P5_strange.col(0)[i_ens], "");
  fetas_distr_heavy= 2.0*L_info.ms_M*Corr.decay_constant_t(corr_P5P5_strange_heavy.col(0)[i_ens], "");
 
  //L
  V_strange_1_distr_L = Corr.corr_t(V_strange_1_L.col(0)[i_ens], "../data/gm2/strange/tm/corr_1_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
  V_strange_2_distr_L = Corr.corr_t(V_strange_2_L.col(0)[i_ens], "../data/gm2/strange/tm/corr_2_"+V_strange_2_L.Tag[i_ens]+"_L.dat");
  V_strange_3_distr_L = Corr.corr_t(V_strange_3_L.col(0)[i_ens], "../data/gm2/strange/tm/corr_3_"+V_strange_3_L.Tag[i_ens]+"_L.dat");
  /*
  MK_distr_L = Corr.effective_mass_t(pt2_K_L.col(0)[i_ens], "../data/gm2/strange/tm/MK_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
  FK_distr_L = (L_info.ml + L_info.ms_L)*Corr.decay_constant_t(pt2_K_L.col(0)[i_ens], "../data/gm2/strange/tm/FK_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
  */
  
  //M
  V_strange_1_distr_M = Corr.corr_t(V_strange_1_M.col(0)[i_ens], "../data/gm2/strange/tm/corr_1_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
  V_strange_2_distr_M = Corr.corr_t(V_strange_2_M.col(0)[i_ens], "../data/gm2/strange/tm/corr_2_"+V_strange_2_M.Tag[i_ens]+"_M.dat");
  V_strange_3_distr_M = Corr.corr_t(V_strange_3_M.col(0)[i_ens], "../data/gm2/strange/tm/corr_3_"+V_strange_3_M.Tag[i_ens]+"_M.dat");
  /*
  MK_distr_M = Corr.effective_mass_t(pt2_K_M.col(0)[i_ens], "../data/gm2/strange/tm/MK_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
  FK_distr_M = (L_info.ml + L_info.ms_M)*Corr.decay_constant_t(pt2_K_M.col(0)[i_ens], "../data/gm2/strange/tm/FK_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
  */
  
  //H
  /*
  V_strange_1_distr_H = Corr.corr_t(V_strange_1_H.col(0)[i_ens], "../data/gm2/strange/tm/corr_1_"+V_strange_1_H.Tag[i_ens]+"_H.dat");
  V_strange_2_distr_H = Corr.corr_t(V_strange_2_H.col(0)[i_ens], "../data/gm2/strange/tm/corr_2_"+V_strange_2_H.Tag[i_ens]+"_H.dat");
  V_strange_3_distr_H = Corr.corr_t(V_strange_3_H.col(0)[i_ens], "../data/gm2/strange/tm/corr_3_"+V_strange_3_H.Tag[i_ens]+"_H.dat");
  MK_distr_H = Corr.effective_mass_t(pt2_K_H.col(0)[i_ens], "../data/gm2/strange/tm/MK_"+V_strange_1_H.Tag[i_ens]+"_H.dat");
  FK_distr_H = (L_info.ml + L_info.ms_H)*Corr.decay_constant_t(pt2_K_H.col(0)[i_ens], "../data/gm2/strange/tm/FK_"+V_strange_1_H.Tag[i_ens]+"_H.dat");
  */
 
  //OS
 
  P5P5_OS_distr = Corr.corr_t(corr_P5P5_OS_strange.col(0)[i_ens], "");
  P5P5_OS_distr_heavy = Corr.corr_t(corr_P5P5_OS_strange_heavy.col(0)[i_ens], "");
  M_etas_OS_distr= Corr.effective_mass_t(corr_P5P5_OS_strange.col(0)[i_ens], "../data/gm2/strange/OS/P5P5_ss_mass_"+V_strange_1_L.Tag[i_ens]+".dat");
  M_etas_OS_distr_heavy = Corr.effective_mass_t(corr_P5P5_OS_strange_heavy.col(0)[i_ens], "../data/gm2/strange/OS/P5P5_ss_mass_heavy"+V_strange_1_L.Tag[i_ens]+".dat");

 
  overlap_P5P5_OS_distr= Corr.residue_t(corr_P5P5_OS_strange.col(0)[i_ens], "");
  overlap_P5P5_OS_distr_heavy = Corr.residue_t(corr_P5P5_OS_strange_heavy.col(0)[i_ens], "");

  //L
  V_strange_OS_1_distr_L = Corr.corr_t(V_strange_OS_1_L.col(0)[i_ens], "../data/gm2/strange/OS/corr_1_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
  V_strange_OS_2_distr_L = Corr.corr_t(V_strange_OS_2_L.col(0)[i_ens], "../data/gm2/strange/OS/corr_2_"+V_strange_2_L.Tag[i_ens]+"_L.dat");
  V_strange_OS_3_distr_L = Corr.corr_t(V_strange_OS_3_L.col(0)[i_ens], "../data/gm2/strange/OS/corr_3_"+V_strange_3_L.Tag[i_ens]+"_L.dat");
  //M
  V_strange_OS_1_distr_M = Corr.corr_t(V_strange_OS_1_M.col(0)[i_ens], "../data/gm2/strange/OS/corr_1_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
  V_strange_OS_2_distr_M = Corr.corr_t(V_strange_OS_2_M.col(0)[i_ens], "../data/gm2/strange/OS/corr_2_"+V_strange_2_M.Tag[i_ens]+"_M.dat");
  V_strange_OS_3_distr_M = Corr.corr_t(V_strange_OS_3_M.col(0)[i_ens], "../data/gm2/strange/OS/corr_3_"+V_strange_3_M.Tag[i_ens]+"_M.dat");
  //H
  /*
  V_strange_OS_1_distr_H = Corr.corr_t(V_strange_OS_1_H.col(0)[i_ens], "../data/gm2/strange/OS/corr_1_"+V_strange_1_H.Tag[i_ens]+"_H.dat");
  V_strange_OS_2_distr_H = Corr.corr_t(V_strange_OS_2_H.col(0)[i_ens], "../data/gm2/strange/OS/corr_2_"+V_strange_2_H.Tag[i_ens]+"_H.dat");
  V_strange_OS_3_distr_H = Corr.corr_t(V_strange_OS_3_H.col(0)[i_ens], "../data/gm2/strange/OS/corr_3_"+V_strange_3_H.Tag[i_ens]+"_H.dat");
  */



   
  //#######################################  COMPUTATION OF ZV AND ZA (Hadronic method) ##################################
  //define lambda functions to be used
  auto sqr= [=](double a, double b) {return sqrt(a);};
  auto SINH= [](double m) -> double  {return sinh(m);};

  
  //take ratio between OS and tm pion amplitude to compute Zp/Zs RC.
  ratio_P5P5_overlap_OS_tm= overlap_P5P5_OS_distr/overlap_P5P5_distr;
  Zp_ov_Zs_distr = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm);
  ratio_P5P5_overlap_OS_tm_heavy= overlap_P5P5_OS_distr_heavy/overlap_P5P5_distr_heavy;
  Zp_ov_Zs_distr_heavy = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm_heavy);
  

  
  //antysymmetrize w.r.t. t -> T-t for A0P5 correlators
  //Corr.Reflection_sign = -1;
  A0P5_distr= Corr.corr_t(corr_A0P5_strange.col(0)[i_ens], "");
  A0P5_OS_distr = Corr.corr_t(corr_A0P5_OS_strange.col(0)[i_ens], "");
  A0P5_distr_heavy = Corr.corr_t(corr_A0P5_strange_heavy.col(0)[i_ens], "");
  A0P5_OS_distr_heavy = Corr.corr_t(corr_A0P5_OS_strange_heavy.col(0)[i_ens],"");


  
  //restore symmetrization
  Corr.Reflection_sign = 1;

  //compute RV (estimator for Zv)


  RV= 2.0*L_info.ms_L*etas_corr/distr_t_list::derivative(A0P5_distr, 0); //central derivative
  RV_heavy = 2.0*L_info.ms_M*etas_corr_heavy/distr_t_list::derivative(A0P5_distr_heavy, 0); //central derivative

   
  //tm and OS P5P5
  M_etas = Corr.Fit_distr(M_etas_distr);
  M_etas_heavy = Corr.Fit_distr(M_etas_distr_heavy);
  M_etas_OS= Corr.Fit_distr(M_etas_OS_distr);
  M_etas_OS_heavy = Corr.Fit_distr(M_etas_OS_distr_heavy);
  //fit obs to compute Zv and Za (hadronic method)
  Zp_ov_Zs = Corr.Fit_distr(Zp_ov_Zs_distr);
  Zp_ov_Zs_heavy = Corr.Fit_distr(Zp_ov_Zs_distr_heavy);
  //set plateaux for RV
  Corr.Tmin= 20;
  Corr.Tmax= Corr.Nt/2 -2;
  if(V_strange_1_L.Tag[i_ens].substr(1,1) == "C") { Corr.Tmin=30; Corr.Tmax=60;}
  if(V_strange_1_L.Tag[i_ens].substr(1,1) == "B") {Corr.Tmin=25;}
  if(V_strange_1_L.Tag[i_ens].substr(1,1) == "D") { Corr.Tmin=30; Corr.Tmax=60;}
  Zv_hadr= Corr.Fit_distr(RV);
  Zv_hadr_heavy = Corr.Fit_distr(RV_heavy);
  RA = 2.0*L_info.ms_L*(P5P5_OS_distr/distr_t_list::derivative(A0P5_OS_distr, 0))*(M_etas_OS/M_etas)*(distr_t::f_of_distr(SINH, M_etas_OS)/distr_t::f_of_distr(SINH, M_etas))*(1.0/Zp_ov_Zs);
  RA_heavy = 2.0*L_info.ms_M*(P5P5_OS_distr_heavy/distr_t_list::derivative(A0P5_OS_distr_heavy, 0))*(M_etas_OS_heavy/M_etas_heavy)*(distr_t::f_of_distr(SINH, M_etas_OS_heavy)/distr_t::f_of_distr(SINH, M_etas_heavy))*(1.0/Zp_ov_Zs_heavy);
  //set plateaux for RA
  int Tmin_RA=0;
  int Tmax_RA=0;
  //set time intervals for RA
  if(V_strange_1_L.Tag[i_ens].substr(1,1) == "C") {
    if(V_strange_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_RA=25; Tmax_RA=50;}
    else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
  }
  else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "B") {
    if(V_strange_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_RA=18; Tmax_RA=26;}
    else if(V_strange_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_RA=10;Tmax_RA=21;}
    else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_RA=20; Tmax_RA= 50;}
    else if(V_strange_1_L.Tag[i_ens] == "cB211b.072.96") {Tmin_RA=20; Tmax_RA= 50;}
    else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
  }
  else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "A") {
    if(V_strange_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_RA=15; Tmax_RA=25;}
    else if(V_strange_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_RA=12; Tmax_RA=20;}
    else if(V_strange_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_RA=11; Tmax_RA=20;}
    else if(V_strange_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_RA=14;Tmax_RA=22;}
    else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
  }
  else if(V_strange_1_L.Tag[i_ens].substr(1,1) == "D") {
    if(V_strange_1_L.Tag[i_ens] == "cD211a.054.96") {Tmin_RA=32; Tmax_RA=90;}
    else crash("Cannot find ensemble tag: "+V_strange_1_L.Tag[i_ens]);
  }
  else crash("Ensemble tag not valid");

  Corr.Tmin=Tmin_RA;
  Corr.Tmax=Tmax_RA;
  
  Za_hadr= Corr.Fit_distr(RA);
  Za_hadr_heavy = Corr.Fit_distr(RA_heavy);

  //print Rv and RA
  //print RV
  Print_To_File({}, {RV.ave(), RV.err(), RV_heavy.ave(), RV_heavy.err()}, "../data/gm2/strange/RV_ss"+V_strange_1_L.Tag[i_ens]+".dat.t", "", "");
  //print RA
  Print_To_File({}, {RA.ave(), RA.err(), RA_heavy.ave(), RA_heavy.err()}, "../data/gm2/strange/RA_ss"+V_strange_1_L.Tag[i_ens]+".dat.t", "", "");
  //print Zp_ov_Zs
  Print_To_File({}, {Zp_ov_Zs_distr.ave(), Zp_ov_Zs_distr.err(), Zp_ov_Zs_distr_heavy.ave(), Zp_ov_Zs_distr_heavy.err()}, "../data/gm2/strange/Zp_ov_Zs_ss"+V_strange_1_L.Tag[i_ens]+".dat.t", "", "");

  //push_back

  Zv_fit_strange.distr_list.push_back(Zv_hadr);
  Za_fit_strange.distr_list.push_back(Za_hadr);
  Zv_fit_strange_heavy.distr_list.push_back(Zv_hadr_heavy);
  Za_fit_strange_heavy.distr_list.push_back(Za_hadr_heavy);

  if(Use_Za_Zv_from_strange_run) { Zv= Zv_hadr; Za=Za_hadr; Zv_heavy = Zv_hadr_heavy; Za_heavy= Za_hadr_heavy;}

  //################################################ END OF COMPUTATION OF ZV AND ZA (Hadronic method) #####################################à



 
  //sum over the Lorenz indices of the e.m. current
  //L
  //tm
  V_strange_distr_L= (pow(qs,2)/3.0)*(V_strange_1_distr_L+ V_strange_2_distr_L + V_strange_3_distr_L);
  //OS
  V_strange_OS_distr_L= (pow(qs,2)/3.0)*(V_strange_OS_1_distr_L+ V_strange_OS_2_distr_L + V_strange_OS_3_distr_L);
  
  //M
  //tm
  V_strange_distr_M= (pow(qs,2)/3.0)*(V_strange_1_distr_M+ V_strange_2_distr_M + V_strange_3_distr_M);
  //OS
  V_strange_OS_distr_M= (pow(qs,2)/3.0)*(V_strange_OS_1_distr_M+ V_strange_OS_2_distr_M + V_strange_OS_3_distr_M);

  //H
  /*
  //tm
  V_strange_distr_H= (pow(qs,2)/3.0)*(V_strange_1_distr_H+ V_strange_2_distr_H + V_strange_3_distr_H);
  //OS
  V_strange_OS_distr_H= (pow(qs,2)/3.0)*(V_strange_OS_1_distr_H+ V_strange_OS_2_distr_H + V_strange_OS_3_distr_H);
  */



  //free corr LO artifacts
  Vfloat free_corr_log_art(Corr.Nt, 0.0);
  for(int t=0;t<Corr.Nt;t++) { if( t*a_distr.ave()/fm_to_inv_Gev < 1.0 && t != 0) {   free_corr_log_art[t] = -1.0*(qs*qs)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));} else free_corr_log_art[t] = 0.0;}

  //L
  distr_t_list V_strange_distr_L_pert_sub = V_strange_distr_L + free_corr_log_art;
  distr_t_list V_strange_OS_distr_L_pert_sub = V_strange_OS_distr_L + free_corr_log_art;
  //M
  distr_t_list V_strange_distr_M_pert_sub = V_strange_distr_M + free_corr_log_art;
  distr_t_list V_strange_OS_distr_M_pert_sub = V_strange_OS_distr_M + free_corr_log_art;
  
  // print summed correlators to file
  //L
  //tm
  Print_To_File({}, {V_strange_distr_L.ave(), V_strange_distr_L.err(), (Za*Za*V_strange_distr_L).ave(), (Za*Za*V_strange_distr_L).err()}, "../data/gm2/strange/tm/corr_sum_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "# t    bare   renormalized");
  //OS
  Print_To_File({}, {V_strange_OS_distr_L.ave(), V_strange_OS_distr_L.err(), (Zv*Zv*V_strange_OS_distr_L).ave(), (Zv*Zv*V_strange_OS_distr_L).err()}, "../data/gm2/strange/OS/corr_sum_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare   renormalized");

  //M
  //tm
  Print_To_File({}, {V_strange_distr_M.ave(), V_strange_distr_M.err(), (Za_heavy*Za_heavy*V_strange_distr_M).ave(), (Za_heavy*Za_heavy*V_strange_distr_M).err()}, "../data/gm2/strange/tm/corr_sum_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "# t    bare   renormalized");
  //OS
  Print_To_File({}, {V_strange_OS_distr_M.ave(), V_strange_OS_distr_M.err(), (Zv_heavy*Zv_heavy*V_strange_OS_distr_M).ave(), (Zv_heavy*Zv_heavy*V_strange_OS_distr_M).err()}, "../data/gm2/strange/OS/corr_sum_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#t   bare   renormalized");

  //H
  /*
  //tm
  Print_To_File({}, {V_strange_distr_H.ave(), V_strange_distr_H.err(), (Za*Za*V_strange_distr_H).ave(), (Za*Za*V_strange_distr_H).err()}, "../data/gm2/strange/tm/corr_sum_"+V_strange_1_H.Tag[i_ens]+"_H.dat.t", "", "# t    bare   renormalized");
  //OS
  Print_To_File({}, {V_strange_OS_distr_H.ave(), V_strange_OS_distr_H.err(), (Zv*Zv*V_strange_OS_distr_H).ave(), (Zv*Zv*V_strange_OS_distr_H).err()}, "../data/gm2/strange/OS/corr_sum_"+V_strange_1_H.Tag[i_ens]+"_H.dat.t", "", "#t   bare   renormalized");
  */
  
  
  //extract effective masses, overlap from V and fit

 

  
  //L
  /*
  Corr.Tmin=Tmin_P5P5;
  Corr.Tmax=Tmax_P5P5;
  MK_L= Corr.Fit_distr(MK_distr_L);
  FK_L= Corr.Fit_distr(FK_distr_L);
  */
  Corr.Tmin=Tmin_VV;
  Corr.Tmax=Tmax_VV;
  MV_strange_distr_L= Corr.effective_mass_t(V_strange_distr_L, "../data/gm2/strange/tm/MV_mass_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
  ZV_strange_distr_L= Corr.residue_t(V_strange_distr_L, "../data/gm2/strange/tm/ZV_overlap_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
  MV_strange_L = Corr.Fit_distr(MV_strange_distr_L);
  ZV_strange_L = Corr.Fit_distr(ZV_strange_distr_L);
  
  //M
  /*
  Corr.Tmin=Tmin_P5P5;
  Corr.Tmax=Tmax_P5P5;
  MK_M= Corr.Fit_distr(MK_distr_M);
  FK_M= Corr.Fit_distr(FK_distr_M);
  */
  Corr.Tmin=Tmin_VV;
  Corr.Tmax=Tmax_VV;
  MV_strange_distr_M= Corr.effective_mass_t(V_strange_distr_M, "../data/gm2/strange/tm/MV_mass_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
  ZV_strange_distr_M= Corr.residue_t(V_strange_distr_M, "../data/gm2/strange/tm/ZV_overlap_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
  MV_strange_M = Corr.Fit_distr(MV_strange_distr_M);
  ZV_strange_M = Corr.Fit_distr(ZV_strange_distr_M);
 
  //H
  /*
  Corr.Tmin=Tmin_P5P5;
  Corr.Tmax=Tmax_P5P5;
  MK_H= Corr.Fit_distr(MK_distr_H);
  FK_H= Corr.Fit_distr(FK_distr_H);
  Corr.Tmin=Tmin_VV;
  Corr.Tmax=Tmax_VV;
  MV_strange_distr_H= Corr.effective_mass_t(V_strange_distr_H, "../data/gm2/strange/tm/MV_mass_"+V_strange_1_H.Tag[i_ens]+"_H.dat");
  ZV_strange_distr_H= Corr.residue_t(V_strange_distr_H, "../data/gm2/strange/tm/ZV_overlap_"+V_strange_1_H.Tag[i_ens]+"_H.dat");
  MV_strange_H = Corr.Fit_distr(MV_strange_distr_H);
  ZV_strange_H = Corr.Fit_distr(ZV_strange_distr_H);
  */
  

  //push_back MV_strange, ZV_strange, Mpi and MK
  //Mpi_fit_strange.distr_list.push_back(Mpi);
  //L
  MV_fit_strange_L.distr_list.push_back(MV_strange_L);
  ZV_fit_strange_L.distr_list.push_back(Za*Za*ZV_strange_L); 
  //MK_fit_L.distr_list.push_back(MK_L);
  //M
  MV_fit_strange_M.distr_list.push_back(MV_strange_M);
  ZV_fit_strange_M.distr_list.push_back(Za_heavy*Za_heavy*ZV_strange_M); 
  //MK_fit_M.distr_list.push_back(MK_M);
  //H
  /*
  MV_fit_strange_H.distr_list.push_back(MV_strange_H);
  ZV_fit_strange_H.distr_list.push_back(Za*Za*ZV_strange_H); 
  MK_fit_H.distr_list.push_back(MK_H);
  */


  //OS
  
  //L
  Corr.Tmin= Tmin_VV_OS;
  Corr.Tmax= Tmax_VV_OS;
  MV_strange_OS_distr_L= Corr.effective_mass_t(V_strange_OS_distr_L, "../data/gm2/strange/OS/MV_mass_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
  ZV_strange_OS_distr_L= Corr.residue_t(V_strange_OS_distr_L, "../data/gm2/strange/OS/ZV_overlap_"+V_strange_1_L.Tag[i_ens]+"_L.dat");
  MV_strange_OS_L = Corr.Fit_distr(MV_strange_OS_distr_L);
  ZV_strange_OS_L = Corr.Fit_distr(ZV_strange_OS_distr_L);
  //M
  MV_strange_OS_distr_M= Corr.effective_mass_t(V_strange_OS_distr_M, "../data/gm2/strange/OS/MV_mass_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
  ZV_strange_OS_distr_M= Corr.residue_t(V_strange_OS_distr_M, "../data/gm2/strange/OS/ZV_overlap_"+V_strange_1_M.Tag[i_ens]+"_M.dat");
  MV_strange_OS_M = Corr.Fit_distr(MV_strange_OS_distr_M);
  ZV_strange_OS_M = Corr.Fit_distr(ZV_strange_OS_distr_M);
  //H
  /*
  MV_strange_OS_distr_H= Corr.effective_mass_t(V_strange_OS_distr_H, "../data/gm2/strange/OS/MV_mass_"+V_strange_1_H.Tag[i_ens]+"_H.dat");
  ZV_strange_OS_distr_H= Corr.residue_t(V_strange_OS_distr_H, "../data/gm2/strange/OS/ZV_overlap_"+V_strange_1_H.Tag[i_ens]+"_H.dat");
  MV_strange_OS_H = Corr.Fit_distr(MV_strange_OS_distr_H);
  ZV_strange_OS_H = Corr.Fit_distr(ZV_strange_OS_distr_H);
  */

 
  //define Csi_K vector
  //vector<distr_t> Csi_K({MK_L*MK_L/(FK_L*FK_L), MK_M*MK_M/(FK_M*FK_M), MK_H*MK_H/(FK_H*FK_H)});
  vector<distr_t> Mphi2( {MV_strange_L*MV_strange_L/(a_distr*a_distr), MV_strange_M*MV_strange_M/(a_distr*a_distr)});
  if(!Use_phi_etas) { Mphi2 = { M_etas*M_etas/(a_distr*a_distr), M_etas_heavy*M_etas_heavy/(a_distr*a_distr)};}

  //Generate fake M_etas jack distribution

  distr_t m_etas_phys_distr;
  distr_t m_phi_phys_distr;
  for(int ijack=0;ijack<Njacks;ijack++) m_etas_phys_distr.distr.push_back( m_etas + GM()*m_etas_err/sqrt(Njacks-1.0));
  for(int ijack=0;ijack<Njacks;ijack++) m_phi_phys_distr.distr.push_back( m_phi+ GM()*m_phi_err/sqrt(Njacks-1.0));
  distr_t m2_etas_distr= m_etas_phys_distr*m_etas_phys_distr;
  distr_t m2_phi_distr = m_phi_phys_distr*m_phi_phys_distr;
  
  
  //push_back MV_strange and ZV_strange
  
  //Mpi_OS_fit_strange.distr_list.push_back(Mpi_OS);
  //L
  MV_fit_strange_OS_L.distr_list.push_back(MV_strange_OS_L);
  ZV_fit_strange_OS_L.distr_list.push_back(Zv*Zv*ZV_strange_OS_L);
  //M
  MV_fit_strange_OS_M.distr_list.push_back(MV_strange_OS_M);
  ZV_fit_strange_OS_M.distr_list.push_back(Zv*Zv*ZV_strange_OS_M);
  //H
  /*
  MV_fit_strange_OS_H.distr_list.push_back(MV_strange_OS_H);
  ZV_fit_strange_OS_H.distr_list.push_back(Zv*Zv*ZV_strange_OS_H);
  */
 


 
  int Tdata_min= 8;
  int Tdata_max = Corr.Nt/2.0 -2;
  int Tdata_fit = (Tdata_min + Tdata_max)/2;
  
  //compute kernel distribution
  distr_t_list Kernel_distr_No_ELM= distr_t_list::f_of_distr(K,a_distr, Upper_Limit_Time_Integral_strange+1);
  //tm
  //L
  distr_t_list Kernel_distr_L = distr_t_list::f_of_distr(K,MV_strange_L/m_phi, Upper_Limit_Time_Integral_strange+1);
  //M
  distr_t_list Kernel_distr_M = distr_t_list::f_of_distr(K,MV_strange_M/m_phi, Upper_Limit_Time_Integral_strange+1);
  //H
  //distr_t_list Kernel_distr_H = distr_t_list::f_of_distr(K,MV_strange_H/m_phi, Upper_Limit_Time_Integral_strange+1);
  
  //OS
  //L
  distr_t_list Kernel_OS_distr_L = distr_t_list::f_of_distr(K, MV_strange_OS_L/m_phi, Upper_Limit_Time_Integral_strange +1);
  //M
  distr_t_list Kernel_OS_distr_M = distr_t_list::f_of_distr(K, MV_strange_OS_M/m_phi, Upper_Limit_Time_Integral_strange +1);
  //H
  //distr_t_list Kernel_OS_distr_H = distr_t_list::f_of_distr(K, MV_strange_OS_H/m_phi, Upper_Limit_Time_Integral_strange +1);
  
  //compute exp(-Mv*t) distribution
  //tm
  //L
  distr_t_list exp_MVs_L = distr_t_list::f_of_distr(exp_MV, MV_strange_L, Upper_Limit_Time_Integral_strange+1);
  //M
  distr_t_list exp_MVs_M = distr_t_list::f_of_distr(exp_MV, MV_strange_M, Upper_Limit_Time_Integral_strange+1);
  //H
  //distr_t_list exp_MVs_H = distr_t_list::f_of_distr(exp_MV, MV_strange_H, Upper_Limit_Time_Integral_strange+1);

  //OS
  //L
  distr_t_list exp_OS_MVs_L = distr_t_list::f_of_distr(exp_MV, MV_strange_OS_L, Upper_Limit_Time_Integral_strange+1);
  //M
  distr_t_list exp_OS_MVs_M = distr_t_list::f_of_distr(exp_MV, MV_strange_OS_M, Upper_Limit_Time_Integral_strange+1);
  //H
  //distr_t_list exp_OS_MVs_H = distr_t_list::f_of_distr(exp_MV, MV_strange_OS_H, Upper_Limit_Time_Integral_strange+1);

  
  //#######################################################################################################################################à
  //Print single-exponential prediction to file
  //tm
  //L
  Print_To_File({}, {(exp_MVs_L*(ZV_strange_L/(2.0*MV_strange_L))).ave(), (exp_MVs_L*(ZV_strange_L/(2.0*MV_strange_L))).err(), (Za*Za*exp_MVs_L*(ZV_strange_L/(2.0*MV_strange_L))).ave(), (Za*Za*exp_MVs_L*(ZV_strange_L/(2.0*MV_strange_L))).err() }, "../data/gm2/strange/tm/corr_gsd_sum_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare   renormalized");
  //M
  Print_To_File({}, {(exp_MVs_M*(ZV_strange_M/(2.0*MV_strange_M))).ave(), (exp_MVs_M*(ZV_strange_M/(2.0*MV_strange_M))).err(), (Za_heavy*Za_heavy*exp_MVs_M*(ZV_strange_M/(2.0*MV_strange_M))).ave(), (Za_heavy*Za_heavy*exp_MVs_M*(ZV_strange_M/(2.0*MV_strange_M))).err() }, "../data/gm2/strange/tm/corr_gsd_sum_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#t   bare   renormalized");
  //H
  //Print_To_File({}, {(exp_MVs_H*(ZV_strange_H/(2.0*MV_strange_H))).ave(), (exp_MVs_H*(ZV_strange_H/(2.0*MV_strange_H))).err(), (Za*Za*exp_MVs_H*(ZV_strange_H/(2.0*MV_strange_H))).ave(), (Za*Za*exp_MVs_H*(ZV_strange_H/(2.0*MV_strange_H))).err() }, "../data/gm2/strange/tm/corr_gsd_sum_"+V_strange_1_M.Tag[i_ens]+"_H.dat.t", "", "#t   bare   renormalized");


   
  //OS
  //L
  Print_To_File({}, {(exp_OS_MVs_L*(ZV_strange_OS_L/(2.0*MV_strange_OS_L))).ave(), (exp_OS_MVs_L*(ZV_strange_OS_L/(2.0*MV_strange_OS_L))).err(), (Zv*Zv*exp_OS_MVs_L*(ZV_strange_OS_L/(2.0*MV_strange_OS_L))).ave(), (Zv*Zv*exp_OS_MVs_L*(ZV_strange_OS_L/(2.0*MV_strange_OS_L))).err()}, "../data/gm2/strange/OS/corr_gsd_sum_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare   renormalized");
  //M
  Print_To_File({}, {(exp_OS_MVs_M*(ZV_strange_OS_M/(2.0*MV_strange_OS_M))).ave(), (exp_OS_MVs_M*(ZV_strange_OS_M/(2.0*MV_strange_OS_M))).err(), (Zv_heavy*Zv_heavy*exp_OS_MVs_M*(ZV_strange_OS_M/(2.0*MV_strange_OS_M))).ave(), (Zv_heavy*Zv_heavy*exp_OS_MVs_M*(ZV_strange_OS_M/(2.0*MV_strange_OS_M))).err()}, "../data/gm2/strange/OS/corr_gsd_sum_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#t   bare   renormalized");
  //H
  //Print_To_File({}, {(exp_OS_MVs_H*(ZV_strange_OS_H/(2.0*MV_strange_OS_H))).ave(), (exp_OS_MVs_H*(ZV_strange_OS_H/(2.0*MV_strange_OS_H))).err(), (Zv*Zv*exp_OS_MVs_H*(ZV_strange_OS_H/(2.0*MV_strange_OS_H))).ave(), (Zv*Zv*exp_OS_MVs_H*(ZV_strange_OS_H/(2.0*MV_strange_OS_H))).err()}, "../data/gm2/strange/OS/corr_gsd_sum_"+V_strange_1_H.Tag[i_ens]+"_H.dat.t", "", "#t   bare   renormalized");
  //########################################################################################################################################
  
  distr_t_list agm2_distr_Tdata_L(UseJack), agm2_OS_distr_Tdata_L(UseJack);
  distr_t_list agm2_distr_Tdata_No_ELM_L(UseJack), agm2_OS_distr_Tdata_No_ELM_L(UseJack);
  distr_t_list agm2_distr_Tdata_M(UseJack), agm2_OS_distr_Tdata_M(UseJack);
  distr_t_list agm2_distr_Tdata_No_ELM_M(UseJack), agm2_OS_distr_Tdata_No_ELM_M(UseJack);
  distr_t_list agm2_distr_Tdata_H(UseJack), agm2_OS_distr_Tdata_H(UseJack);
  Vfloat Tdata_vec;
  bool Find_Tdata_fit=false;

   
  for(int Tdata=Tdata_min;Tdata<Tdata_max;Tdata++) {
    //compute 4\pia^2 using lattice data up to Tdata (included)
   
    distr_t agm2_L(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_No_ELM_L(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_OS_L(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_OS_No_ELM_L(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_M(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_No_ELM_M(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_OS_M(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_OS_No_ELM_M(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_H(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    distr_t agm2_OS_H(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
    
    for(int t=1;t<=Upper_Limit_Time_Integral_strange;t++) {
      if(t<=Tdata) {
	//L
	agm2_L = agm2_L + 4.0*pow(alpha,2)*V_strange_distr_L_pert_sub.distr_list[t]*Kernel_distr_L.distr_list[t];
	agm2_No_ELM_L = agm2_No_ELM_L + 4.0*pow(alpha,2)*V_strange_distr_L_pert_sub.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	agm2_OS_L = agm2_OS_L + 4.0*pow(alpha,2)*V_strange_OS_distr_L_pert_sub.distr_list[t]*Kernel_OS_distr_L.distr_list[t];
	agm2_OS_No_ELM_L = agm2_OS_No_ELM_L + 4.0*pow(alpha,2)*V_strange_OS_distr_L_pert_sub.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	//M
	agm2_M = agm2_M + 4.0*pow(alpha,2)*V_strange_distr_M_pert_sub.distr_list[t]*Kernel_distr_M.distr_list[t];
	agm2_No_ELM_M = agm2_No_ELM_M + 4.0*pow(alpha,2)*V_strange_distr_M_pert_sub.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	agm2_OS_M = agm2_OS_M + 4.0*pow(alpha,2)*V_strange_OS_distr_M_pert_sub.distr_list[t]*Kernel_OS_distr_M.distr_list[t];
	agm2_OS_No_ELM_M = agm2_OS_No_ELM_M + 4.0*pow(alpha,2)*V_strange_OS_distr_M_pert_sub.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	//H
	/*
	agm2_H = agm2_H + 4.0*pow(alpha,2)*V_strange_distr_H.distr_list[t]*Kernel_distr_H.distr_list[t];
	agm2_OS_H = agm2_OS_H + 4.0*pow(alpha,2)*V_strange_OS_distr_H.distr_list[t]*Kernel_OS_distr_H.distr_list[t];
	*/
	
      }
      else {
	//L
	agm2_L= agm2_L + 4.0*pow(alpha,2)*(ZV_strange_L/(2.0*MV_strange_L))*exp_MVs_L.distr_list[t]*Kernel_distr_L.distr_list[t];
	agm2_No_ELM_L= agm2_No_ELM_L + 4.0*pow(alpha,2)*(ZV_strange_L/(2.0*MV_strange_L))*exp_MVs_L.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	agm2_OS_L= agm2_OS_L + 4.0*pow(alpha,2)*(ZV_strange_OS_L/(2.0*MV_strange_OS_L))*exp_OS_MVs_L.distr_list[t]*Kernel_OS_distr_L.distr_list[t];
	agm2_OS_No_ELM_L= agm2_OS_No_ELM_L + 4.0*pow(alpha,2)*(ZV_strange_OS_L/(2.0*MV_strange_OS_L))*exp_OS_MVs_L.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	//M
	agm2_M= agm2_M + 4.0*pow(alpha,2)*(ZV_strange_M/(2.0*MV_strange_M))*exp_MVs_M.distr_list[t]*Kernel_distr_M.distr_list[t];
	agm2_No_ELM_M= agm2_No_ELM_M + 4.0*pow(alpha,2)*(ZV_strange_M/(2.0*MV_strange_M))*exp_MVs_M.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	agm2_OS_M= agm2_OS_M + 4.0*pow(alpha,2)*(ZV_strange_OS_M/(2.0*MV_strange_OS_M))*exp_OS_MVs_M.distr_list[t]*Kernel_OS_distr_M.distr_list[t];
	agm2_OS_No_ELM_M= agm2_OS_No_ELM_M + 4.0*pow(alpha,2)*(ZV_strange_OS_M/(2.0*MV_strange_OS_M))*exp_OS_MVs_M.distr_list[t]*Kernel_distr_No_ELM.distr_list[t];
	//H
	/*
	agm2_H= agm2_H + 4.0*pow(alpha,2)*(ZV_strange_H/(2.0*MV_strange_H))*exp_MVs_H.distr_list[t]*Kernel_distr_H.distr_list[t];
	agm2_OS_H= agm2_OS_H + 4.0*pow(alpha,2)*(ZV_strange_OS_H/(2.0*MV_strange_OS_H))*exp_OS_MVs_H.distr_list[t]*Kernel_OS_distr_H.distr_list[t];
	*/
      }
    }
    
    Tdata_vec.push_back((double)Tdata);
    //L
    agm2_distr_Tdata_L.distr_list.push_back(Za*Za*agm2_L);
    agm2_OS_distr_Tdata_L.distr_list.push_back(Zv*Zv*agm2_OS_L);
    agm2_distr_Tdata_No_ELM_L.distr_list.push_back(Za*Za*agm2_No_ELM_L);
    agm2_OS_distr_Tdata_No_ELM_L.distr_list.push_back(Zv*Zv*agm2_OS_No_ELM_L);
    //M
    agm2_distr_Tdata_M.distr_list.push_back(Za_heavy*Za_heavy*agm2_M);
    agm2_OS_distr_Tdata_M.distr_list.push_back(Zv_heavy*Zv_heavy*agm2_OS_M);
    agm2_distr_Tdata_No_ELM_M.distr_list.push_back(Za_heavy*Za_heavy*agm2_No_ELM_M);
    agm2_OS_distr_Tdata_No_ELM_M.distr_list.push_back(Zv_heavy*Zv_heavy*agm2_OS_No_ELM_M);
    //H
    /*
    agm2_distr_Tdata_H.distr_list.push_back(Za*Za*agm2_H);
    agm2_OS_distr_Tdata_H.distr_list.push_back(Zv*Zv*agm2_OS_H);
    */


    
    
    if(Tdata==Tdata_fit){
      //push back L
      agm2_strange_L.distr_list.push_back(Za*Za*agm2_L);
      agm2_strange_No_ELM_L.distr_list.push_back(Za*Za*agm2_No_ELM_L);
      agm2_strange_OS_L.distr_list.push_back(Zv*Zv*agm2_OS_L);
      agm2_strange_OS_No_ELM_L.distr_list.push_back(Zv*Zv*agm2_OS_No_ELM_L);
      //push back M
      agm2_strange_M.distr_list.push_back(Za_heavy*Za_heavy*agm2_M);
      agm2_strange_No_ELM_M.distr_list.push_back(Za_heavy*Za_heavy*agm2_No_ELM_M);
      agm2_strange_OS_M.distr_list.push_back(Zv_heavy*Zv_heavy*agm2_OS_M);
      agm2_strange_OS_No_ELM_M.distr_list.push_back(Zv_heavy*Zv_heavy*agm2_OS_No_ELM_M);
      //push back H
      /*
      agm2_strange_H.distr_list.push_back(Za*Za*agm2_H);
      agm2_strange_OS_H.distr_list.push_back(Zv*Zv*agm2_OS_H);
      */
      //extrapolate to the physical kaon point
      vector<distr_t> agm2s_strange({Za*Za*agm2_L, Za_heavy*Za_heavy*agm2_M});
      vector<distr_t> agm2s_strange_OS({ Zv*Zv*agm2_OS_L, Zv_heavy*Zv_heavy*agm2_OS_M});
      vector<distr_t> agm2s_strange_No_ELM({ Za*Za*agm2_No_ELM_L, Za_heavy*Za_heavy*agm2_No_ELM_M});
      vector<distr_t> agm2s_strange_OS_No_ELM({ Zv*Zv*agm2_OS_No_ELM_L, Zv_heavy*Zv_heavy*agm2_OS_No_ELM_M});
      agm2_strange_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange, Mphi2, Use_phi_etas?m2_phi_distr:m2_etas_distr,  "../data/gm2/strange", "agm2_ELM_"+V_strange_1_L.Tag[i_ens], UseJack));
      agm2_strange_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_OS, Mphi2, Use_phi_etas?m2_phi_distr:m2_etas_distr,  "../data/gm2/strange", "agm2_OS_ELM_"+V_strange_1_L.Tag[i_ens], UseJack));
      agm2_strange_No_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_No_ELM, Mphi2, Use_phi_etas?m2_phi_distr:m2_etas_distr,  "../data/gm2/strange", "agm2_"+V_strange_1_L.Tag[i_ens], UseJack));
      agm2_strange_OS_No_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_OS_No_ELM, Mphi2, Use_phi_etas?m2_phi_distr:m2_etas_distr,  "../data/gm2/strange", "agm2_OS_"+V_strange_1_L.Tag[i_ens], UseJack));
      
      Find_Tdata_fit=true;
    }
  }

  if(!Find_Tdata_fit) crash("Cannot find Tdata fit value: "+to_string(Tdata_fit));


  
  //#######################  INTERMEDIATE AND SHORT-DISTANCE WINDOW ###################################

    
  //############################   TWISTED MASS ######################################################

  //L
  distr_t agm2_W_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
  distr_t agm2_SD_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
  distr_t agm2_W_ELM_L(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
  distr_t agm2_SD_ELM_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default
  //M
  distr_t agm2_W_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
  distr_t agm2_SD_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
  distr_t agm2_W_ELM_M(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
  distr_t agm2_SD_ELM_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default
  //H
  distr_t agm2_W_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
  distr_t agm2_SD_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
  distr_t agm2_W_ELM_H(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
  distr_t agm2_SD_ELM_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default


  //#################################################################################################


  
  distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);
  //L
  distr_t_list Ker_ELM_tm_L = distr_t_list::f_of_distr(K, MV_strange_L/m_phi, Corr.Nt/2);
  distr_t_list Ker_ELM_OS_L = distr_t_list::f_of_distr(K, MV_strange_OS_L/m_phi, Corr.Nt/2);
  //M
  distr_t_list Ker_ELM_tm_M = distr_t_list::f_of_distr(K, MV_strange_M/m_phi, Corr.Nt/2);
  distr_t_list Ker_ELM_OS_M = distr_t_list::f_of_distr(K, MV_strange_OS_M/m_phi, Corr.Nt/2);
  //H
  /*
  distr_t_list Ker_ELM_tm_H = distr_t_list::f_of_distr(K, MV_strange_H/m_phi, Corr.Nt/2);
  distr_t_list Ker_ELM_OS_H = distr_t_list::f_of_distr(K, MV_strange_OS_H/m_phi, Corr.Nt/2);
  */
    
  //define lambdas for the theta func
  auto th0 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t0)/Delta));};
  auto th1 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t1)/Delta));};

  
  for(int t=1; t< Corr.Nt/2; t++) {
    //L
    agm2_W_L = agm2_W_L + 4.0*pow(alpha,2)*Za*Za*V_strange_distr_L.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_L = agm2_SD_L + 4.0*pow(alpha,2)*Za*Za*(V_strange_distr_L_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_L = agm2_W_ELM_L + 4.0*pow(alpha,2)*Za*Za*V_strange_distr_L.distr_list[t]*Ker_ELM_tm_L.distr_list[t]*( distr_t::f_of_distr(th0, t*MV_strange_L/m_phi) - distr_t::f_of_distr(th1, t*MV_strange_L/m_phi));
    agm2_SD_ELM_L = agm2_SD_ELM_L + 4.0*pow(alpha,2)*Za*Za*(V_strange_distr_L_pert_sub.distr_list[t])*Ker_ELM_tm_L.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*MV_strange_L/m_phi));
    //M
    agm2_W_M = agm2_W_M + 4.0*pow(alpha,2)*Za_heavy*Za_heavy*V_strange_distr_M.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_M = agm2_SD_M + 4.0*pow(alpha,2)*Za_heavy*Za_heavy*(V_strange_distr_M_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_M = agm2_W_ELM_M + 4.0*pow(alpha,2)*Za_heavy*Za_heavy*V_strange_distr_M.distr_list[t]*Ker_ELM_tm_M.distr_list[t]*( distr_t::f_of_distr(th0, t*MV_strange_M/m_phi) - distr_t::f_of_distr(th1, t*MV_strange_M/m_phi));
    agm2_SD_ELM_M = agm2_SD_ELM_M + 4.0*pow(alpha,2)*Za_heavy*Za_heavy*(V_strange_distr_M_pert_sub.distr_list[t])*Ker_ELM_tm_M.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*MV_strange_M/m_phi));
    //H
    /*
    agm2_W_H = agm2_W_H + 4.0*pow(alpha,2)*Za*Za*V_strange_distr_H.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_H = agm2_SD_H + 4.0*pow(alpha,2)*Za*Za*(V_strange_distr_H.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_H = agm2_W_ELM_H + 4.0*pow(alpha,2)*Za*Za*V_strange_distr_H.distr_list[t]*Ker_ELM_tm_H.distr_list[t]*( distr_t::f_of_distr(th0, t*MV_strange_H/m_phi) - distr_t::f_of_distr(th1, t*MV_strange_H/m_phi));
    agm2_SD_ELM_H = agm2_SD_ELM_H + 4.0*pow(alpha,2)*Za*Za*(V_strange_distr_H.distr_list[t])*Ker_ELM_tm_H.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*MV_strange_H/m_phi));
    */

  }
  
  //push_back the result

  //L
  agm2_strange_W_L.distr_list.push_back(agm2_W_L);
  agm2_strange_SD_L.distr_list.push_back(agm2_SD_L);
  agm2_strange_W_ELM_L.distr_list.push_back(agm2_W_ELM_L);
  agm2_strange_SD_ELM_L.distr_list.push_back(agm2_SD_ELM_L);
  //M
  agm2_strange_W_M.distr_list.push_back(agm2_W_M);
  agm2_strange_SD_M.distr_list.push_back(agm2_SD_M);
  agm2_strange_W_ELM_M.distr_list.push_back(agm2_W_ELM_M);
  agm2_strange_SD_ELM_M.distr_list.push_back(agm2_SD_ELM_M);
  //H
  /*
  agm2_strange_W_H.distr_list.push_back(agm2_W_H);
  agm2_strange_SD_H.distr_list.push_back(agm2_SD_H);
  agm2_strange_W_ELM_H.distr_list.push_back(agm2_W_ELM_H);
  agm2_strange_SD_ELM_H.distr_list.push_back(agm2_SD_ELM_H);
  */


  //extrapolate the result to the physical kaon point

  vector<distr_t> agm2s_strange_W({agm2_W_L, agm2_W_M});
  vector<distr_t> agm2s_strange_SD({agm2_SD_L, agm2_SD_M});
  vector<distr_t> agm2s_strange_W_ELM({agm2_W_ELM_L, agm2_W_ELM_M});
  vector<distr_t> agm2s_strange_SD_ELM({agm2_SD_ELM_L, agm2_SD_ELM_M});

  
  agm2_strange_W_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_W, Mphi2, Use_phi_etas?m2_phi_distr:m2_etas_distr, "../data/gm2/strange", "agm2_W_"+V_strange_1_L.Tag[i_ens], UseJack));
  agm2_strange_SD_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_SD, Mphi2, Use_phi_etas?m2_phi_distr:m2_etas_distr, "../data/gm2/strange", "agm2_SD_"+V_strange_1_L.Tag[i_ens], UseJack));
  agm2_strange_W_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_W_ELM, Mphi2,Use_phi_etas?m2_phi_distr:m2_etas_distr, "../data/gm2/strange", "agm2_W_ELM_"+V_strange_1_L.Tag[i_ens], UseJack));
  agm2_strange_SD_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_SD_ELM, Mphi2, Use_phi_etas?m2_phi_distr:m2_etas_distr,  "../data/gm2/strange", "agm2_SD_ELM_"+V_strange_1_L.Tag[i_ens], UseJack));
 
  
  

  //####################################################################################################



  // ######################################## OSTERWALDER-SEILER #######################################

  //L
  distr_t agm2_W_OS_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
  distr_t agm2_SD_OS_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
  distr_t agm2_W_ELM_OS_L(UseJack, UseJack?Njacks:Nboots); //constructor sets  to zero by default
  distr_t agm2_SD_ELM_OS_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
  //M
  distr_t agm2_W_OS_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
  distr_t agm2_SD_OS_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
  distr_t agm2_W_ELM_OS_M(UseJack, UseJack?Njacks:Nboots); //constructor sets  to zero by default
  distr_t agm2_SD_ELM_OS_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
  //H
  distr_t agm2_W_OS_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
  distr_t agm2_SD_OS_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default
  distr_t agm2_W_ELM_OS_H(UseJack, UseJack?Njacks:Nboots); //constructor sets  to zero by default
  distr_t agm2_SD_ELM_OS_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets  to zero by default


  //#################################################################################################

  for(int t=1; t< Corr.Nt/2; t++) {
    //L
    agm2_W_OS_L = agm2_W_OS_L + 4.0*pow(alpha,2)*Zv*Zv*V_strange_OS_distr_L.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_OS_L = agm2_SD_OS_L + 4.0*pow(alpha,2)*Zv*Zv*(V_strange_OS_distr_L_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_OS_L = agm2_W_ELM_OS_L + 4.0*pow(alpha,2)*Zv*Zv*V_strange_OS_distr_L.distr_list[t]*Ker_ELM_OS_L.distr_list[t]*( distr_t::f_of_distr(th0, t*MV_strange_OS_L/m_phi) - distr_t::f_of_distr(th1, t*MV_strange_OS_L/m_phi));
    agm2_SD_ELM_OS_L = agm2_SD_ELM_OS_L + 4.0*pow(alpha,2)*Zv*Zv*V_strange_OS_distr_L_pert_sub.distr_list[t]*Ker_ELM_OS_L.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*MV_strange_OS_L/m_phi));
    //M
    agm2_W_OS_M = agm2_W_OS_M + 4.0*pow(alpha,2)*Zv_heavy*Zv_heavy*V_strange_OS_distr_M.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_OS_M = agm2_SD_OS_M + 4.0*pow(alpha,2)*Zv_heavy*Zv_heavy*(V_strange_OS_distr_M_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_OS_M = agm2_W_ELM_OS_M + 4.0*pow(alpha,2)*Zv_heavy*Zv_heavy*V_strange_OS_distr_M.distr_list[t]*Ker_ELM_OS_M.distr_list[t]*( distr_t::f_of_distr(th0, t*MV_strange_OS_M/m_phi) - distr_t::f_of_distr(th1, t*MV_strange_OS_M/m_phi));
    agm2_SD_ELM_OS_M = agm2_SD_ELM_OS_M + 4.0*pow(alpha,2)*Zv_heavy*Zv_heavy*V_strange_OS_distr_M_pert_sub.distr_list[t]*Ker_ELM_OS_M.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*MV_strange_OS_M/m_phi));
    //H
    /*
    agm2_W_OS_H = agm2_W_OS_H + 4.0*pow(alpha,2)*Zv*Zv*V_strange_OS_distr_H.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_OS_H = agm2_SD_OS_H + 4.0*pow(alpha,2)*Zv*Zv*(V_strange_OS_distr_H.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_OS_H = agm2_W_ELM_OS_H + 4.0*pow(alpha,2)*Zv*Zv*V_strange_OS_distr_H.distr_list[t]*Ker_ELM_OS_H.distr_list[t]*( distr_t::f_of_distr(th0, t*MV_strange_OS_H/m_phi) - distr_t::f_of_distr(th1, t*MV_strange_OS_H/m_phi));
    agm2_SD_ELM_OS_H = agm2_SD_ELM_OS_H + 4.0*pow(alpha,2)*Zv*Zv*V_strange_OS_distr_H.distr_list[t]*Ker_ELM_OS_H.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*MV_strange_OS_H/m_phi));
    */
  }
  
  //push_back the result

  //L
  agm2_strange_W_OS_L.distr_list.push_back(agm2_W_OS_L);
  agm2_strange_SD_OS_L.distr_list.push_back(agm2_SD_OS_L);
  agm2_strange_W_ELM_OS_L.distr_list.push_back(agm2_W_ELM_OS_L);
  agm2_strange_SD_ELM_OS_L.distr_list.push_back(agm2_SD_ELM_OS_L);
  //M
  agm2_strange_W_OS_M.distr_list.push_back(agm2_W_OS_M);
  agm2_strange_SD_OS_M.distr_list.push_back(agm2_SD_OS_M);
  agm2_strange_W_ELM_OS_M.distr_list.push_back(agm2_W_ELM_OS_M);
  agm2_strange_SD_ELM_OS_M.distr_list.push_back(agm2_SD_ELM_OS_M);
  //H
  /*
  agm2_strange_W_OS_H.distr_list.push_back(agm2_W_OS_H);
  agm2_strange_SD_OS_H.distr_list.push_back(agm2_SD_OS_H);
  agm2_strange_W_ELM_OS_H.distr_list.push_back(agm2_W_ELM_OS_H);
  agm2_strange_SD_ELM_OS_H.distr_list.push_back(agm2_SD_ELM_OS_H);
  */


  //extrapolate the result to the physical kaon point
  vector<distr_t> agm2s_strange_W_OS({agm2_W_OS_L, agm2_W_OS_M});
  vector<distr_t> agm2s_strange_SD_OS({agm2_SD_OS_L, agm2_SD_OS_M});
  vector<distr_t> agm2s_strange_W_ELM_OS({agm2_W_ELM_OS_L, agm2_W_ELM_OS_M});
  vector<distr_t> agm2s_strange_SD_ELM_OS({agm2_SD_ELM_OS_L, agm2_SD_ELM_OS_M});

  
  agm2_strange_W_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_W_OS, Mphi2, Use_phi_etas?m2_phi_distr:m2_etas_distr, "../data/gm2/strange", "agm2_W_OS_"+V_strange_1_L.Tag[i_ens], UseJack));
  agm2_strange_SD_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_SD_OS, Mphi2, Use_phi_etas?m2_phi_distr:m2_etas_distr, "../data/gm2/strange", "agm2_SD_OS_"+V_strange_1_L.Tag[i_ens], UseJack));
  agm2_strange_W_ELM_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_W_ELM_OS, Mphi2, Use_phi_etas?m2_phi_distr:m2_etas_distr, "../data/gm2/strange", "agm2_W_ELM_OS_"+V_strange_1_L.Tag[i_ens], UseJack));
  agm2_strange_SD_ELM_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_strange_SD_ELM_OS, Mphi2, Use_phi_etas?m2_phi_distr:m2_etas_distr, "../data/gm2/strange", "agm2_SD_ELM_OS_"+V_strange_1_L.Tag[i_ens], UseJack));

  //####################################################################################################

  
  
  //print to file
  //L
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_L.ave(), agm2_distr_Tdata_L.err()}, "../data/gm2/strange/tm/agm2_Tdata_ELM_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_L.ave(), agm2_OS_distr_Tdata_L.err()}, "../data/gm2/strange/OS/agm2_Tdata_ELM_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_No_ELM_L.ave(), agm2_distr_Tdata_No_ELM_L.err()}, "../data/gm2/strange/tm/agm2_Tdata_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_No_ELM_L.ave(), agm2_OS_distr_Tdata_No_ELM_L.err()}, "../data/gm2/strange/OS/agm2_Tdata_"+V_strange_1_L.Tag[i_ens]+"_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
  //M
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_M.ave(), agm2_distr_Tdata_M.err()}, "../data/gm2/strange/tm/agm2_Tdata_ELM_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_M.ave(), agm2_OS_distr_Tdata_M.err()}, "../data/gm2/strange/OS/agm2_Tdata_ELM_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_No_ELM_M.ave(), agm2_distr_Tdata_No_ELM_M.err()}, "../data/gm2/strange/tm/agm2_Tdata_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_No_ELM_M.ave(), agm2_OS_distr_Tdata_No_ELM_M.err()}, "../data/gm2/strange/OS/agm2_Tdata_"+V_strange_1_M.Tag[i_ens]+"_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
  //H
  /*
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_H.ave(), agm2_distr_Tdata_H.err()}, "../data/gm2/strange/tm/agm2_Tdata_"+V_strange_1_H.Tag[i_ens]+"_ELM_H.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_H.ave(), agm2_OS_distr_Tdata_H.err()}, "../data/gm2/strange/OS/agm2_Tdata_"+V_strange_1_H.Tag[i_ens]+"_ELM_H.dat.t", "", "#id  Tdata   ag2m agm2_err");
  */

  
  }

  cout<<"strange quark correlator analyzed!"<<endl;

  //charm
  channel="c";
  for(int i_ens=0;i_ens<Nens_charm;i_ens++) { 
    
  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = V_charm_1_L.nrows[i_ens];


  

  //resample lattice spacing
  distr_t a_distr(UseJack), Za_L(UseJack), Zv_L(UseJack);
  LatticeInfo L_info;
  L_info.LatInfo_new_ens(V_charm_1_L.Tag[i_ens]);
  //generate jackknife sample of input parameters
  
  if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
      Za_L.distr.push_back(  L_info.Za + GM()*L_info.Za_err*(1.0/sqrt(Njacks-1.0)));
      Zv_L.distr.push_back(  L_info.Zv + GM()*L_info.Zv_err*(1.0/sqrt(Njacks-1.0)));
    }
  }
  else {
    for (int iboot=0; iboot<Nboots;iboot++) {
      Za_L.distr.push_back(  L_info.Za + GM()*L_info.Za_err);
      Zv_L.distr.push_back(  L_info.Zv + GM()*L_info.Zv_err);
    }
  }

  distr_t Za_M = Za_L;
  distr_t Za_H = Za_L;
  distr_t Zv_M = Zv_L;
  distr_t Zv_H = Zv_H;

   //push_back lattice info
  L_charm_list.push_back(L_info.L);
  a_charm_list.push_back(L_info.a);
  ml_charm_list.push_back(L_info.ml);


  //Use Za and Zv from strange run if required
  //if(Use_Za_Zv_from_strange_run) { Zv= Zv_fit_strange[i_ens]; Za= Za_fit_strange[i_ens];}

  int Tmin_P5P5;
  int Tmax_P5P5;
  int Tmin_VV;
  int Tmax_VV;

  
  //set time intervals for pseudoscalar obs
  if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
    a_distr = a_C;
    if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Corr.Tmin=23; Corr.Tmax=39;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
    a_distr = a_B;
    if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Corr.Tmin=19; Corr.Tmax=36;}
    else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Corr.Tmin=21;Corr.Tmax=34;}
    else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Corr.Tmin=16; Corr.Tmax= 27;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
    a_distr = a_A;
    if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Corr.Tmin=18; Corr.Tmax=28;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Corr.Tmin=16; Corr.Tmax=22;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Corr.Tmin=15; Corr.Tmax=22;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Corr.Tmin=14; Corr.Tmax=25;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else crash("Ensemble tag not valid");

  a_distr_list_charm.distr_list.push_back(a_distr);

  //set Tmin_P5P5 and Tmax_P5P5 to the values Corr.Tmin and Corr.Tmax
  Tmin_P5P5 = Corr.Tmin;
  Tmax_P5P5 = Corr.Tmax;

  
  //set time intervals for vector obs
  if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
    if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_VV=36; Tmax_VV=66;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
    if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_VV=31; Tmax_VV=43;}
    else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_VV=32;Tmax_VV=42;}
    else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_VV=28; Tmax_VV= 40;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
    if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_VV=27; Tmax_VV=42;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_VV=19; Tmax_VV=24;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_VV=20; Tmax_VV=24;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_VV=24;Tmax_VV=31;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }

  else crash("Ensemble tag not valid");

  

  //set Tmin and Tmax for the eta_C
  int Tmin_etaC, Tmax_etaC;
  if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
    if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_etaC=33; Tmax_etaC=60;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
    if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_etaC=33; Tmax_etaC=52;}
    else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_etaC=32;Tmax_etaC=46;}
    else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_etaC=40; Tmax_etaC= 62;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
    if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_etaC=28; Tmax_etaC=46;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_etaC=19; Tmax_etaC=24;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_etaC=21; Tmax_etaC=24;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_etaC=22;Tmax_etaC=31;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else crash("Ensemble tag not valid");
  
  

  //tm
  distr_t_list Mpi_distr;
  distr_t Mpi;
  //L
  distr_t_list  V_charm_1_distr_L, V_charm_2_distr_L, V_charm_3_distr_L;
  distr_t_list  V_charm_distr_L, MV_charm_distr_L, ZV_charm_distr_L;
  distr_t_list   MD_distr_L, FD_distr_L, M_etaC_distr_L, M_etaC_distr_OS_L;
  distr_t MV_charm_L , ZV_charm_L, MD_L, FD_L, M_etaC_L;
  //M
  distr_t_list  V_charm_1_distr_M, V_charm_2_distr_M, V_charm_3_distr_M;
  distr_t_list  V_charm_distr_M, MV_charm_distr_M, ZV_charm_distr_M;
  distr_t_list  MD_distr_M, FD_distr_M, M_etaC_distr_M, M_etaC_distr_OS_M;
  distr_t MV_charm_M, ZV_charm_M, MD_M, FD_M, M_etaC_M;
  //H
  distr_t_list  V_charm_1_distr_H, V_charm_2_distr_H, V_charm_3_distr_H;
  distr_t_list  V_charm_distr_H, MV_charm_distr_H, ZV_charm_distr_H;
  distr_t_list   MD_distr_H, FD_distr_H, M_etaC_distr_H, M_etaC_distr_OS_H;
  distr_t MV_charm_H , ZV_charm_H, MD_H, FD_H, M_etaC_H;


  //OS
  distr_t_list Mpi_OS_distr;
  distr_t Mpi_OS;
  //L
  distr_t_list  V_charm_OS_1_distr_L, V_charm_OS_2_distr_L, V_charm_OS_3_distr_L;
  distr_t_list  V_charm_OS_distr_L, MV_charm_OS_distr_L, ZV_charm_OS_distr_L;
  distr_t MV_charm_OS_L , ZV_charm_OS_L;
  //M
  distr_t_list  V_charm_OS_1_distr_M, V_charm_OS_2_distr_M, V_charm_OS_3_distr_M;
  distr_t_list  V_charm_OS_distr_M, MV_charm_OS_distr_M, ZV_charm_OS_distr_M;
  distr_t MV_charm_OS_M, ZV_charm_OS_M;
  //H
  distr_t_list  V_charm_OS_1_distr_H, V_charm_OS_2_distr_H, V_charm_OS_3_distr_H;
  distr_t_list  V_charm_OS_distr_H, MV_charm_OS_distr_H, ZV_charm_OS_distr_H;
  distr_t MV_charm_OS_H , ZV_charm_OS_H;


  //observables to extract Zv and Za (Hadronic method)  (L)
  distr_t_list cbar_c_distr_L, cbar_c_OS_distr_L;
  distr_t_list overlap_P5P5_distr_L, overlap_P5P5_OS_distr_L;
  distr_t_list ratio_P5P5_overlap_OS_tm_L,  Zp_ov_Zs_distr_L;
  distr_t_list  A0P5_distr_L, A0P5_OS_distr_L;
  distr_t_list RA_L, RV_L;
  distr_t cbar_c_mass_L, cbar_c_OS_mass_L;
  distr_t Zp_ov_Zs_L, Zv_hadr_L, Za_hadr_L;

  //observables to extract Zv and Za (Hadronic method)  (M)
  distr_t_list cbar_c_distr_M, cbar_c_OS_distr_M;
  distr_t_list overlap_P5P5_distr_M, overlap_P5P5_OS_distr_M;
  distr_t_list ratio_P5P5_overlap_OS_tm_M,  Zp_ov_Zs_distr_M;
  distr_t_list  A0P5_distr_M, A0P5_OS_distr_M;
  distr_t_list RA_M, RV_M;
  distr_t cbar_c_mass_M, cbar_c_OS_mass_M;
  distr_t Zp_ov_Zs_M, Zv_hadr_M, Za_hadr_M;

  //observables to extract Zv and Za (Hadronic method)  (H)
  distr_t_list cbar_c_distr_H, cbar_c_OS_distr_H;
  distr_t_list overlap_P5P5_distr_H, overlap_P5P5_OS_distr_H;
  distr_t_list ratio_P5P5_overlap_OS_tm_H,  Zp_ov_Zs_distr_H;
  distr_t_list  A0P5_distr_H, A0P5_OS_distr_H;
  distr_t_list RA_H, RV_H;
  distr_t cbar_c_mass_H, cbar_c_OS_mass_H;
  distr_t Zp_ov_Zs_H, Zv_hadr_H, Za_hadr_H;


  
  //Analyze correlators

  
  //tm
  Mpi_distr = Corr.effective_mass_t(pt2_pion_charm.col(0)[i_ens], "../data/gm2/charm/tm/Mpi_"+V_charm_1_L.Tag[i_ens]+".dat");
  //L
  V_charm_1_distr_L = Corr.corr_t(V_charm_1_L.col(0)[i_ens], "../data/gm2/charm/tm/corr_1_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
  V_charm_2_distr_L = Corr.corr_t(V_charm_2_L.col(0)[i_ens], "../data/gm2/charm/tm/corr_2_"+V_charm_2_L.Tag[i_ens]+"_L.dat");
  V_charm_3_distr_L = Corr.corr_t(V_charm_3_L.col(0)[i_ens], "../data/gm2/charm/tm/corr_3_"+V_charm_3_L.Tag[i_ens]+"_L.dat");
  MD_distr_L = Corr.effective_mass_t(pt2_D_L.col(0)[i_ens], "../data/gm2/charm/tm/MD_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
  FD_distr_L = (L_info.ml + L_info.mc_L)*Corr.decay_constant_t(pt2_D_L.col(0)[i_ens], "../data/gm2/charm/tm/FD_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
  Corr.Tmin= Tmin_etaC;
  Corr.Tmax= Tmax_etaC;
  cbar_c_distr_L = Corr.corr_t(pt2_etaC_L.col(0)[i_ens], "");
  cbar_c_distr_M = Corr.corr_t(pt2_etaC_M.col(0)[i_ens], "");
  cbar_c_distr_H = Corr.corr_t(pt2_etaC_H.col(0)[i_ens], "");
  M_etaC_distr_L = Corr.effective_mass_t(cbar_c_distr_L, "../data/gm2/charm/tm/M_etaC_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
  M_etaC_distr_M = Corr.effective_mass_t(cbar_c_distr_M, "../data/gm2/charm/tm/M_etaC_"+V_charm_1_L.Tag[i_ens]+"_M.dat");
  M_etaC_distr_H = Corr.effective_mass_t(cbar_c_distr_H, "../data/gm2/charm/tm/M_etaC_"+V_charm_1_L.Tag[i_ens]+"_H.dat");
  overlap_P5P5_distr_L = Corr.residue_t(cbar_c_distr_L, "");
  overlap_P5P5_distr_M = Corr.residue_t(cbar_c_distr_M, "");
  overlap_P5P5_distr_H = Corr.residue_t(cbar_c_distr_H, "");
  cbar_c_OS_distr_L = Corr.corr_t(pt2_etaC_OS_L.col(0)[i_ens], "");
  cbar_c_OS_distr_M = Corr.corr_t(pt2_etaC_OS_M.col(0)[i_ens], "");
  cbar_c_OS_distr_H = Corr.corr_t(pt2_etaC_OS_H.col(0)[i_ens], "");
  M_etaC_distr_OS_L = Corr.effective_mass_t(cbar_c_OS_distr_L, "../data/gm2/charm/OS/M_etaC_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
  M_etaC_distr_OS_M = Corr.effective_mass_t(cbar_c_OS_distr_M, "../data/gm2/charm/OS/M_etaC_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
  M_etaC_distr_OS_H = Corr.effective_mass_t(cbar_c_OS_distr_H, "../data/gm2/charm/OS/M_etaC_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
  overlap_P5P5_OS_distr_L = Corr.residue_t(cbar_c_OS_distr_L, "");
  overlap_P5P5_OS_distr_M = Corr.residue_t(cbar_c_OS_distr_M, "");
  overlap_P5P5_OS_distr_H = Corr.residue_t(cbar_c_OS_distr_H, "");
  Corr.Tmin = Tmin_P5P5;
  Corr.Tmax = Tmax_P5P5;
  //M
  V_charm_1_distr_M = Corr.corr_t(V_charm_1_M.col(0)[i_ens], "../data/gm2/charm/tm/corr_1_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
  V_charm_2_distr_M = Corr.corr_t(V_charm_2_M.col(0)[i_ens], "../data/gm2/charm/tm/corr_2_"+V_charm_2_M.Tag[i_ens]+"_M.dat");
  V_charm_3_distr_M = Corr.corr_t(V_charm_3_M.col(0)[i_ens], "../data/gm2/charm/tm/corr_3_"+V_charm_3_M.Tag[i_ens]+"_M.dat");
  MD_distr_M = Corr.effective_mass_t(pt2_D_M.col(0)[i_ens], "../data/gm2/charm/tm/MD_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
  FD_distr_M = (L_info.ml + L_info.mc_M)*Corr.decay_constant_t(pt2_D_M.col(0)[i_ens], "../data/gm2/charm/tm/FD_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
  //H
  V_charm_1_distr_H = Corr.corr_t(V_charm_1_H.col(0)[i_ens], "../data/gm2/charm/tm/corr_1_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
  V_charm_2_distr_H = Corr.corr_t(V_charm_2_H.col(0)[i_ens], "../data/gm2/charm/tm/corr_2_"+V_charm_2_H.Tag[i_ens]+"_H.dat");
  V_charm_3_distr_H = Corr.corr_t(V_charm_3_H.col(0)[i_ens], "../data/gm2/charm/tm/corr_3_"+V_charm_3_H.Tag[i_ens]+"_H.dat");
  MD_distr_H = Corr.effective_mass_t(pt2_D_H.col(0)[i_ens], "../data/gm2/charm/tm/MD_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
  FD_distr_H = (L_info.ml + L_info.mc_H)*Corr.decay_constant_t(pt2_D_H.col(0)[i_ens], "../data/gm2/charm/tm/FD_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
  
  //OS
  Mpi_OS_distr = Corr.effective_mass_t(pt2_pion_OS_charm.col(0)[i_ens], "../data/gm2/charm/OS/Mpi_"+V_charm_1_L.Tag[i_ens]+".dat");
  //L
  V_charm_OS_1_distr_L = Corr.corr_t(V_charm_OS_1_L.col(0)[i_ens], "../data/gm2/charm/OS/corr_1_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
  V_charm_OS_2_distr_L = Corr.corr_t(V_charm_OS_2_L.col(0)[i_ens], "../data/gm2/charm/OS/corr_2_"+V_charm_2_L.Tag[i_ens]+"_L.dat");
  V_charm_OS_3_distr_L = Corr.corr_t(V_charm_OS_3_L.col(0)[i_ens], "../data/gm2/charm/OS/corr_3_"+V_charm_3_L.Tag[i_ens]+"_L.dat");
  //M
  V_charm_OS_1_distr_M = Corr.corr_t(V_charm_OS_1_M.col(0)[i_ens], "../data/gm2/charm/OS/corr_1_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
  V_charm_OS_2_distr_M = Corr.corr_t(V_charm_OS_2_M.col(0)[i_ens], "../data/gm2/charm/OS/corr_2_"+V_charm_2_M.Tag[i_ens]+"_M.dat");
  V_charm_OS_3_distr_M = Corr.corr_t(V_charm_OS_3_M.col(0)[i_ens], "../data/gm2/charm/OS/corr_3_"+V_charm_3_M.Tag[i_ens]+"_M.dat");
  //H
  V_charm_OS_1_distr_H = Corr.corr_t(V_charm_OS_1_H.col(0)[i_ens], "../data/gm2/charm/OS/corr_1_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
  V_charm_OS_2_distr_H = Corr.corr_t(V_charm_OS_2_H.col(0)[i_ens], "../data/gm2/charm/OS/corr_2_"+V_charm_2_H.Tag[i_ens]+"_H.dat");
  V_charm_OS_3_distr_H = Corr.corr_t(V_charm_OS_3_H.col(0)[i_ens], "../data/gm2/charm/OS/corr_3_"+V_charm_3_H.Tag[i_ens]+"_H.dat");


  //to compute Zv and Za
  ratio_P5P5_overlap_OS_tm_L = A0P5_OS_distr_L/A0P5_distr_L;
  ratio_P5P5_overlap_OS_tm_M = A0P5_OS_distr_M/A0P5_distr_M;
  ratio_P5P5_overlap_OS_tm_H = A0P5_OS_distr_H/A0P5_distr_H;

  
  //sum over the Lorenz indices of the e.m. current
  //tm
  //L
  V_charm_distr_L= (pow(qc,2)/3.0)*(V_charm_1_distr_L+ V_charm_2_distr_L + V_charm_3_distr_L);
  //M
  V_charm_distr_M= (pow(qc,2)/3.0)*(V_charm_1_distr_M+ V_charm_2_distr_M + V_charm_3_distr_M);
  //H
  V_charm_distr_H= (pow(qc,2)/3.0)*(V_charm_1_distr_H+ V_charm_2_distr_H + V_charm_3_distr_H);
  
  //OS
  //L
  V_charm_OS_distr_L= (pow(qc,2)/3.0)*(V_charm_OS_1_distr_L+ V_charm_OS_2_distr_L + V_charm_OS_3_distr_L);
  //M
  V_charm_OS_distr_M= (pow(qc,2)/3.0)*(V_charm_OS_1_distr_M+ V_charm_OS_2_distr_M + V_charm_OS_3_distr_M);
  //H
  V_charm_OS_distr_H= (pow(qc,2)/3.0)*(V_charm_OS_1_distr_H+ V_charm_OS_2_distr_H + V_charm_OS_3_distr_H);


  
  
  //#######################################  COMPUTATION OF ZV AND ZA (Hadronic method) ##################################
  //define lambda functions to be used
  auto sqr= [=](double a, double b) {return sqrt(a);};
  auto SINH= [](double m) -> double  {return sinh(m);};

  
  //take ratio between OS and tm pion amplitude to compute Zp/Zs RC.
  ratio_P5P5_overlap_OS_tm_L= overlap_P5P5_OS_distr_L/overlap_P5P5_distr_L;
  Zp_ov_Zs_distr_L = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm_L);
  ratio_P5P5_overlap_OS_tm_M= overlap_P5P5_OS_distr_M/overlap_P5P5_distr_M;
  Zp_ov_Zs_distr_M = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm_M);
  ratio_P5P5_overlap_OS_tm_H= overlap_P5P5_OS_distr_H/overlap_P5P5_distr_H;
  Zp_ov_Zs_distr_H = distr_t_list::f_of_distr_list(sqr, ratio_P5P5_overlap_OS_tm_H);



  //antysymmetrize w.r.t. t -> T-t for A0P5 correlators
  Corr.Reflection_sign = -1;
  A0P5_distr_L= Corr.corr_t(corr_A0P5_charm_L.col(0)[i_ens], "");
  A0P5_OS_distr_L = Corr.corr_t(corr_A0P5_OS_charm_L.col(0)[i_ens], "");
  A0P5_distr_M= Corr.corr_t(corr_A0P5_charm_M.col(0)[i_ens], "");
  A0P5_OS_distr_M = Corr.corr_t(corr_A0P5_OS_charm_M.col(0)[i_ens], "");
  A0P5_distr_H= Corr.corr_t(corr_A0P5_charm_H.col(0)[i_ens], "");
  A0P5_OS_distr_H = Corr.corr_t(corr_A0P5_OS_charm_H.col(0)[i_ens], "");

  
  //restore symmetrization
  Corr.Reflection_sign = 1;

  //compute RV (estimator for Zv)


  RV_L= 2.0*L_info.mc_L*cbar_c_distr_L/distr_t_list::derivative(A0P5_distr_L, 0); //central derivative
  RV_M= 2.0*L_info.mc_M*cbar_c_distr_M/distr_t_list::derivative(A0P5_distr_M, 0); //central derivative
  RV_H= 2.0*L_info.mc_H*cbar_c_distr_H/distr_t_list::derivative(A0P5_distr_H, 0); //central derivative
  
  //tm and OS P5P5
  cbar_c_mass_L = Corr.Fit_distr(M_etaC_distr_L);
  cbar_c_OS_mass_L= Corr.Fit_distr(M_etaC_distr_OS_L);
  cbar_c_mass_M = Corr.Fit_distr(M_etaC_distr_M);
  cbar_c_OS_mass_M= Corr.Fit_distr(M_etaC_distr_OS_M);
  cbar_c_mass_H = Corr.Fit_distr(M_etaC_distr_H);
  cbar_c_OS_mass_H= Corr.Fit_distr(M_etaC_distr_OS_H);
  //fit obs to compute Zv and Za (hadronic method)
  Zp_ov_Zs_L = Corr.Fit_distr(Zp_ov_Zs_distr_L);
  Zp_ov_Zs_M = Corr.Fit_distr(Zp_ov_Zs_distr_M);
  Zp_ov_Zs_H = Corr.Fit_distr(Zp_ov_Zs_distr_H);
  //set plateaux for RV
  int Tmin_RV, Tmax_RV;
  if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
    if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_RV=36; Tmax_RV=73;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
    if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_RV=28; Tmax_RV=57;}
    else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_RV=34;Tmax_RV=46;}
    else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_RV=35; Tmax_RV= 61;}
    else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.96") {Tmin_RV=35; Tmax_RV= 61;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
    if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_RV=28; Tmax_RV=46;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_RV=17; Tmax_RV=23;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_RV=16; Tmax_RV=23;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_RV=21;Tmax_RV=31;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "D") {
    if(V_charm_1_L.Tag[i_ens] == "cD211a.054.96") {Tmin_RV=32; Tmax_RV=70;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else crash("Ensemble tag not valid");

  Corr.Tmin = Tmin_RV;
  Corr.Tmax = Tmax_RV;
  

  Zv_hadr_L= Corr.Fit_distr(RV_L);
  Zv_hadr_M= Corr.Fit_distr(RV_M);
  Zv_hadr_H= Corr.Fit_distr(RV_H);
  
  RA_L = 2.0*L_info.mc_L*(cbar_c_OS_distr_L/distr_t_list::derivative(A0P5_OS_distr_L, 0))*(cbar_c_OS_mass_L/cbar_c_mass_L)*(distr_t::f_of_distr(SINH, cbar_c_OS_mass_L)/distr_t::f_of_distr(SINH, cbar_c_mass_L))*(1.0/Zp_ov_Zs_L);
  RA_M = 2.0*L_info.mc_M*(cbar_c_OS_distr_M/distr_t_list::derivative(A0P5_OS_distr_M, 0))*(cbar_c_OS_mass_M/cbar_c_mass_M)*(distr_t::f_of_distr(SINH, cbar_c_OS_mass_M)/distr_t::f_of_distr(SINH, cbar_c_mass_M))*(1.0/Zp_ov_Zs_M);
  RA_H = 2.0*L_info.mc_H*(cbar_c_OS_distr_H/distr_t_list::derivative(A0P5_OS_distr_H, 0))*(cbar_c_OS_mass_H/cbar_c_mass_H)*(distr_t::f_of_distr(SINH, cbar_c_OS_mass_H)/distr_t::f_of_distr(SINH, cbar_c_mass_H))*(1.0/Zp_ov_Zs_H);
  //set plateaux for RA
  int Tmin_RA=0;
  int Tmax_RA=0;
  //set time intervals for RA
  Corr.Tmin= Tmin_etaC;
  Corr.Tmax= Tmax_etaC;
  /*
  if(V_charm_1_L.Tag[i_ens].substr(1,1) == "C") {
    if(V_charm_1_L.Tag[i_ens]=="cC211a.06.80") { Tmin_RA=25; Tmax_RA=50;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "B") {
    if(V_charm_1_L.Tag[i_ens]== "cB211a.14.64") {Tmin_RA=18; Tmax_RA=26;}
    else if(V_charm_1_L.Tag[i_ens] == "cB211a.25.48") {Tmin_RA=10;Tmax_RA=21;}
    else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.64") {Tmin_RA=20; Tmax_RA= 50;}
    else if(V_charm_1_L.Tag[i_ens] == "cB211b.072.96") {Tmin_RA=20; Tmax_RA= 50;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "A") {
    if(V_charm_1_L.Tag[i_ens] == "cA211a.12.48") {Tmin_RA=15; Tmax_RA=25;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211a.40.24") {Tmin_RA=12; Tmax_RA=20;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211a.53.24") {Tmin_RA=11; Tmax_RA=20;}
    else if(V_charm_1_L.Tag[i_ens] == "cA211ab.30.32") {Tmin_RA=14;Tmax_RA=22;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else if(V_charm_1_L.Tag[i_ens].substr(1,1) == "D") {
    if(V_charm_1_L.Tag[i_ens] == "cD211a.054.96") {Tmin_RA=32; Tmax_RA=70;}
    else crash("Cannot find ensemble tag: "+V_charm_1_L.Tag[i_ens]);
  }
  else crash("Ensemble tag not valid");
  

  Corr.Tmin=Tmin_RA;
  Corr.Tmax=Tmax_RA;
  */
  
  Za_hadr_L= Corr.Fit_distr(RA_L);
  Za_hadr_M= Corr.Fit_distr(RA_M);
  Za_hadr_H= Corr.Fit_distr(RA_H);

  //print Rv and RA
  //print RV
  Print_To_File({}, {RV_L.ave(), RV_L.err(), RV_M.ave(), RV_M.err(), RV_H.ave(), RV_H.err()}, "../data/gm2/charm/RV_cc"+V_charm_1_L.Tag[i_ens]+".dat.t", "", "");
  //print RA
  Print_To_File({}, {RA_L.ave(), RA_L.err(), RA_M.ave(), RA_M.err(), RA_H.ave(), RA_H.err()}, "../data/gm2/charm/RA_cc"+V_charm_1_L.Tag[i_ens]+".dat.t", "", "");
  //print Zp_ov_Zs
  Print_To_File({}, {Zp_ov_Zs_distr_L.ave(), Zp_ov_Zs_distr_L.err(), Zp_ov_Zs_distr_M.ave(), Zp_ov_Zs_distr_M.err(), Zp_ov_Zs_distr_H.ave(), Zp_ov_Zs_distr_H.err()}, "../data/gm2/charm/Zp_ov_Zs_cc"+V_charm_1_L.Tag[i_ens]+".dat.t", "", "");

  //push_back

  Zv_fit_charm_L.distr_list.push_back(Zv_hadr_L);
  Za_fit_charm_L.distr_list.push_back(Za_hadr_L);
  Zv_fit_charm_M.distr_list.push_back(Zv_hadr_M);
  Za_fit_charm_M.distr_list.push_back(Za_hadr_M);
  Zv_fit_charm_H.distr_list.push_back(Zv_hadr_H);
  Za_fit_charm_H.distr_list.push_back(Za_hadr_H);

  if(Use_Za_Zv_from_charm_run) { Zv_L= Zv_hadr_L; Za_L=Za_hadr_L; Zv_M = Zv_hadr_M; Za_M = Za_hadr_M; Zv_H = Zv_hadr_H; Za_H = Za_hadr_H;}


  //################################################ END OF COMPUTATION OF ZV AND ZA (Hadronic method) #####################################à




  

  //free corr LO artifacts
  Vfloat free_corr_log_art(Corr.Nt, 0.0);
  for(int t=0;t<Corr.Nt;t++) {  if( t*a_distr.ave()/fm_to_inv_Gev < 1.0 && t != 0) {free_corr_log_art[t] = -1.0*(qc*qc)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));}}
  
  //L
  distr_t_list V_charm_distr_L_pert_sub = V_charm_distr_L + free_corr_log_art;
  distr_t_list V_charm_OS_distr_L_pert_sub = V_charm_OS_distr_L + free_corr_log_art;
  //M
  distr_t_list V_charm_distr_M_pert_sub = V_charm_distr_M + free_corr_log_art;
  distr_t_list V_charm_OS_distr_M_pert_sub = V_charm_OS_distr_M + free_corr_log_art;
  //H
  distr_t_list V_charm_distr_H_pert_sub = V_charm_distr_H + free_corr_log_art;
  distr_t_list V_charm_OS_distr_H_pert_sub = V_charm_OS_distr_H + free_corr_log_art;

 

  // print summed correlators to file
  //tm
  //L
  Print_To_File({}, {V_charm_distr_L.ave(), V_charm_distr_L.err(), (Za_L*Za_L*V_charm_distr_L).ave(), (Za_L*Za_L*V_charm_distr_L).err()}, "../data/gm2/charm/tm/corr_sum_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare  renormalized");
  //M
  Print_To_File({}, {V_charm_distr_M.ave(), V_charm_distr_M.err(), (Za_M*Za_M*V_charm_distr_M).ave(), (Za_M*Za_M*V_charm_distr_M).err()}, "../data/gm2/charm/tm/corr_sum_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#t   bare  renormalized");
  //H
  Print_To_File({}, {V_charm_distr_H.ave(), V_charm_distr_H.err(), (Za_H*Za_H*V_charm_distr_H).ave(), (Za_H*Za_H*V_charm_distr_H).err()}, "../data/gm2/charm/tm/corr_sum_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#t   bare  renormalized");

  
  //OS
  //L
  Print_To_File({}, {V_charm_OS_distr_L.ave(), V_charm_OS_distr_L.err(), (Zv_L*Zv_L*V_charm_OS_distr_L).ave(), (Zv_H*Zv_H*V_charm_OS_distr_L).err()}, "../data/gm2/charm/OS/corr_sum_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#t   bare  renormalized");
  //M
  Print_To_File({}, {V_charm_OS_distr_M.ave(), V_charm_OS_distr_M.err(), (Zv_M*Zv_M*V_charm_OS_distr_M).ave(), (Zv_H*Zv_H*V_charm_OS_distr_M).err()}, "../data/gm2/charm/OS/corr_sum_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#t   bare  renormalized");
  //H
  Print_To_File({}, {V_charm_OS_distr_H.ave(), V_charm_OS_distr_H.err(), (Zv_H*Zv_H*V_charm_OS_distr_H).ave(), (Zv_H*Zv_H*V_charm_OS_distr_H).err()}, "../data/gm2/charm/OS/corr_sum_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#t   bare  renormalized");

 //extract effective masses, overlap from V and fit

  


  //tm

  //L
  Corr.Tmin=Tmin_P5P5;
  Corr.Tmax=Tmax_P5P5;
  Mpi = Corr.Fit_distr(Mpi_distr);
  MD_L= Corr.Fit_distr(MD_distr_L);
  FD_L= Corr.Fit_distr(FD_distr_L);
  Corr.Tmin=Tmin_etaC;
  Corr.Tmax=Tmax_etaC;
  M_etaC_L = Corr.Fit_distr(M_etaC_distr_L);
  Corr.Tmin= Tmin_VV;
  Corr.Tmax= Tmax_VV;
  MV_charm_distr_L= Corr.effective_mass_t(V_charm_distr_L, "../data/gm2/charm/tm/MV_mass_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
  ZV_charm_distr_L= Corr.residue_t(V_charm_distr_L, "../data/gm2/charm/tm/ZV_overlap_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
  MV_charm_L = Corr.Fit_distr(MV_charm_distr_L);
  ZV_charm_L = Corr.Fit_distr(ZV_charm_distr_L);
 
  //M
  Corr.Tmin=Tmin_P5P5;
  Corr.Tmax=Tmax_P5P5;
  MD_M= Corr.Fit_distr(MD_distr_M);
  FD_M= Corr.Fit_distr(FD_distr_M);
  Corr.Tmin=Tmin_etaC;
  Corr.Tmax=Tmax_etaC;
  M_etaC_M = Corr.Fit_distr(M_etaC_distr_M);
  Corr.Tmin=Tmin_VV;
  Corr.Tmax=Tmax_VV;
  MV_charm_distr_M= Corr.effective_mass_t(V_charm_distr_M, "../data/gm2/charm/tm/MV_mass_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
  ZV_charm_distr_M= Corr.residue_t(V_charm_distr_M, "../data/gm2/charm/tm/ZV_overlap_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
  MV_charm_M = Corr.Fit_distr(MV_charm_distr_M);
  ZV_charm_M = Corr.Fit_distr(ZV_charm_distr_M);
 
  //H
  Corr.Tmin=Tmin_P5P5;
  Corr.Tmax=Tmax_P5P5;
  MD_H= Corr.Fit_distr(MD_distr_H);
  FD_H= Corr.Fit_distr(FD_distr_H);
  Corr.Tmin=Tmin_etaC;
  Corr.Tmax=Tmax_etaC;
  M_etaC_H = Corr.Fit_distr(M_etaC_distr_H);
  Corr.Tmin=Tmin_VV;
  Corr.Tmax=Tmax_VV;
  MV_charm_distr_H= Corr.effective_mass_t(V_charm_distr_H, "../data/gm2/charm/tm/MV_mass_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
  ZV_charm_distr_H= Corr.residue_t(V_charm_distr_H, "../data/gm2/charm/tm/ZV_overlap_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
  MV_charm_H = Corr.Fit_distr(MV_charm_distr_H);
  ZV_charm_H = Corr.Fit_distr(ZV_charm_distr_H);



  //push back to Csi_D vector
  vector<distr_t> Csi_D({MD_L*MD_L/(FD_L*FD_L), MD_M*MD_M/(FD_M*FD_M), MD_H*MD_H/(FD_H*FD_H)});
  //push back to Md^2 vector (physical units)
  vector<distr_t> M2_D({ MD_L*MD_L/(a_distr*a_distr), MD_M*MD_M/(a_distr*a_distr), MD_H*MD_H/(a_distr*a_distr)});
  //push back to M_etac^2 vector (physical units)
  vector<distr_t> M2_etaC({ M_etaC_L*M_etaC_L/(a_distr*a_distr), M_etaC_M*M_etaC_M/(a_distr*a_distr), M_etaC_H*M_etaC_H/(a_distr*a_distr)});
  vector<distr_t> Jpsi2( { (MV_charm_L*MV_charm_L)/(a_distr*a_distr),  (MV_charm_M*MV_charm_M)/(a_distr*a_distr),  (MV_charm_H*MV_charm_H)/(a_distr*a_distr)});


  vector<distr_t> X_2_fit = M2_etaC; // M2_D; //Jpsi2; //M2_D; //M2_etaC;
  double X_2_phys= m_etac*m_etac; //m_etac*m_etac; //m_etac*m_etac; //m_d*m_d; //m_d*m_d; //m_etac*m_etac;

  
  //push_back MV_charm, ZV_charm, Mpi and MD
 
  Mpi_fit_charm.distr_list.push_back(Mpi);
  //L
  MV_fit_charm_L.distr_list.push_back(MV_charm_L);
  ZV_fit_charm_L.distr_list.push_back(Za_L*Za_L*ZV_charm_L);
  MD_fit_L.distr_list.push_back(MD_L);
  //M
  MV_fit_charm_M.distr_list.push_back(MV_charm_M);
  ZV_fit_charm_M.distr_list.push_back(Za_M*Za_M*ZV_charm_M);
  MD_fit_M.distr_list.push_back(MD_M);
  //H
  MV_fit_charm_H.distr_list.push_back(MV_charm_H);
  ZV_fit_charm_H.distr_list.push_back(Za_H*Za_H*ZV_charm_H);
  MD_fit_H.distr_list.push_back(MD_H);

  //OS
  Corr.Tmin=Tmin_P5P5;
  Corr.Tmax=Tmax_P5P5;
  Mpi_OS = Corr.Fit_distr(Mpi_OS_distr);
  Corr.Tmin=Tmin_VV;
  Corr.Tmax=Tmax_VV;
  //L
  MV_charm_OS_distr_L= Corr.effective_mass_t(V_charm_OS_distr_L, "../data/gm2/charm/OS/MV_mass_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
  ZV_charm_OS_distr_L= Corr.residue_t(V_charm_OS_distr_L, "../data/gm2/charm/OS/ZV_overlap_"+V_charm_1_L.Tag[i_ens]+"_L.dat");
  MV_charm_OS_L = Corr.Fit_distr(MV_charm_OS_distr_L);
  ZV_charm_OS_L = Corr.Fit_distr(ZV_charm_OS_distr_L);
  //M
  MV_charm_OS_distr_M= Corr.effective_mass_t(V_charm_OS_distr_M, "../data/gm2/charm/OS/MV_mass_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
  ZV_charm_OS_distr_M= Corr.residue_t(V_charm_OS_distr_M, "../data/gm2/charm/OS/ZV_overlap_"+V_charm_1_M.Tag[i_ens]+"_M.dat");
  MV_charm_OS_M = Corr.Fit_distr(MV_charm_OS_distr_M);
  ZV_charm_OS_M = Corr.Fit_distr(ZV_charm_OS_distr_M);
  //H
  MV_charm_OS_distr_H= Corr.effective_mass_t(V_charm_OS_distr_H, "../data/gm2/charm/OS/MV_mass_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
  ZV_charm_OS_distr_H= Corr.residue_t(V_charm_OS_distr_H, "../data/gm2/charm/OS/ZV_overlap_"+V_charm_1_H.Tag[i_ens]+"_H.dat");
  MV_charm_OS_H = Corr.Fit_distr(MV_charm_OS_distr_H);
  ZV_charm_OS_H = Corr.Fit_distr(ZV_charm_OS_distr_H);

  

  //push_back MV_charm, ZV_charm Mpi_OS
  Mpi_OS_fit_charm.distr_list.push_back(Mpi_OS);
  //L
  MV_fit_charm_OS_L.distr_list.push_back(MV_charm_OS_L);
  ZV_fit_charm_OS_L.distr_list.push_back(Zv_L*Zv_L*ZV_charm_OS_L);
  //M
  MV_fit_charm_OS_M.distr_list.push_back(MV_charm_OS_M);
  ZV_fit_charm_OS_M.distr_list.push_back(Zv_M*Zv_M*ZV_charm_OS_M);
  //H
  MV_fit_charm_OS_H.distr_list.push_back(MV_charm_OS_H);
  ZV_fit_charm_OS_H.distr_list.push_back(Zv_H*Zv_H*ZV_charm_OS_H);
 


  int Tdata_min= 8;
  int Tdata_max = Corr.Nt/2.0 -1;
  int Tdata_fit = (Tmin_VV+Tmax_VV)/2;
  
  //compute kernel distribution
  //NON-ELM
  distr_t_list Kernel_distr = distr_t_list::f_of_distr(K, a_distr, Upper_Limit_Time_Integral_charm+1);

  
  //tm
  //L
  distr_t_list Kernel_distr_L = distr_t_list::f_of_distr(K,MV_charm_L/m_Jpsi, Upper_Limit_Time_Integral_charm+1);
  //M
  distr_t_list Kernel_distr_M = distr_t_list::f_of_distr(K,MV_charm_M/m_Jpsi, Upper_Limit_Time_Integral_charm+1);
  //H
  distr_t_list Kernel_distr_H = distr_t_list::f_of_distr(K,MV_charm_H/m_Jpsi, Upper_Limit_Time_Integral_charm+1);

  
  //OS
  //L
  distr_t_list Kernel_OS_distr_L = distr_t_list::f_of_distr(K,MV_charm_OS_L/m_Jpsi, Upper_Limit_Time_Integral_charm+1);
  //M
  distr_t_list Kernel_OS_distr_M = distr_t_list::f_of_distr(K,MV_charm_OS_M/m_Jpsi, Upper_Limit_Time_Integral_charm+1);
  //H
  distr_t_list Kernel_OS_distr_H = distr_t_list::f_of_distr(K,MV_charm_OS_H/m_Jpsi, Upper_Limit_Time_Integral_charm+1);


  
  //compute exp(-Mv*t) distribution
  //tm
  //L
  distr_t_list exp_MVc_L = distr_t_list::f_of_distr(exp_MV, MV_charm_L, Upper_Limit_Time_Integral_charm+1);
  //M
  distr_t_list exp_MVc_M = distr_t_list::f_of_distr(exp_MV, MV_charm_M, Upper_Limit_Time_Integral_charm+1);
  //H
  distr_t_list exp_MVc_H = distr_t_list::f_of_distr(exp_MV, MV_charm_H, Upper_Limit_Time_Integral_charm+1);
 
  //OS
  //L
  distr_t_list exp_OS_MVc_L = distr_t_list::f_of_distr(exp_MV, MV_charm_OS_L, Upper_Limit_Time_Integral_charm+1);
  //M
  distr_t_list exp_OS_MVc_M = distr_t_list::f_of_distr(exp_MV, MV_charm_OS_M, Upper_Limit_Time_Integral_charm+1);
  //H
  distr_t_list exp_OS_MVc_H = distr_t_list::f_of_distr(exp_MV, MV_charm_OS_H, Upper_Limit_Time_Integral_charm+1);

  
  //Print single-exponential prediction to file
  //tm
  //L
  Print_To_File({}, {(exp_MVc_L*(ZV_charm_L/(2.0*MV_charm_L))).ave(), (exp_MVc_L*(ZV_charm_L/(2.0*MV_charm_L))).err(), (Za_L*Za_L*exp_MVc_L*(ZV_charm_L/(2.0*MV_charm_L))).ave(), (Za_L*Za_L*exp_MVc_L*(ZV_charm_L/(2.0*MV_charm_L))).err() }, "../data/gm2/charm/tm/corr_gsd_sum_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#t  bare  renormalized");
  //M
  Print_To_File({}, {(exp_MVc_M*(ZV_charm_M/(2.0*MV_charm_M))).ave(), (exp_MVc_M*(ZV_charm_M/(2.0*MV_charm_M))).err(), (Za_M*Za_M*exp_MVc_M*(ZV_charm_M/(2.0*MV_charm_M))).ave(), (Za_M*Za_M*exp_MVc_M*(ZV_charm_M/(2.0*MV_charm_M))).err() }, "../data/gm2/charm/tm/corr_gsd_sum_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#t  bare  renormalized");
  //H
  Print_To_File({}, {(exp_MVc_H*(ZV_charm_H/(2.0*MV_charm_H))).ave(), (exp_MVc_H*(ZV_charm_H/(2.0*MV_charm_H))).err(), (Za_H*Za_H*exp_MVc_H*(ZV_charm_H/(2.0*MV_charm_H))).ave(), (Za_H*Za_H*exp_MVc_H*(ZV_charm_H/(2.0*MV_charm_H))).err() }, "../data/gm2/charm/tm/corr_gsd_sum_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#t  bare  renormalized");

  
  //OS
  //L
  Print_To_File({}, {(exp_OS_MVc_L*(ZV_charm_OS_L/(2.0*MV_charm_OS_L))).ave(), (exp_OS_MVc_L*(ZV_charm_OS_L/(2.0*MV_charm_OS_L))).err(), (Zv_L*Zv_L*exp_OS_MVc_L*(ZV_charm_OS_L/(2.0*MV_charm_OS_L))).ave(), (Zv_L*Zv_L*exp_OS_MVc_L*(ZV_charm_OS_L/(2.0*MV_charm_OS_L))).err() }, "../data/gm2/charm/OS/corr_gsd_sum_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#t  bare  renormalized");
  //M
  Print_To_File({}, {(exp_OS_MVc_M*(ZV_charm_OS_M/(2.0*MV_charm_OS_M))).ave(), (exp_OS_MVc_M*(ZV_charm_OS_M/(2.0*MV_charm_OS_M))).err(), (Zv_M*Zv_M*exp_OS_MVc_M*(ZV_charm_OS_M/(2.0*MV_charm_OS_M))).ave(), (Zv_M*Zv_M*exp_OS_MVc_M*(ZV_charm_OS_M/(2.0*MV_charm_OS_M))).err() }, "../data/gm2/charm/OS/corr_gsd_sum_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#t  bare  renormalized");
  //H
  Print_To_File({}, {(exp_OS_MVc_H*(ZV_charm_OS_H/(2.0*MV_charm_OS_H))).ave(), (exp_OS_MVc_H*(ZV_charm_OS_H/(2.0*MV_charm_OS_H))).err(), (Zv_H*Zv_H*exp_OS_MVc_H*(ZV_charm_OS_H/(2.0*MV_charm_OS_H))).ave(), (Zv_H*Zv_H*exp_OS_MVc_H*(ZV_charm_OS_H/(2.0*MV_charm_OS_H))).err() }, "../data/gm2/charm/OS/corr_gsd_sum_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#t  bare  renormalized");
  

  distr_t_list agm2_distr_Tdata_L(UseJack), agm2_OS_distr_Tdata_L(UseJack);
  distr_t_list agm2_distr_Tdata_M(UseJack), agm2_OS_distr_Tdata_M(UseJack);
  distr_t_list agm2_distr_Tdata_H(UseJack), agm2_OS_distr_Tdata_H(UseJack);
  distr_t_list agm2_distr_Tdata_No_ELM_L(UseJack), agm2_OS_distr_Tdata_No_ELM_L(UseJack);
  distr_t_list agm2_distr_Tdata_No_ELM_M(UseJack), agm2_OS_distr_Tdata_No_ELM_M(UseJack);
  distr_t_list agm2_distr_Tdata_No_ELM_H(UseJack), agm2_OS_distr_Tdata_No_ELM_H(UseJack);
  Vfloat Tdata_vec;
  bool Find_Tdata_fit=false;
  
  for(int Tdata=Tdata_min;Tdata<Tdata_max;Tdata++) {
    //compute 4\pia^2 using lattice data up to Tdata (included)

    //L
    distr_t agm2_L(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_OS_L(UseJack, UseJack?Njacks:Nboots);
    //M
    distr_t agm2_M(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_OS_M(UseJack, UseJack?Njacks:Nboots);
    //H
    distr_t agm2_H(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_OS_H(UseJack, UseJack?Njacks:Nboots);

    //NON-ELM
    //L
    distr_t agm2_No_ELM_L(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_OS_No_ELM_L(UseJack, UseJack?Njacks:Nboots);
    //M
    distr_t agm2_No_ELM_M(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_OS_No_ELM_M(UseJack, UseJack?Njacks:Nboots);
    //H
    distr_t agm2_No_ELM_H(UseJack, UseJack?Njacks:Nboots);
    distr_t agm2_OS_No_ELM_H(UseJack, UseJack?Njacks:Nboots);

    
    for(int t=1;t<=Upper_Limit_Time_Integral_charm;t++) {
      if(t<=Tdata) {
	//L
	agm2_L = agm2_L + 4.0*pow(alpha,2)*V_charm_distr_L_pert_sub.distr_list[t]*Kernel_distr_L.distr_list[t];
	agm2_OS_L = agm2_OS_L + 4.0*pow(alpha,2)*V_charm_OS_distr_L_pert_sub.distr_list[t]*Kernel_OS_distr_L.distr_list[t];
	//M
	agm2_M = agm2_M + 4.0*pow(alpha,2)*V_charm_distr_M_pert_sub.distr_list[t]*Kernel_distr_M.distr_list[t];
	agm2_OS_M = agm2_OS_M + 4.0*pow(alpha,2)*V_charm_OS_distr_M_pert_sub.distr_list[t]*Kernel_OS_distr_M.distr_list[t];
	//H
	agm2_H = agm2_H + 4.0*pow(alpha,2)*V_charm_distr_H_pert_sub.distr_list[t]*Kernel_distr_H.distr_list[t];
	agm2_OS_H = agm2_OS_H + 4.0*pow(alpha,2)*V_charm_OS_distr_H_pert_sub.distr_list[t]*Kernel_OS_distr_H.distr_list[t];

	//NON_ELM
	//L
	agm2_No_ELM_L = agm2_No_ELM_L + 4.0*pow(alpha,2)*V_charm_distr_L_pert_sub.distr_list[t]*Kernel_distr.distr_list[t];
	agm2_OS_No_ELM_L = agm2_OS_No_ELM_L + 4.0*pow(alpha,2)*V_charm_OS_distr_L_pert_sub.distr_list[t]*Kernel_distr.distr_list[t];
	//M
	agm2_No_ELM_M = agm2_No_ELM_M + 4.0*pow(alpha,2)*V_charm_distr_M_pert_sub.distr_list[t]*Kernel_distr.distr_list[t];
	agm2_OS_No_ELM_M = agm2_OS_No_ELM_M + 4.0*pow(alpha,2)*V_charm_OS_distr_M_pert_sub.distr_list[t]*Kernel_distr.distr_list[t];
	//H
	agm2_No_ELM_H = agm2_No_ELM_H + 4.0*pow(alpha,2)*V_charm_distr_H_pert_sub.distr_list[t]*Kernel_distr.distr_list[t];
	agm2_OS_No_ELM_H = agm2_OS_No_ELM_H + 4.0*pow(alpha,2)*V_charm_OS_distr_H_pert_sub.distr_list[t]*Kernel_distr.distr_list[t];
	
      }
      else {
	//L
	agm2_L= agm2_L + 4.0*pow(alpha,2)*(ZV_charm_L/(2.0*MV_charm_L))*exp_MVc_L.distr_list[t]*Kernel_distr_L.distr_list[t];
	agm2_OS_L= agm2_OS_L + 4.0*pow(alpha,2)*(ZV_charm_OS_L/(2.0*MV_charm_OS_L))*exp_OS_MVc_L.distr_list[t]*Kernel_OS_distr_L.distr_list[t];
	//M
	agm2_M= agm2_M + 4.0*pow(alpha,2)*(ZV_charm_M/(2.0*MV_charm_M))*exp_MVc_M.distr_list[t]*Kernel_distr_M.distr_list[t];
	agm2_OS_M= agm2_OS_M + 4.0*pow(alpha,2)*(ZV_charm_OS_M/(2.0*MV_charm_OS_M))*exp_OS_MVc_M.distr_list[t]*Kernel_OS_distr_M.distr_list[t];
	//H
	agm2_H= agm2_H + 4.0*pow(alpha,2)*(ZV_charm_H/(2.0*MV_charm_H))*exp_MVc_H.distr_list[t]*Kernel_distr_H.distr_list[t];
	agm2_OS_H= agm2_OS_H + 4.0*pow(alpha,2)*(ZV_charm_OS_H/(2.0*MV_charm_OS_H))*exp_OS_MVc_H.distr_list[t]*Kernel_OS_distr_H.distr_list[t];

	//NON_ELM
	//L
	agm2_No_ELM_L= agm2_No_ELM_L + 4.0*pow(alpha,2)*(ZV_charm_L/(2.0*MV_charm_L))*exp_MVc_L.distr_list[t]*Kernel_distr.distr_list[t];
	agm2_OS_No_ELM_L= agm2_OS_No_ELM_L + 4.0*pow(alpha,2)*(ZV_charm_OS_L/(2.0*MV_charm_OS_L))*exp_OS_MVc_L.distr_list[t]*Kernel_distr.distr_list[t];
	//M
	agm2_No_ELM_M= agm2_No_ELM_M + 4.0*pow(alpha,2)*(ZV_charm_M/(2.0*MV_charm_M))*exp_MVc_M.distr_list[t]*Kernel_distr.distr_list[t];
	agm2_OS_No_ELM_M= agm2_OS_No_ELM_M + 4.0*pow(alpha,2)*(ZV_charm_OS_M/(2.0*MV_charm_OS_M))*exp_OS_MVc_M.distr_list[t]*Kernel_distr.distr_list[t];
	//H
	agm2_No_ELM_H= agm2_No_ELM_H + 4.0*pow(alpha,2)*(ZV_charm_H/(2.0*MV_charm_H))*exp_MVc_H.distr_list[t]*Kernel_distr.distr_list[t];
	agm2_OS_No_ELM_H= agm2_OS_No_ELM_H + 4.0*pow(alpha,2)*(ZV_charm_OS_H/(2.0*MV_charm_OS_H))*exp_OS_MVc_H.distr_list[t]*Kernel_distr.distr_list[t];
	
      }
    }
    
    Tdata_vec.push_back((double)Tdata);
    //L
    agm2_distr_Tdata_L.distr_list.push_back(Za_L*Za_L*agm2_L);
    agm2_OS_distr_Tdata_L.distr_list.push_back(Zv_L*Zv_L*agm2_OS_L);
    //M
    agm2_distr_Tdata_M.distr_list.push_back(Za_M*Za_M*agm2_M);
    agm2_OS_distr_Tdata_M.distr_list.push_back(Zv_M*Zv_M*agm2_OS_M);
    //H
    agm2_distr_Tdata_H.distr_list.push_back(Za_H*Za_H*agm2_H);
    agm2_OS_distr_Tdata_H.distr_list.push_back(Zv_H*Zv_H*agm2_OS_H);

    //NON_ELM
    //L
    agm2_distr_Tdata_No_ELM_L.distr_list.push_back(Za_L*Za_L*agm2_No_ELM_L);
    agm2_OS_distr_Tdata_No_ELM_L.distr_list.push_back(Zv_L*Zv_L*agm2_OS_No_ELM_L);
    //M
    agm2_distr_Tdata_No_ELM_M.distr_list.push_back(Za_M*Za_M*agm2_No_ELM_M);
    agm2_OS_distr_Tdata_No_ELM_M.distr_list.push_back(Zv_M*Zv_M*agm2_OS_No_ELM_M);
    //H
    agm2_distr_Tdata_No_ELM_H.distr_list.push_back(Za_H*Za_H*agm2_No_ELM_H);
    agm2_OS_distr_Tdata_No_ELM_H.distr_list.push_back(Zv_H*Zv_H*agm2_OS_No_ELM_H);

    
    if(Tdata==Tdata_fit) {
      //L
      agm2_charm_L.distr_list.push_back(Za_L*Za_L*agm2_L);
      agm2_charm_OS_L.distr_list.push_back(Zv_L*Zv_L*agm2_OS_L);
      //M
      agm2_charm_M.distr_list.push_back(Za_M*Za_M*agm2_M);
      agm2_charm_OS_M.distr_list.push_back(Zv_M*Zv_M*agm2_OS_M);
      //H
      agm2_charm_H.distr_list.push_back(Za_H*Za_H*agm2_H);
      agm2_charm_OS_H.distr_list.push_back(Zv_H*Zv_H*agm2_OS_H);
      
      //NON_ELM
      //L
      agm2_charm_No_ELM_L.distr_list.push_back(Za_L*Za_L*agm2_No_ELM_L);
      agm2_charm_OS_No_ELM_L.distr_list.push_back(Zv_L*Zv_L*agm2_OS_No_ELM_L);
      //M
      agm2_charm_No_ELM_M.distr_list.push_back(Za_M*Za_M*agm2_No_ELM_M);
      agm2_charm_OS_No_ELM_M.distr_list.push_back(Zv_M*Zv_M*agm2_OS_No_ELM_M);
      //H
      agm2_charm_No_ELM_H.distr_list.push_back(Za_H*Za_H*agm2_No_ELM_H);
      agm2_charm_OS_No_ELM_H.distr_list.push_back(Zv_H*Zv_H*agm2_OS_No_ELM_H);

   
      
      //extrapolate to the physical D-meson point
      vector<distr_t> agm2s_charm({Za_L*Za_L*agm2_L, Za_M*Za_M*agm2_M, Za_H*Za_H*agm2_H});
      vector<distr_t> agm2s_charm_OS({ Zv_L*Zv_L*agm2_OS_L, Zv_M*Zv_M*agm2_OS_M, Zv_H*Zv_H*agm2_OS_H});
      vector<distr_t> agm2s_charm_No_ELM({Za_L*Za_L*agm2_No_ELM_L, Za_M*Za_M*agm2_No_ELM_M, Za_H*Za_H*agm2_No_ELM_H});
      vector<distr_t> agm2s_charm_OS_No_ELM({ Zv_L*Zv_L*agm2_OS_No_ELM_L, Zv_M*Zv_M*agm2_OS_No_ELM_M, Zv_H*Zv_H*agm2_OS_No_ELM_H});
      agm2_charm_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_ELM_"+V_charm_1_L.Tag[i_ens], UseJack));
      agm2_charm_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_OS, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_OS_ELM_"+V_charm_1_L.Tag[i_ens], UseJack));
      agm2_charm_No_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_No_ELM, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_"+V_charm_1_L.Tag[i_ens], UseJack));
      agm2_charm_OS_No_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_OS_No_ELM, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_OS_"+V_charm_1_L.Tag[i_ens], UseJack));
      Find_Tdata_fit=true;
    }
  }

  if(!Find_Tdata_fit) crash("Cannot find the Tdata value: "+to_string(Tdata_fit));



  
   //#######################  INTERMEDIATE AND SHORT-DISTANCE WINDOW ###################################

    
  //############################   TWISTED MASS ######################################################

  //L
  distr_t agm2_W_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
  distr_t agm2_SD_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
  distr_t agm2_W_ELM_L(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
  distr_t agm2_SD_ELM_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default
  //M
  distr_t agm2_W_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
  distr_t agm2_SD_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
  distr_t agm2_W_ELM_M(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
  distr_t agm2_SD_ELM_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default
  //H
  distr_t agm2_W_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
  distr_t agm2_SD_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
  distr_t agm2_W_ELM_H(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
  distr_t agm2_SD_ELM_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default


  //#################################################################################################


  
  distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);
  //L
  distr_t_list Ker_ELM_tm_L = distr_t_list::f_of_distr(K, MV_charm_L/m_Jpsi, Corr.Nt/2);
  distr_t_list Ker_ELM_OS_L = distr_t_list::f_of_distr(K, MV_charm_OS_L/m_Jpsi, Corr.Nt/2);
  //M
  distr_t_list Ker_ELM_tm_M = distr_t_list::f_of_distr(K, MV_charm_M/m_Jpsi, Corr.Nt/2);
  distr_t_list Ker_ELM_OS_M = distr_t_list::f_of_distr(K, MV_charm_OS_M/m_Jpsi, Corr.Nt/2);
  //H
  distr_t_list Ker_ELM_tm_H = distr_t_list::f_of_distr(K, MV_charm_H/m_Jpsi, Corr.Nt/2);
  distr_t_list Ker_ELM_OS_H = distr_t_list::f_of_distr(K, MV_charm_OS_H/m_Jpsi, Corr.Nt/2);
    
  //define lambdas for the theta func
  auto th0 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t0)/Delta));};
  auto th1 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t1)/Delta));};

  
  for(int t=1; t< Corr.Nt/2; t++) {
    //L
    agm2_W_L = agm2_W_L + 4.0*pow(alpha,2)*Za_L*Za_L*V_charm_distr_L.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_L = agm2_SD_L + 4.0*pow(alpha,2)*Za_L*Za_L*(V_charm_distr_L_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_L = agm2_W_ELM_L + 4.0*pow(alpha,2)*Za_L*Za_L*V_charm_distr_L.distr_list[t]*Ker_ELM_tm_L.distr_list[t]*( distr_t::f_of_distr(th0, t*MV_charm_L/m_Jpsi) - distr_t::f_of_distr(th1, t*MV_charm_L/m_Jpsi));
    agm2_SD_ELM_L = agm2_SD_ELM_L + 4.0*pow(alpha,2)*(Za_L*Za_L*V_charm_distr_L_pert_sub.distr_list[t])*Ker_ELM_tm_L.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*MV_charm_L/m_Jpsi));
    //M
    agm2_W_M = agm2_W_M + 4.0*pow(alpha,2)*Za_M*Za_M*V_charm_distr_M.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_M = agm2_SD_M + 4.0*pow(alpha,2)*Za_M*Za_M*(V_charm_distr_M_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_M = agm2_W_ELM_M + 4.0*pow(alpha,2)*Za_M*Za_M*V_charm_distr_M.distr_list[t]*Ker_ELM_tm_M.distr_list[t]*( distr_t::f_of_distr(th0, t*MV_charm_M/m_Jpsi) - distr_t::f_of_distr(th1, t*MV_charm_M/m_Jpsi));
    agm2_SD_ELM_M = agm2_SD_ELM_M + 4.0*pow(alpha,2)*(Za_M*Za_M*V_charm_distr_M_pert_sub.distr_list[t])*Ker_ELM_tm_M.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*MV_charm_M/m_Jpsi));
    //H
    agm2_W_H = agm2_W_H + 4.0*pow(alpha,2)*Za_H*Za_H*V_charm_distr_H.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_H = agm2_SD_H + 4.0*pow(alpha,2)*Za_H*Za_H*(V_charm_distr_H_pert_sub.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_H = agm2_W_ELM_H + 4.0*pow(alpha,2)*Za_H*Za_H*V_charm_distr_H.distr_list[t]*Ker_ELM_tm_H.distr_list[t]*( distr_t::f_of_distr(th0, t*MV_charm_H/m_Jpsi) - distr_t::f_of_distr(th1, t*MV_charm_H/m_Jpsi));
    agm2_SD_ELM_H = agm2_SD_ELM_H + 4.0*pow(alpha,2)*(Za_H*Za_H*V_charm_distr_H_pert_sub.distr_list[t])*Ker_ELM_tm_H.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*MV_charm_H/m_Jpsi));
  }
  
  //push_back the result

  //L
  agm2_charm_W_L.distr_list.push_back(agm2_W_L);
  agm2_charm_SD_L.distr_list.push_back(agm2_SD_L);
  agm2_charm_W_ELM_L.distr_list.push_back(agm2_W_ELM_L);
  agm2_charm_SD_ELM_L.distr_list.push_back(agm2_SD_ELM_L);
  //M
  agm2_charm_W_M.distr_list.push_back(agm2_W_M);
  agm2_charm_SD_M.distr_list.push_back(agm2_SD_M);
  agm2_charm_W_ELM_M.distr_list.push_back(agm2_W_ELM_M);
  agm2_charm_SD_ELM_M.distr_list.push_back(agm2_SD_ELM_M);
  //H
  agm2_charm_W_H.distr_list.push_back(agm2_W_H);
  agm2_charm_SD_H.distr_list.push_back(agm2_SD_H);
  agm2_charm_W_ELM_H.distr_list.push_back(agm2_W_ELM_H);
  agm2_charm_SD_ELM_H.distr_list.push_back(agm2_SD_ELM_H);

  
  //extrapolate the result to the physical D-meson point

  vector<distr_t> agm2s_charm_W({agm2_W_L, agm2_W_M, agm2_W_H});
  vector<distr_t> agm2s_charm_SD({agm2_SD_L, agm2_SD_M, agm2_SD_H});
  vector<distr_t> agm2s_charm_W_ELM({agm2_W_ELM_L, agm2_W_ELM_M, agm2_W_ELM_H});
  vector<distr_t> agm2s_charm_SD_ELM({agm2_SD_ELM_L, agm2_SD_ELM_M, agm2_SD_ELM_H});

  
  agm2_charm_W_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_W, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_W_"+V_charm_1_L.Tag[i_ens], UseJack));
  agm2_charm_SD_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_SD, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_SD_"+V_charm_1_L.Tag[i_ens], UseJack));
  agm2_charm_W_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_W_ELM, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_W_ELM_"+V_charm_1_L.Tag[i_ens], UseJack));
  agm2_charm_SD_ELM_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_SD_ELM, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_SD_ELM_"+V_charm_1_L.Tag[i_ens], UseJack));
  //####################################################################################################



  // ######################################## OSTERWALDER-SEILER #######################################

  //L
  distr_t agm2_W_OS_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W_OS to zero by default
  distr_t agm2_SD_OS_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_OS to zero by default
  distr_t agm2_W_ELM_OS_L(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM_OS to zero by default
  distr_t agm2_SD_ELM_OS_L(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM_OS to zero by default
  //M
  distr_t agm2_W_OS_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W_OS to zero by default
  distr_t agm2_SD_OS_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_OS to zero by default
  distr_t agm2_W_ELM_OS_M(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM_OS to zero by default
  distr_t agm2_SD_ELM_OS_M(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM_OS to zero by default
  //H
  distr_t agm2_W_OS_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W_OS to zero by default
  distr_t agm2_SD_OS_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_OS to zero by default
  distr_t agm2_W_ELM_OS_H(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM_OS to zero by default
  distr_t agm2_SD_ELM_OS_H(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM_OS to zero by default


  //#################################################################################################

  for(int t=1; t< Corr.Nt/2; t++) {
    //L
    agm2_W_OS_L = agm2_W_OS_L + 4.0*pow(alpha,2)*Zv_L*Zv_L*V_charm_OS_distr_L.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_OS_L = agm2_SD_OS_L + 4.0*pow(alpha,2)*Zv_L*Zv_L*V_charm_OS_distr_L_pert_sub.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_OS_L = agm2_W_ELM_OS_L + 4.0*pow(alpha,2)*Zv_L*Zv_L*V_charm_OS_distr_L.distr_list[t]*Ker_ELM_OS_L.distr_list[t]*( distr_t::f_of_distr(th0, t*MV_charm_OS_L/m_Jpsi) - distr_t::f_of_distr(th1, t*MV_charm_OS_L/m_Jpsi));
    agm2_SD_ELM_OS_L = agm2_SD_ELM_OS_L + 4.0*pow(alpha,2)*Zv_L*Zv_L*V_charm_OS_distr_L_pert_sub.distr_list[t]*Ker_ELM_OS_L.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*MV_charm_OS_L/m_Jpsi));
    //M
    agm2_W_OS_M = agm2_W_OS_M + 4.0*pow(alpha,2)*Zv_M*Zv_M*V_charm_OS_distr_M.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_OS_M = agm2_SD_OS_M + 4.0*pow(alpha,2)*Zv_M*Zv_M*V_charm_OS_distr_M_pert_sub.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_OS_M = agm2_W_ELM_OS_M + 4.0*pow(alpha,2)*Zv_M*Zv_M*V_charm_OS_distr_M.distr_list[t]*Ker_ELM_OS_M.distr_list[t]*( distr_t::f_of_distr(th0, t*MV_charm_OS_M/m_Jpsi) - distr_t::f_of_distr(th1, t*MV_charm_OS_M/m_Jpsi));
    agm2_SD_ELM_OS_M = agm2_SD_ELM_OS_M + 4.0*pow(alpha,2)*Zv_M*Zv_M*V_charm_OS_distr_M_pert_sub.distr_list[t]*Ker_ELM_OS_M.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*MV_charm_OS_M/m_Jpsi));
    //H
    agm2_W_OS_H = agm2_W_OS_H + 4.0*pow(alpha,2)*Zv_H*Zv_H*V_charm_OS_distr_H.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_OS_H = agm2_SD_OS_H + 4.0*pow(alpha,2)*Zv_H*Zv_H*V_charm_OS_distr_H_pert_sub.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_OS_H = agm2_W_ELM_OS_H + 4.0*pow(alpha,2)*Zv_H*Zv_H*V_charm_OS_distr_H.distr_list[t]*Ker_ELM_OS_H.distr_list[t]*( distr_t::f_of_distr(th0, t*MV_charm_OS_H/m_Jpsi) - distr_t::f_of_distr(th1, t*MV_charm_OS_H/m_Jpsi));
    agm2_SD_ELM_OS_H = agm2_SD_ELM_OS_H + 4.0*pow(alpha,2)*Zv_H*Zv_H*V_charm_OS_distr_H_pert_sub.distr_list[t]*Ker_ELM_OS_H.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*MV_charm_OS_H/m_Jpsi));
  }
  
  //push_back the result
  //L
  agm2_charm_W_OS_L.distr_list.push_back(agm2_W_OS_L);
  agm2_charm_SD_OS_L.distr_list.push_back(agm2_SD_OS_L);
  agm2_charm_W_ELM_OS_L.distr_list.push_back(agm2_W_ELM_OS_L);
  agm2_charm_SD_ELM_OS_L.distr_list.push_back(agm2_SD_ELM_OS_L);
  //M
  agm2_charm_W_OS_M.distr_list.push_back(agm2_W_OS_M);
  agm2_charm_SD_OS_M.distr_list.push_back(agm2_SD_OS_M);
  agm2_charm_W_ELM_OS_M.distr_list.push_back(agm2_W_ELM_OS_M);
  agm2_charm_SD_ELM_OS_M.distr_list.push_back(agm2_SD_ELM_OS_M);
  //H
  agm2_charm_W_OS_H.distr_list.push_back(agm2_W_OS_H);
  agm2_charm_SD_OS_H.distr_list.push_back(agm2_SD_OS_H);
  agm2_charm_W_ELM_OS_H.distr_list.push_back(agm2_W_ELM_OS_H);
  agm2_charm_SD_ELM_OS_H.distr_list.push_back(agm2_SD_ELM_OS_H);


  
  //extrapolate the result to the physical D-meson point
  vector<distr_t> agm2s_charm_W_OS({agm2_W_OS_L, agm2_W_OS_M, agm2_W_OS_H});
  vector<distr_t> agm2s_charm_SD_OS({agm2_SD_OS_L, agm2_SD_OS_M, agm2_SD_OS_H});
  vector<distr_t> agm2s_charm_W_ELM_OS({agm2_W_ELM_OS_L, agm2_W_ELM_OS_M, agm2_W_ELM_OS_H});
  vector<distr_t> agm2s_charm_SD_ELM_OS({agm2_SD_ELM_OS_L, agm2_SD_ELM_OS_M, agm2_SD_ELM_OS_H});

  
  agm2_charm_W_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_W_OS, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_W_OS_"+V_charm_1_L.Tag[i_ens], UseJack));
  agm2_charm_SD_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_SD_OS, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_SD_OS_"+V_charm_1_L.Tag[i_ens], UseJack));
  agm2_charm_W_ELM_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_W_ELM_OS,X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_W_ELM_OS_"+V_charm_1_L.Tag[i_ens], UseJack));
  agm2_charm_SD_ELM_OS_Extr.distr_list.push_back( Obs_extrapolation_meson_mass(agm2s_charm_SD_ELM_OS, X_2_fit, X_2_phys, "../data/gm2/charm", "agm2_SD_ELM_OS_"+V_charm_1_L.Tag[i_ens], UseJack));

  //####################################################################################################

  
  //print to file
  //L
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_L.ave(), agm2_distr_Tdata_L.err()}, "../data/gm2/charm/tm/agm2_Tdata_"+V_charm_1_L.Tag[i_ens]+"_ELM_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_L.ave(), agm2_OS_distr_Tdata_L.err()}, "../data/gm2/charm/OS/agm2_Tdata_"+V_charm_1_L.Tag[i_ens]+"_ELM_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
  //M
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_M.ave(), agm2_distr_Tdata_M.err()}, "../data/gm2/charm/tm/agm2_Tdata_"+V_charm_1_M.Tag[i_ens]+"_ELM_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_M.ave(), agm2_OS_distr_Tdata_M.err()}, "../data/gm2/charm/OS/agm2_Tdata_"+V_charm_1_M.Tag[i_ens]+"_ELM_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
  //H
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_H.ave(), agm2_distr_Tdata_H.err()}, "../data/gm2/charm/tm/agm2_Tdata_"+V_charm_1_H.Tag[i_ens]+"_ELM_H.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_H.ave(), agm2_OS_distr_Tdata_H.err()}, "../data/gm2/charm/OS/agm2_Tdata_"+V_charm_1_H.Tag[i_ens]+"_ELM_H.dat.t", "", "#id  Tdata   ag2m agm2_err");

  //NON_ELM
  //L
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_No_ELM_L.ave(), agm2_distr_Tdata_No_ELM_L.err()}, "../data/gm2/charm/tm/agm2_Tdata_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_No_ELM_L.ave(), agm2_OS_distr_Tdata_No_ELM_L.err()}, "../data/gm2/charm/OS/agm2_Tdata_"+V_charm_1_L.Tag[i_ens]+"_L.dat.t", "", "#id  Tdata   ag2m agm2_err");
  //M
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_No_ELM_M.ave(), agm2_distr_Tdata_No_ELM_M.err()}, "../data/gm2/charm/tm/agm2_Tdata_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_No_ELM_M.ave(), agm2_OS_distr_Tdata_No_ELM_M.err()}, "../data/gm2/charm/OS/agm2_Tdata_"+V_charm_1_M.Tag[i_ens]+"_M.dat.t", "", "#id  Tdata   ag2m agm2_err");
  //H
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata_No_ELM_H.ave(), agm2_distr_Tdata_No_ELM_H.err()}, "../data/gm2/charm/tm/agm2_Tdata_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_OS_distr_Tdata_No_ELM_H.ave(), agm2_OS_distr_Tdata_No_ELM_H.err()}, "../data/gm2/charm/OS/agm2_Tdata_"+V_charm_1_H.Tag[i_ens]+"_H.dat.t", "", "#id  Tdata   ag2m agm2_err");
  
  }


  cout<<"charm quark correlator analyzed!"<<endl;

  //light
  channel="l";
  

 
  
  for(int i_ens=0;i_ens<Nens_light;i_ens++) {
    
  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = V_light_1.nrows[i_ens];

  //resample lattice spacing
  distr_t a_distr(UseJack);
  LatticeInfo L_info;
  L_info.LatInfo_new_ens(V_light_1.Tag[i_ens]);
 
   
  distr_t_list  V_light_1_distr, V_light_2_distr, V_light_3_distr, disco_distr;
  distr_t_list  V_light_OS_1_distr, V_light_OS_2_distr, V_light_OS_3_distr, V_light_OS_distr;  
  distr_t_list  V_light_distr, MV_light_distr, ZV_light_distr;
  distr_t_list  MV_light_OS_distr, ZV_light_OS_distr;
  distr_t_list Mpi_distr, Mpi_OS_distr, fp_distr, overlap_P5P5_distr, overlap_P5P5_OS_distr, pion_corr;
  distr_t_list A0P5_distr, A0P5_OS_distr;
  distr_t MV_light, ZV_light, Mpi, Mpi_OS, fp;
  distr_t MV_light_OS, ZV_light_OS;
  distr_t_list P5P5_OS_distr,  RV, RA,RA0, ratio_P5P5_overlap_OS_tm, Zp_ov_Zs_distr;
  distr_t Zv, Za, Zp_ov_Zs, fp_ov_Mpi;
  

  //Analyze correlators
  if(V_light_1.Tag[i_ens].substr(1,1)=="A") {Corr.Tmin=14;Corr.Tmax=19; a_distr=a_A;}
  else if(V_light_1.Tag[i_ens].substr(1,12)=="B211b.072.96") {Corr.Tmin=26; Corr.Tmax=60; a_distr=a_B;}
  else if(V_light_1.Tag[i_ens].substr(1,12)=="B211b.072.64") { Corr.Tmin=24; Corr.Tmax=36; a_distr = a_B;}
  else if(V_light_1.Tag[i_ens].substr(1,11)=="B211a.25.24") {Corr.Tmin=15; Corr.Tmax=20; a_distr=a_B;}
  else if(V_light_1.Tag[i_ens].substr(1,11)=="B211a.25.32") {Corr.Tmin=21; Corr.Tmax=30; a_distr = a_B;}
  else if(V_light_1.Tag[i_ens].substr(1,11)=="B211a.25.48") {Corr.Tmin=21; Corr.Tmax=40; a_distr = a_B;}
  else if(V_light_1.Tag[i_ens].substr(1,1)=="C") {Corr.Tmin=40; Corr.Tmax=60; a_distr=a_C;}
  else if(V_light_1.Tag[i_ens].substr(1,1)=="D") {Corr.Tmin=41; Corr.Tmax=80; a_distr=a_D;}
  else crash("Cannot find ensemble: "+V_light_1.Tag[i_ens]+" while modifying fit interval for MV_light");



  //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR OPPOSITE R ####################################
  string Pt_free_oppor= "../free_theory_VV/"+to_string(Corr.Nt/2)+"_OPPOR";
  Vfloat VV_free_oppor= Read_From_File(Pt_free_oppor, 1, 2);
  if(VV_free_oppor.size() != Corr.Nt) crash("Failed to read properly free VV correlator w opposite r");
  //for(auto &el: VV_free_oppor) el*=-1;
  //################## READ FREE THEORY VECTOR-VECTOR CORRELATOR SAME R ####################################
  string Pt_free_samer= "../free_theory_VV/"+to_string(Corr.Nt/2)+"_SAMER";
  Vfloat VV_free_samer= Read_From_File(Pt_free_samer, 1, 2);
  if(VV_free_samer.size() != Corr.Nt) crash("Failed to read properly free VV correlator w same r");


  //define lambda function to compute function of distributions
  auto sqr= [=](double a, double b) {return sqrt(a);};
  auto SINH= [](double m) -> double  {return sinh(m);};
  auto SINH2 = [](double m, double t) -> double {return sinh(m);};

  
  //vector light sector
  V_light_1_distr = Corr.corr_t(V_light_1.col(0)[i_ens], "../data/gm2/light/tm/corr_1_"+V_light_1.Tag[i_ens]+".dat");
  V_light_2_distr = Corr.corr_t(V_light_2.col(0)[i_ens], "../data/gm2/light/tm/corr_2_"+V_light_2.Tag[i_ens]+".dat");
  V_light_3_distr = Corr.corr_t(V_light_3.col(0)[i_ens], "../data/gm2/light/tm/corr_3_"+V_light_3.Tag[i_ens]+".dat");
  V_light_OS_1_distr = Corr.corr_t(V_light_OS_1.col(0)[i_ens], "../data/gm2/light/OS/corr_1_"+V_light_1.Tag[i_ens]+".dat");
  V_light_OS_2_distr = Corr.corr_t(V_light_OS_2.col(0)[i_ens], "../data/gm2/light/OS/corr_2_"+V_light_2.Tag[i_ens]+".dat");
  V_light_OS_3_distr = Corr.corr_t(V_light_OS_3.col(0)[i_ens], "../data/gm2/light/OS/corr_3_"+V_light_3.Tag[i_ens]+".dat");
  //sum over the Lorentz indices of the e.m. current
  V_light_distr= ((pow(qu,2)+pow(qd,2))/3.0)*(V_light_1_distr+ V_light_2_distr + V_light_3_distr);
  V_light_OS_distr =  ((pow(qu,2)+pow(qd,2))/3.0)*(V_light_OS_1_distr+ V_light_OS_2_distr + V_light_OS_3_distr);
  bool Found_disco_ens=false;
  if(Include_light_disco) {
    int i_ens_disco=0;
    for(int j=0;j<Nens_disco_light;j++) if(disco_light.Tag[j] == V_light_1.Tag[i_ens]) { Found_disco_ens=true; i_ens_disco=j;disco_Tags.push_back(disco_light.Tag[j]) ;break;}
    if(Found_disco_ens) {
    disco_distr = Corr.corr_t(disco_light.col(0)[i_ens_disco], "");
    disco_distr = disco_distr*(pow(qu+qd,2));
    }
  }
  MV_light_distr= Corr.effective_mass_t(V_light_distr, "../data/gm2/light/tm/MV_"+V_light_1.Tag[i_ens]+".dat");
  ZV_light_distr= Corr.residue_t(V_light_distr, "../data/gm2/light/tm/ZV_overlap_"+V_light_1.Tag[i_ens]+"dat");
  MV_light_OS_distr= Corr.effective_mass_t(V_light_OS_distr, "../data/gm2/light/OS/MV_"+V_light_1.Tag[i_ens]+".dat");
  ZV_light_OS_distr= Corr.residue_t(V_light_OS_distr, "../data/gm2/light/OS/ZV_overlap_"+V_light_1.Tag[i_ens]+"dat");
  


  //tm pion sector
  pion_corr= Corr.corr_t(pt2_pion.col(0)[i_ens], "");
  Mpi_distr= Corr.effective_mass_t(pion_corr, "../data/gm2/light/Mpi_"+V_light_1.Tag[i_ens]+".dat");
  overlap_P5P5_distr = Corr.residue_t(pion_corr, "");
  fp_distr = Corr.decay_constant_t(pow(2.0*L_info.ml,2)*pion_corr, "../data/gm2/light/fp_"+V_light_1.Tag[i_ens]+".dat");


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
  else if(V_light_1.Tag[i_ens].substr(1,1)=="B" || V_light_1.Tag[i_ens].substr(1,1) =="A") {Corr.Tmin=16; Corr.Tmax=22;}
  else crash("Cannot find ensemble: "+V_light_1.Tag[i_ens]+" while modifying fit interval for MV_light");
  //################################################

  
  MV_light = Corr.Fit_distr(MV_light_distr);
  ZV_light = Corr.Fit_distr(ZV_light_distr);
  MV_light_OS = Corr.Fit_distr(MV_light_OS_distr);
  ZV_light_OS = Corr.Fit_distr(ZV_light_OS_distr);
  //push_back MV_light and ZV_light
  MV_fit_light.distr_list.push_back(MV_light);
  ZV_fit_light.distr_list.push_back(ZV_light);
  MV_fit_light_OS.distr_list.push_back(MV_light_OS);
  ZV_fit_light_OS.distr_list.push_back(ZV_light_OS);


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

  //set plateaux for RA
  
  //Analyze correlators
  if(V_light_1.Tag[i_ens].substr(1,1)=="A") {Corr.Tmin=14;Corr.Tmax=19; a_distr=a_A;}
  else if(V_light_1.Tag[i_ens].substr(1,12)=="B211b.072.96") {Corr.Tmin=12; Corr.Tmax=35; a_distr=a_B;}
  else if(V_light_1.Tag[i_ens].substr(1,12)=="B211b.072.64") { Corr.Tmin=15; Corr.Tmax=38; a_distr = a_B;}
  else if(V_light_1.Tag[i_ens].substr(1,11)=="B211a.25.24") {Corr.Tmin=12; Corr.Tmax=22; a_distr=a_B;}
  else if(V_light_1.Tag[i_ens].substr(1,11)=="B211a.25.32") {Corr.Tmin=19; Corr.Tmax=29; a_distr = a_B;}
  else if(V_light_1.Tag[i_ens].substr(1,11)=="B211a.25.48") {Corr.Tmin=16; Corr.Tmax=32; a_distr = a_B;}
  else if(V_light_1.Tag[i_ens].substr(1,1)=="C") {Corr.Tmin=21; Corr.Tmax=43; a_distr=a_C;}
  else if(V_light_1.Tag[i_ens].substr(1,1)=="D") {Corr.Tmin=28; Corr.Tmax=51; a_distr=a_D;}
  else crash("Cannot find ensemble: "+V_light_1.Tag[i_ens]+" while modifying fit interval for MV_light");
  Za = Corr.Fit_distr(RA);

  a_distr_list.distr_list.push_back(a_distr);
  //#################################################################

  Vfloat free_corr_log_art(Corr.Nt);
  for(int t=0;t<Corr.Nt;t++) {free_corr_log_art[t] = -1.0*(qu*qu +qd*qd)*(t !=0)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));}

  distr_t_list V_light_distr_tm_corr = (V_light_distr); //free_corr_log_art
  distr_t_list V_light_distr_OS_corr = (V_light_OS_distr);
  //Print To File additional observables
  // print summed connected correlators to file
  Print_To_File({}, {V_light_distr.ave(), V_light_distr.err(), (Za*Za*V_light_distr).ave(), (Za*Za*V_light_distr).err(),  ( (Za*Za*V_light_distr- Zv*Zv*V_light_OS_distr)/(Za*Za*V_light_distr)).ave(),  ( (Za*Za*V_light_distr- Zv*Zv*V_light_OS_distr)/(Za*Za*V_light_distr)).err(), (Za*Za*V_light_distr_tm_corr).ave(), (Za*Za*V_light_distr_tm_corr).err()}, "../data/gm2/light/tm/corr_sum_"+V_light_1.Tag[i_ens]+".dat.t", "", "#time   V(t)^tm    V(t)^tm(renormalized)    DV(t)(tm-os/tm)  V(t)^tm(VV_free corrected) ");
  Print_To_File({}, {V_light_OS_distr.ave(), V_light_OS_distr.err(), (Zv*Zv*V_light_OS_distr).ave(), (Zv*Zv*V_light_OS_distr).err(), (Zv*Zv*V_light_distr_OS_corr).ave() , (Zv*Zv*V_light_distr_OS_corr).err()}, "../data/gm2/light/OS/corr_sum_"+V_light_1.Tag[i_ens]+".dat.t", "", "#time V(t)^OS   V(t)^OS(renormalized) V(t)^OS(VV_free_corrected)");

  if(Include_light_disco && Found_disco_ens) Print_To_File({}, {(Zv*Zv*disco_distr).ave(), (Zv*Zv*disco_distr).err()}, "../data/gm2/light/disco/disc_"+V_light_1.Tag[i_ens]+".dat.t","",""); 
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
  cout<<"lat spacing from Mpi: "<<Mpi.ave()/(m_pi*fm_to_inv_Gev)<<" +- "<<Mpi.err()/(m_pi*fm_to_inv_Gev)<<" (fm)"<<endl;
  cout<<"lat spacing from fp: "<<fp.ave()/(fp_phys*fm_to_inv_Gev)<<" +- "<<fp.err()/(fp_phys*fm_to_inv_Gev)<<" (fm)"<<endl;
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

  //push_back Mpi, Mpi_OS, fp, Zv a and L if Include_light_disco && Found_disco
  if(Include_light_disco && Found_disco_ens) {
    Mpi_fit_disco.distr_list.push_back(Mpi);
    Mpi_OS_fit_disco.distr_list.push_back(Mpi_OS);
    fp_fit_disco.distr_list.push_back(fp);
    Zv_fit_disco.distr_list.push_back(Zv);
    ml_list_disco.push_back(L_info.ml);
    a_list_disco.push_back(L_info.a);
    L_list_disco.push_back((double)L_info.L);
  }
  

  


  int Tdata_min= 8;
  int Tdata_max = Corr.Nt/2.0 -4;
  int Tdata_fit = (Corr.Tmax+Corr.Tmin)/2;
  
  //compute kernel distribution
  distr_t_list Kernel_distr = distr_t_list::f_of_distr(K,Mpi/m_pi, Upper_Limit_Time_Integral_light+1);
  //compute exp(-Mv*t) distribution
  distr_t_list exp_MVl = distr_t_list::f_of_distr(exp_MV, MV_light, Upper_Limit_Time_Integral_light+1);
  distr_t_list exp_MVl_OS = distr_t_list::f_of_distr(exp_MV, MV_light_OS, Upper_Limit_Time_Integral_light+1);


  //Print single-exponential prediction to file
  Print_To_File({}, {(exp_MVl*(ZV_light/(2.0*MV_light))).ave(), (exp_MVl*(ZV_light/(2.0*MV_light))).err()}, "../data/gm2/light/tm/corr_gsd_sum_"+V_light_1.Tag[i_ens]+".dat.t", "", "");
  Print_To_File({}, {(exp_MVl_OS*(ZV_light_OS/(2.0*MV_light_OS))).ave(), (exp_MVl_OS*(ZV_light_OS/(2.0*MV_light_OS))).err()}, "../data/gm2/light/OS/corr_gsd_sum_"+V_light_1.Tag[i_ens]+".dat.t", "", "");

  distr_t_list agm2_distr_Tdata(UseJack);
  distr_t_list agm2_distr_OS_Tdata(UseJack);
  Vfloat Tdata_vec;

  
  for(int Tdata=Tdata_min;Tdata<Tdata_max;Tdata++) {
    //compute 4\pia^2 using lattice data up to Tdata (included)
   
    distr_t agm2(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2 to zero by default
    distr_t agm2_OS(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_OS to zero by default 
    
    
    for(int t=1;t<=Upper_Limit_Time_Integral_light;t++) {
      if(t<=Tdata) {
	agm2 = agm2 + 4.0*pow(alpha,2)*V_light_distr.distr_list[t]*Kernel_distr.distr_list[t];
	agm2_OS = agm2_OS + 4.0*pow(alpha,2)*V_light_OS_distr.distr_list[t]*Kernel_distr.distr_list[t];
      }
      else {
	agm2= agm2 + 4.0*pow(alpha,2)*(ZV_light/(2.0*MV_light))*exp_MVl.distr_list[t]*Kernel_distr.distr_list[t];
	agm2_OS= agm2_OS + 4.0*pow(alpha,2)*(ZV_light_OS/(2.0*MV_light_OS))*exp_MVl_OS.distr_list[t]*Kernel_distr.distr_list[t];
      }
    }
    Tdata_vec.push_back((double)Tdata);
    agm2 = agm2*Za*Za;
    agm2_OS = agm2_OS*Zv*Zv;
    agm2_distr_Tdata.distr_list.push_back(agm2);
    agm2_distr_OS_Tdata.distr_list.push_back(agm2_OS);
    if(Tdata==Tdata_fit) {
      agm2_light.distr_list.push_back(agm2);
      agm2_light_OS.distr_list.push_back(agm2_OS);
    }
  }
  //print to file
  Print_To_File({}, {Tdata_vec, agm2_distr_Tdata.ave(), agm2_distr_Tdata.err()}, "../data/gm2/light/tm/agm2_Tdata_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m agm2_err");
  Print_To_File({}, {Tdata_vec, agm2_distr_OS_Tdata.ave(), agm2_distr_OS_Tdata.err()}, "../data/gm2/light/OS/agm2_Tdata_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id  Tdata   ag2m agm2_err");
  
  
 

  //#######################  INTERMEDIATE AND SHORT-DISTANCE WINDOW ###################################

    
  //############################   TWISTED MASS ######################################################

  distr_t agm2_W(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W to zero by default
  distr_t agm2_SD(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD to zero by default
  distr_t agm2_W_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM to zero by default
  distr_t agm2_SD_ELM(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM to zero by default


  //#################################################################################################
    
  distr_t X = Mpi/( (Mpi/a_distr).ave()); //variable to be used in ELT

  distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);
  distr_t_list Ker_ELM = distr_t_list::f_of_distr(K, X, Corr.Nt/2);
    
  //define lambdas for the theta func
  auto th0 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t0)/Delta));};
  auto th1 = [](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t1)/Delta));};
  for(int t=1; t< Corr.Nt/2; t++) {
    agm2_W = agm2_W + 4.0*pow(alpha,2)*Za*Za*V_light_distr_tm_corr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD = agm2_SD + 4.0*pow(alpha,2)*Za*Za*(V_light_distr_tm_corr.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM = agm2_W_ELM + 4.0*pow(alpha,2)*Za*Za*V_light_distr_tm_corr.distr_list[t]*Ker_ELM.distr_list[t]*( distr_t::f_of_distr(th0, t*X) - distr_t::f_of_distr(th1, t*X));
    agm2_SD_ELM = agm2_SD_ELM + 4.0*pow(alpha,2)*(Za*Za*V_light_distr_tm_corr.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
  }
  
  //push_back the result

  agm2_light_W.distr_list.push_back(agm2_W);
  agm2_light_SD.distr_list.push_back(agm2_SD);
  agm2_light_W_ELM.distr_list.push_back(agm2_W_ELM);
  agm2_light_SD_ELM.distr_list.push_back(agm2_SD_ELM);

  //####################################################################################################



  // ######################################## OSTERWALDER-SEILER #######################################


  distr_t agm2_W_OS(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_W_OS to zero by default
  distr_t agm2_SD_OS(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_OS to zero by default
  distr_t agm2_W_ELM_OS(UseJack, UseJack?Njacks:Nboots); //constructor sets agm2_W_ELM_OS to zero by default
  distr_t agm2_SD_ELM_OS(UseJack, UseJack?Njacks:Nboots);  //constructor sets agm2_SD_ELM_OS to zero by default


  //#################################################################################################

  for(int t=1; t< Corr.Nt/2; t++) {
    agm2_W_OS = agm2_W_OS + 4.0*pow(alpha,2)*Zv*Zv*V_light_distr_OS_corr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_SD_OS = agm2_SD_OS + 4.0*pow(alpha,2)*Zv*Zv*(V_light_distr_OS_corr.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_W_ELM_OS = agm2_W_ELM_OS + 4.0*pow(alpha,2)*Zv*Zv*V_light_distr_OS_corr.distr_list[t]*Ker_ELM.distr_list[t]*( distr_t::f_of_distr(th0, t*X) - distr_t::f_of_distr(th1, t*X));
    agm2_SD_ELM_OS = agm2_SD_ELM_OS + 4.0*pow(alpha,2)*(Zv*Zv*V_light_distr_OS_corr.distr_list[t])*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
  }
  
  //push_back the result

  agm2_light_W_OS.distr_list.push_back(agm2_W_OS);
  agm2_light_SD_OS.distr_list.push_back(agm2_SD_OS);
  agm2_light_W_ELM_OS.distr_list.push_back(agm2_W_ELM_OS);
  agm2_light_SD_ELM_OS.distr_list.push_back(agm2_SD_ELM_OS);

  //####################################################################################################




  //####################################### DISCO LIGHT ###############################################
  distr_t agm2_disco_W(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
  distr_t agm2_disco_W_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
  distr_t agm2_disco_SD(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0
  distr_t agm2_disco_SD_ELM(UseJack, UseJack?Njacks:Nboots); //constructor sets it to 0

  if(Include_light_disco && Found_disco_ens) {
    
    for(int t=1; t< Corr.Nt/2; t++) {
    agm2_disco_W = agm2_disco_W + 4.0*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr) - distr_t::f_of_distr(th1, t*a_distr));
    agm2_disco_SD = agm2_disco_SD + 4.0*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*a_distr));
    agm2_disco_W_ELM= agm2_disco_W_ELM + 4.0*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker_ELM.distr_list[t]*( distr_t::f_of_distr(th0, t*X) - distr_t::f_of_distr(th1, t*X));
    agm2_disco_SD_ELM= agm2_disco_SD_ELM + 4.0*pow(alpha,2)*Zv*Zv*disco_distr.distr_list[t]*Ker_ELM.distr_list[t]*( 1.0 - distr_t::f_of_distr(th0, t*X));
    }
  
  //push_back the result

  agm2_disco_light_W.distr_list.push_back(agm2_disco_W);
  agm2_disco_light_SD.distr_list.push_back(agm2_disco_SD);
  agm2_disco_light_W_ELM.distr_list.push_back(agm2_disco_W_ELM);
  agm2_disco_light_SD_ELM.distr_list.push_back(agm2_disco_SD_ELM);

  }



  
  










  


  
  //#################################################################################
  //#################################################################################
  //#################################################################################
  //#################################################################################
  //fit lattice data using analytic representation for V(t)_light


  //#########################################  TWISTED MASS CORRELATOR ###################################################################################################


  if(!Skip_total_light_calc) {

  bootstrap_fit<fit_par,ipar> bf(Njacks);
  bf.set_warmup_lev(4); //sets warmup

  //######################FIT INTERVALS####################
  //##############INCLUDE only times t for which dC(t)/C(t) < 0.1#####################

  //Tfit_min always starts at 0.2fm
  //double t02fm= 0.2/L_info.a +1.0;
  int Tfit_min= 5;
  if(V_light_1.Tag[i_ens].substr(1,1) == "D") Tfit_min++;
  int Tfit_max= Tfit_min;  
  bool Found_error_less_10_percent=false;
  while(!Found_error_less_10_percent && Tfit_max < Corr.Nt/2) {
   
    if( (Za*Za*V_light_distr.distr_list[Tfit_max]).err()/fabs( (Za*Za*V_light_distr.distr_list[Tfit_max]).ave()) <  0.1) Tfit_max++;
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
    bf.Set_limits("Mrho", 2.3 , 7.0);
  }
  else {
    bf.Add_par("Ed", 1.3, 0.2);
    bf.Set_limits("Ed", 0.5, 3.0);
    bf.Add_par("Mrho",3.0, 0.05);
    bf.Set_limits("Mrho", 2.4, 5.0);
  }
  bf.Add_par("gpi", 1.0, 0.01);
  bf.Add_prior_par("kappa", -3.0, 0.1);
  bf.Set_limits("gpi",0.6, 1.5);
  //bf.Set_limits("kappa", -5.0, 1.0);


  //bf.Fix_par("gpi",1.0);
  bf.Fix_par("kappa", 0.0);

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
		     return ip.V_light;
		   };
  bf.error =  [&](const fit_par& p,const ipar &ip) -> double {
		return ip.V_light_err;
	      };

  //fill the data
  vector<vector<ipar>> data(Njacks);
  //allocate space for output result
  boot_fit_data<fit_par> Bt_fit;

  

  for(auto &data_iboot: data) data_iboot.resize(Tfit_points);
  
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int t=Tfit_min;t<Tfit_points+Tfit_min;t++) {
      int tt=t-Tfit_min;
      data[ijack][tt].V_light = pow(Za.distr[ijack],2)*V_light_distr.distr_list[t].distr[ijack];
      data[ijack][tt].V_light_err = (Za*Za*V_light_distr.distr_list[t]).err();
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
    //bf.Append_to_prior("gpi", 1.0, 0.2);
    bf.Append_to_prior("kappa", 0.0, 3.0);
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
      func.distr_list[istep].distr.push_back( bf.ansatz( Bt_fit.par[ijack], my_ipar ));
    }
  }

  //print to file
  Print_To_File({}, {times, func.ave(), func.err()}, "../data/gm2/light/tm/Corr_anal_rep_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id_col  time  val   err  (fit interval: ["+to_string(Tfit_min)+","+to_string(Tfit_points+Tfit_min-1)+"])  a = "+to_string_with_precision(L_info.a,6)+" +- "+to_string_with_precision(L_info.a_err, 6)+" fm");
  //####################################################
    
    
  //push_back jack distribution of fitted params
  par_list_anal_repr.push_back(Bt_fit.par);
  cout<<"SIZE OF ENERGY LEVEL MAP FOR TM-FIT IS: "<<Energy_lev_list.size()<<endl;
  cout<<"#########################END ENSEMBLE FIT"<<endl;

  //##################################################################################################################################################################################################


  
  //#########################################################  OSTERWALDER-SEILER CORRELATOR  ########################################################################################################




  bootstrap_fit<fit_par,ipar> bf_OS(Njacks);
  bf_OS.set_warmup_lev(4); //sets warmup

  //######################FIT INTERVALS####################
  //##############INCLUDE only times t for which dC(t)/C(t) < 0.1#####################

  //Tfit_min always starts at 0.2fm
  //double t02fm= 0.2/L_info.a +1.0;
  Tfit_min= 6;
  if(V_light_1.Tag[i_ens].substr(1,1) == "D") Tfit_min++;
  Tfit_max= Tfit_min;  
  Found_error_less_10_percent=false;
  while(!Found_error_less_10_percent && Tfit_max < Corr.Nt/2) {
   
    if( (Zv*Zv*V_light_OS_distr.distr_list[Tfit_max]).err()/fabs( (Zv*Zv*V_light_OS_distr.distr_list[Tfit_max]).ave()) <  0.1) Tfit_max++;
    else Found_error_less_10_percent=true;
  }
  Tfit_points = Tfit_max +1 - Tfit_min;
  //###################################################################################

  
  bf_OS.Set_number_of_measurements(Tfit_points);
  bf_OS.Set_verbosity(verbosity);

  cout<<"FITTING WITH PIPI+DUAL REPRESENTATION ENSEMBLE  (OS correlator): "<<V_light_1.Tag[i_ens]<<endl;
  cout<<"Tmin: "<<Tfit_min<<", Tmax: "<<Tfit_max<<endl;

  //initial guesses are based on the results on the old ETMC confs
  bf_OS.Add_par("Rd", 1.1, 0.1);
  bf_OS.Set_limits("Rd", 0.7, 1.8);
  if(!Use_Mpi_OS) {
    bf_OS.Add_par("Ed", 2.5, 0.2);
    bf_OS.Set_limits("Ed", 0.7, 5.0);
    bf_OS.Add_par("Mrho",5.5, 0.1);
    bf_OS.Set_limits("Mrho", 2.3 , 7.0);
  }
  else {
    bf_OS.Add_par("Ed", 1.3, 0.2);
    bf_OS.Set_limits("Ed", 0.5, 3.0);
    bf_OS.Add_par("Mrho",3.0, 0.05);
    bf_OS.Set_limits("Mrho", 2.4, 5.0);
  }
  bf_OS.Add_par("gpi", 1.0, 0.01);
  bf_OS.Add_prior_par("kappa", -3.0, 0.1);
  bf_OS.Set_limits("gpi",0.6, 1.5);
  //bf.Set_limits("kappa", -5.0, 1.0);


  //bf.Fix_par("gpi",1.0);
  bf_OS.Fix_par("kappa", 0.0);

  map<pair<pair<double,double>, pair<double,double>>,Vfloat> Energy_lev_list_OS;
  

  bf_OS.ansatz =  [&](const fit_par &p, const ipar &ip) -> double {


		 double Pi_M = Use_Mpi_OS?ip.Mp_OS:ip.Mp;

		 double GPI = p.gpi*p.Mrho*Pi_M/ip.fp;
		 Vfloat Knpp;
		 pair<double,double> Mass_par= make_pair(p.Mrho*Pi_M, Pi_M);
		 pair<double,double> Couplings = make_pair(GPI, p.kappa);
		 pair< pair<double,double>,pair<double, double>> input_pars = make_pair( Mass_par, Couplings);
		 map<pair<pair<double,double>, pair<double, double>>,Vfloat>::iterator it;
		 it= Energy_lev_list_OS.find(input_pars);
		 if(it != Energy_lev_list_OS.end()) Knpp= it->second;
		 else {
		   LL.Find_pipi_energy_lev(ip.L,p.Mrho*Pi_M, GPI, Pi_M, p.kappa, Knpp);
		   //add to Energy_lev_list
		   Energy_lev_list_OS.insert( make_pair(input_pars, Knpp));
		 }
		
		 return Qfact*LL.V_pipi(ip.t, ip.L, p.Mrho*Pi_M, GPI, Pi_M, p.kappa, Knpp) + LL.Vdual(ip.t, p.Mrho*Pi_M, p.Ed*Pi_M, p.Rd);
	       };
  bf_OS.measurement = [&](const fit_par& p,const ipar& ip) -> double {
		     return ip.V_light;
		   };
  bf_OS.error =  [&](const fit_par& p,const ipar &ip) -> double {
		return ip.V_light_err;
	      };

  //fill the data
  vector<vector<ipar>> data_OS(Njacks);
  //allocate space for output result
  boot_fit_data<fit_par> Bt_fit_OS;

  

  for(auto &data_iboot: data_OS) data_iboot.resize(Tfit_points);
  
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int t=Tfit_min;t<Tfit_points+Tfit_min;t++) {
      int tt=t-Tfit_min;
      data_OS[ijack][tt].V_light = pow(Zv.distr[ijack],2)*V_light_OS_distr.distr_list[t].distr[ijack];
      data_OS[ijack][tt].V_light_err = (Zv*Zv*V_light_OS_distr.distr_list[t]).err();
      data_OS[ijack][tt].L = L_info.L;
      data_OS[ijack][tt].Mp = Mpi.distr[ijack];
      data_OS[ijack][tt].Mp_OS = Mpi_OS.distr[ijack];
      data_OS[ijack][tt].t = t;
      data_OS[ijack][tt].fp = fp.distr[ijack];
    } 
  }

    
  //add prior values
  for(int ijack=0;ijack<Njacks;ijack++) {
    bf_OS.ib= &ijack;
    //bf.Append_to_prior("gpi", 1.0, 0.2);
    bf_OS.Append_to_prior("kappa", 0.0, 3.0);
  }

    
  //append
  bf_OS.Append_to_input_par(data_OS);

    
    
  //fit
  Bt_fit_OS= bf_OS.Perform_bootstrap_fit();
    
  Bt_fit_OS.ch2_ave();

  //###################################################
  //print fitted func

  cout<<"Printing fit function (OS) with: istep: "<<step_size_time<<" and nstep: "<<nsteps<<endl;
  distr_t_list func_OS(UseJack, nsteps);
  Vfloat times_OS;

  //generate data
  for(int istep=0; istep<nsteps;istep++) {
    ipar my_ipar;
    my_ipar.t = step_size_time*(1.0 + istep);
    my_ipar.L = (double)L_info.L;
    times_OS.push_back(my_ipar.t);
    for(int ijack=0;ijack<Njacks;ijack++) {
      my_ipar.Mp=Mpi.distr[ijack];
      my_ipar.Mp_OS=Mpi_OS.distr[ijack];
      my_ipar.fp = fp.distr[ijack];
      func_OS.distr_list[istep].distr.push_back( bf.ansatz( Bt_fit_OS.par[ijack], my_ipar ));
    }
  }

  //print to file
  Print_To_File({}, {times_OS, func_OS.ave(), func_OS.err()}, "../data/gm2/light/OS/Corr_anal_rep_"+V_light_1.Tag[i_ens]+".dat.t", "", "#id_col  time  val   err  (fit interval: ["+to_string(Tfit_min)+","+to_string(Tfit_points+Tfit_min-1)+"])  a = "+to_string_with_precision(L_info.a,6)+" +- "+to_string_with_precision(L_info.a_err, 6)+" fm");
  //####################################################
    
    
  //push_back jack distribution of fitted params
  par_list_anal_repr_OS.push_back(Bt_fit_OS.par);
  cout<<"SIZE OF ENERGY LEVEL MAP FOR OS-FIT IS: "<<Energy_lev_list_OS.size()<<endl;
  cout<<"#########################END ENSEMBLE FIT"<<endl;

  //#################################################################################################################################################################################################

  }
  }
  
  cout<<"light quark correlator analyzed!"<<endl;








  

  distr_t_list Edual(UseJack), Rdual(UseJack), Mrho(UseJack), gpi(UseJack), Kappa(UseJack);
  distr_t_list Edual_OS(UseJack), Rdual_OS(UseJack), Mrho_OS(UseJack), gpi_OS(UseJack), Kappa_OS(UseJack);

  if(!Skip_total_light_calc) {





  cout<<"####### Reconstructing agm2_light using fit paramters Ed, Rd, Mrho, gpi....."<<endl;


  //#############################################   TWISTED MASS CORRELATOR ##########################################################################

  
  cout<<"#########################  TWISTED MASS CORRELATOR ###########################"<<endl;



  
 
  for(int i_ens=0; i_ens<Nens_light;i_ens++) {
    cout<<"######### ENSEMBLE: "<<V_light_1.Tag[i_ens]<<endl;
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_light_1.Tag[i_ens]);
  

    distr_t a_distr = a_distr_list[i_ens];
    
    
    distr_t Ed(UseJack), Rd(UseJack), Mr(UseJack), g(UseJack), kap(UseJack), agm2_dual(UseJack), mp(UseJack), agm2_2L_dual(UseJack), agm2_infL_dual(UseJack), agm2_Lprime_dual(UseJack);
    for(int ijack=0;ijack<Njacks;ijack++) {
      //retrieve fit_parameters
      fit_par my_fit_par = par_list_anal_repr[i_ens][ijack];
      Ed.distr.push_back(my_fit_par.Ed);
      Rd.distr.push_back(my_fit_par.Rd);
      Mr.distr.push_back(my_fit_par.Mrho);
      g.distr.push_back(my_fit_par.gpi);
      kap.distr.push_back(my_fit_par.kappa);
       

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
		  
		    //double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
		    double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		    double func_val =  Qfact*LL.V_pipi(t, L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa, Elev)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

		    double F_int_val = kern_val*func_val;
        

		    return F_int_val;
		  };


      auto F_int_2L = [&](double t) {

			//double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
			//double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
			double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
			double func_val =  Qfact*LL.V_pipi(t, 1.5*L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa, Elev_2L)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

			double F_int_val = kern_val*func_val;
        

			return F_int_val;
		      };

      auto F_int_Lprime = [&](double t) {

			    //  double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
			    //double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
			    double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
			    double func_val =  Qfact*LL.V_pipi(t, Lprime, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa,  Elev_MpiL_4dot2)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

			    double F_int_val = kern_val*func_val;
        

			    return F_int_val;
			  };

      auto F_infL = [&](double t) {

		      //double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
		      //double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
		      double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);

		      double func_val =  Qfact*LL.V_pipi_infL(t, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

		      double F_int_val = kern_val*func_val;
        

		      return F_int_val;

		    };

   

    
    
      //Vfloat Sum_Series_Terms;
      double agm2_summ=0.0;
      double agm2_2L_summ=0.0;
      double agm2_infL_summ=0.0;
      double agm2_Lprime_summ=0.0;
      int Nterms = 2000;
      int Nterms_2L= 2000;
    
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



    //print ensemble infos

    //intermediate and SD window info
    cout<<"%%%%%%% Print info from W and SD windows TM-correlator %%%%%%"<<endl;
    cout<<"SD: "<<agm2_light_SD.ave(i_ens)<<" +- "<<agm2_light_SD.err(i_ens)<<endl;
    cout<<"SD (ELM): "<<agm2_light_SD_ELM.ave(i_ens)<<" +- "<<agm2_light_SD_ELM.err(i_ens)<<endl;
    cout<<"W: "<<agm2_light_W.ave(i_ens)<<" +- "<<agm2_light_W.err(i_ens)<<endl;
    cout<<"W (ELM): "<<agm2_light_W_ELM.ave(i_ens)<<" +- "<<agm2_light_W_ELM.err(i_ens)<<endl;
    cout<<"%%%%%%% End windows info %%%%%%%"<<endl;
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


  //########################################      OSTERWALDER-SEILER CORRELATOR        ########################################################



  cout<<"#########################  OSTERWALDER-SEILER CORRELATOR ###########################"<<endl;



  
 
  for(int i_ens=0; i_ens<Nens_light;i_ens++) {
    cout<<"######### ENSEMBLE: "<<V_light_1.Tag[i_ens]<<endl;
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_light_1.Tag[i_ens]);
  

    distr_t a_distr = a_distr_list[i_ens];
    
    distr_t Ed(UseJack), Rd(UseJack), Mr(UseJack), g(UseJack), kap(UseJack), agm2_dual(UseJack), mp(UseJack), agm2_2L_dual(UseJack), agm2_infL_dual(UseJack), agm2_Lprime_dual(UseJack);
    for(int ijack=0;ijack<Njacks;ijack++) {
      //retrieve fit_parameters
      fit_par my_fit_par = par_list_anal_repr_OS[i_ens][ijack];
      Ed.distr.push_back(my_fit_par.Ed);
      Rd.distr.push_back(my_fit_par.Rd);
      Mr.distr.push_back(my_fit_par.Mrho);
      g.distr.push_back(my_fit_par.gpi);
      kap.distr.push_back(my_fit_par.kappa);

   

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
		  
		    //double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
		    double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		    double func_val =  Qfact*LL.V_pipi(t, L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa, Elev)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

		    double F_int_val = kern_val*func_val;
        

		    return F_int_val;
		  };


      auto F_int_2L = [&](double t) {

			//double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
			//double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
			double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
			double func_val =  Qfact*LL.V_pipi(t, 1.5*L, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa, Elev_2L)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

			double F_int_val = kern_val*func_val;
        

			return F_int_val;
		      };

      auto F_int_Lprime = [&](double t) {

			    //  double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
			    //double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
			    double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
			    double func_val =  Qfact*LL.V_pipi(t, Lprime, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa,  Elev_MpiL_4dot2)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

			    double F_int_val = kern_val*func_val;
        

			    return F_int_val;
			  };

      auto F_infL = [&](double t) {

		      //double kern_val = 4.0*pow(alpha,2)*sqrt(m_rho*a_distr.distr[ijack]/(my_fit_par.Mrho*Mp))*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);
		  
		      //double kern_val = 4.0*pow(alpha,2)*kernel_K(t, a_distr.distr[ijack] );
		      double kern_val= 4.0*pow(alpha,2)*kernel_K(t, my_fit_par.Mrho*Mp/m_rho);

		      double func_val =  Qfact*LL.V_pipi_infL(t, my_fit_par.Mrho*Mp, my_fit_par.gpi*(my_fit_par.Mrho*Mp/fp), Mp, my_fit_par.kappa)+ LL.Vdual(t, my_fit_par.Mrho*Mp, my_fit_par.Ed*Mp, my_fit_par.Rd);

		      double F_int_val = kern_val*func_val;
        

		      return F_int_val;

		    };

   

    
    
      //Vfloat Sum_Series_Terms;
      double agm2_summ=0.0;
      double agm2_2L_summ=0.0;
      double agm2_infL_summ=0.0;
      double agm2_Lprime_summ=0.0;
      int Nterms = 2000;
      int Nterms_2L= 2000;
    
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

  
    

    //push_back
    Edual_OS.distr_list.push_back(Ed);
    Rdual_OS.distr_list.push_back(Rd);
    Mrho_OS.distr_list.push_back(Mr);
    gpi_OS.distr_list.push_back(g);
    Kappa_OS.distr_list.push_back(kap);
    agm2_light_fit_OS.distr_list.push_back(agm2_dual);
    agm2_light_2L_fit_OS.distr_list.push_back(agm2_2L_dual);
    agm2_light_infL_fit_OS.distr_list.push_back(agm2_infL_dual);
    agm2_light_Lprime_fit_OS.distr_list.push_back(agm2_Lprime_dual);



    //print ensemble infos

    //find whether the disconnected diagram corresponding to i_ens has been analyzed
    bool Found_disco_tag=false;
    int disco_ens=0;
    if(Include_light_disco) {
      for(int j=0; j<(signed)disco_Tags.size();j++) { if(disco_Tags[j] == V_light_1.Tag[i_ens]) {Found_disco_tag=true; disco_ens=j; break;}}
    }

    //intermediate and SD window info
    cout<<"%%%%%%% Print info from W and SD windows OS-correlator %%%%%%"<<endl;
    cout<<"SD: "<<agm2_light_SD_OS.ave(i_ens)<<" +- "<<agm2_light_SD_OS.err(i_ens)<<endl;
    cout<<"SD (ELM): "<<agm2_light_SD_ELM_OS.ave(i_ens)<<" +- "<<agm2_light_SD_ELM_OS.err(i_ens)<<endl;
    cout<<"W: "<<agm2_light_W_OS.ave(i_ens)<<" +- "<<agm2_light_W_OS.err(i_ens)<<endl;
    cout<<"W (ELM): "<<agm2_light_W_ELM_OS.ave(i_ens)<<" +- "<<agm2_light_W_ELM_OS.err(i_ens)<<endl;
    if(Found_disco_tag) {
      cout<<"%%%%%% disconnected contribution to the windows  %%%%%%"<<endl;
      cout<<"SD: "<<agm2_disco_light_SD.ave(disco_ens)<<" +- "<<agm2_disco_light_SD.err(disco_ens)<<endl;
      cout<<"SD (ELM): "<<agm2_disco_light_SD_ELM.ave(disco_ens)<<" +- "<<agm2_disco_light_SD_ELM.err(disco_ens)<<endl;
      cout<<"W: "<<agm2_disco_light_W.ave(disco_ens)<<" +- "<<agm2_disco_light_W.err(disco_ens)<<endl;
      cout<<"W (ELM): "<<agm2_disco_light_W_ELM.ave(disco_ens)<<" +- "<<agm2_disco_light_W_ELM.err(disco_ens)<<endl;
    }
    cout<<"%%%%%%% End windows info %%%%%%%"<<endl;
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

  }

  
  //############################################################################################################################################


  
  cout<<"########### DONE #############"<<endl;
  




  //strange
  //tm
  //L
  Print_To_File(V_strange_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  ZV_fit_strange_L.ave(), ZV_fit_strange_L.err(), MV_fit_strange_L.ave(), MV_fit_strange_L.err(), agm2_strange_L.ave(), agm2_strange_L.err(), agm2_strange_No_ELM_L.ave(), agm2_strange_No_ELM_L.err()}, "../data/gm2/strange/tm/agm2_fit_L.list", "", "# Ens L a ml Mpi Mpi_OS MK  ZV   MV    agm2(ELM)  agm2(No ELM)");
  //M
  Print_To_File(V_strange_1_M.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  ZV_fit_strange_M.ave(), ZV_fit_strange_M.err(), MV_fit_strange_M.ave(), MV_fit_strange_M.err(), agm2_strange_M.ave(), agm2_strange_M.err(), agm2_strange_No_ELM_M.ave(), agm2_strange_No_ELM_M.err()}, "../data/gm2/strange/tm/agm2_fit_M.list", "", "# Ens L a ml Mpi Mpi_OS MK  ZV   MV    agm2(ELM), agm2(Non ELM)");
  //H
  /*
  Print_To_File(V_strange_1_H.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(), Mpi_OS_fit_strange.ave(), Mpi_OS_fit_strange.err(), MK_fit_H.ave(), MK_fit_H.err(),  ZV_fit_strange_H.ave(), ZV_fit_strange_H.err(), MV_fit_strange_H.ave(), MV_fit_strange_H.err(), agm2_strange_H.ave(), agm2_strange_H.err()}, "../data/gm2/strange/tm/agm2_fit_H.list", "", "# Ens L a ml Mpi Mpi_OS MK  ZV   MV    agm2");
  */
  //Extr
  Print_To_File(V_strange_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(), agm2_strange_Extr.ave(), agm2_strange_Extr.err(), agm2_strange_No_ELM_Extr.ave(), agm2_strange_No_ELM_Extr.err()}, "../data/gm2/strange/tm/agm2_fit_Extr.list", "", "# Ens L a ml Mpi Mpi_OS agm2(ELM)  agm2(No ELM)");
  
  //OS
  //L
  Print_To_File(V_strange_OS_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(), ZV_fit_strange_OS_L.ave(), ZV_fit_strange_OS_L.err(), MV_fit_strange_OS_L.ave(), MV_fit_strange_OS_L.err(), agm2_strange_OS_L.ave(), agm2_strange_OS_L.err(), agm2_strange_OS_No_ELM_L.ave(), agm2_strange_OS_No_ELM_L.err()}, "../data/gm2/strange/OS/agm2_fit_L.list", "", "# Ens L a ml Mpi Mpi_OS MK  ZV   MV    agm2(ELM)  agm2(No ELM)");
  //M
  Print_To_File(V_strange_OS_1_M.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(), ZV_fit_strange_OS_M.ave(), ZV_fit_strange_OS_M.err(), MV_fit_strange_OS_M.ave(), MV_fit_strange_OS_M.err(), agm2_strange_OS_M.ave(), agm2_strange_OS_M.err(), agm2_strange_OS_No_ELM_M.ave(), agm2_strange_OS_No_ELM_M.err()}, "../data/gm2/strange/OS/agm2_fit_M.list", "", "# Ens L a ml Mpi Mpi_OS MK  ZV   MV    agm2(ELM)  agm2(No ELM)");
  //H
  /*
  Print_To_File(V_strange_OS_1_H.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(), Mpi_OS_fit_strange.ave(), Mpi_OS_fit_strange.err(), MK_fit_H.ave(), MK_fit_H.err(), ZV_fit_strange_OS_H.ave(), ZV_fit_strange_OS_H.err(), MV_fit_strange_OS_H.ave(), MV_fit_strange_OS_H.err(), agm2_strange_OS_H.ave(), agm2_strange_OS_H.err()}, "../data/gm2/strange/OS/agm2_fit_H.list", "", "# Ens L a ml Mpi Mpi_OS MK  ZV   MV    agm2");
  */
  //Extr
  Print_To_File(V_strange_OS_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(), agm2_strange_OS_Extr.ave(), agm2_strange_OS_Extr.err(), agm2_strange_OS_No_ELM_Extr.ave(), agm2_strange_OS_No_ELM_Extr.err()}, "../data/gm2/strange/OS/agm2_fit_Extr.list", "", "# Ens L a ml Mpi Mpi_OS agm2(ELM)  agm2(No ELM)");

     

  //charm
  //tm
  //L
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_L.ave(), MD_fit_L.err(),  ZV_fit_charm_L.ave(), ZV_fit_charm_L.err(), MV_fit_charm_L.ave(), MV_fit_charm_L.err(), agm2_charm_L.ave(), agm2_charm_L.err()}, "../data/gm2/charm/tm/agm2_fit_ELM_L.list", "", "# Ens  L a ml Mpi Mpi_OS MD ZV MV  agm2");
  //M
  Print_To_File(V_charm_1_M.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_M.ave(), MD_fit_M.err(),  ZV_fit_charm_M.ave(), ZV_fit_charm_M.err(), MV_fit_charm_M.ave(), MV_fit_charm_M.err(), agm2_charm_M.ave(), agm2_charm_M.err()}, "../data/gm2/charm/tm/agm2_fit_ELM_M.list", "", "# Ens  L a ml Mpi Mpi_OS MD ZV MV  agm2");
  //H
  Print_To_File(V_charm_1_H.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_H.ave(), MD_fit_H.err(),  ZV_fit_charm_H.ave(), ZV_fit_charm_H.err(), MV_fit_charm_H.ave(), MV_fit_charm_H.err(), agm2_charm_H.ave(), agm2_charm_H.err()}, "../data/gm2/charm/tm/agm2_fit_ELM_H.list", "", "# Ens  L a ml Mpi Mpi_OS MD ZV MV  agm2");
  //Extr
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), agm2_charm_Extr.ave(), agm2_charm_Extr.err()}, "../data/gm2/charm/tm/agm2_fit_ELM_Extr.list", "", "# Ens L a ml Mpi Mpi_OS agm2");

  //OS
  //L
  Print_To_File(V_charm_OS_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_L.ave(), MD_fit_L.err(), ZV_fit_charm_OS_L.ave(), ZV_fit_charm_OS_L.err(), MV_fit_charm_OS_L.ave(), MV_fit_charm_OS_L.err(), agm2_charm_OS_L.ave(), agm2_charm_OS_L.err()}, "../data/gm2/charm/OS/agm2_fit_ELM_L.list", "", "# Ens  L a ml Mpi Mpi_OS MD ZV MV  agm2");
  //M
  Print_To_File(V_charm_OS_1_M.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_M.ave(), MD_fit_M.err(), ZV_fit_charm_OS_M.ave(), ZV_fit_charm_OS_M.err(), MV_fit_charm_OS_M.ave(), MV_fit_charm_OS_M.err(), agm2_charm_OS_M.ave(), agm2_charm_OS_M.err()}, "../data/gm2/charm/OS/agm2_fit_ELM_M.list", "", "# Ens  L a ml Mpi Mpi_OS MD ZV MV  agm2");
  //H
  Print_To_File(V_charm_OS_1_H.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_H.ave(), MD_fit_H.err(), ZV_fit_charm_OS_H.ave(), ZV_fit_charm_OS_H.err(), MV_fit_charm_OS_H.ave(), MV_fit_charm_OS_H.err(), agm2_charm_OS_H.ave(), agm2_charm_OS_H.err()}, "../data/gm2/charm/OS/agm2_fit_ELM_H.list", "", "# Ens  L a ml Mpi Mpi_OS MD ZV MV  agm2");
  //Extr
  Print_To_File(V_charm_OS_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), agm2_charm_OS_Extr.ave(), agm2_charm_OS_Extr.err()}, "../data/gm2/charm/OS/agm2_fit_ELM_Extr.list", "", "# Ens L a ml Mpi Mpi_OS agm2");


  //charm NON_ELM
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_L.ave(), MD_fit_L.err(),  ZV_fit_charm_L.ave(), ZV_fit_charm_L.err(), MV_fit_charm_L.ave(), MV_fit_charm_L.err(), agm2_charm_No_ELM_L.ave(), agm2_charm_No_ELM_L.err()}, "../data/gm2/charm/tm/agm2_fit.list", "", "# Ens  L a ml Mpi Mpi_OS MD ZV MV  agm2");
  //M
  Print_To_File(V_charm_1_M.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_M.ave(), MD_fit_M.err(),  ZV_fit_charm_M.ave(), ZV_fit_charm_M.err(), MV_fit_charm_M.ave(), MV_fit_charm_M.err(), agm2_charm_No_ELM_M.ave(), agm2_charm_No_ELM_M.err()}, "../data/gm2/charm/tm/agm2_fit_M.list", "", "# Ens  L a ml Mpi Mpi_OS MD ZV MV  agm2");
  //H
  Print_To_File(V_charm_1_H.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_H.ave(), MD_fit_H.err(),  ZV_fit_charm_H.ave(), ZV_fit_charm_H.err(), MV_fit_charm_H.ave(), MV_fit_charm_H.err(), agm2_charm_No_ELM_H.ave(), agm2_charm_No_ELM_H.err()}, "../data/gm2/charm/tm/agm2_fit_H.list", "", "# Ens  L a ml Mpi Mpi_OS MD ZV MV  agm2");
  //Extr
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), agm2_charm_No_ELM_Extr.ave(), agm2_charm_No_ELM_Extr.err()}, "../data/gm2/charm/tm/agm2_fit_Extr.list", "", "# Ens L a ml Mpi Mpi_OS agm2");

  //OS
  //L
  Print_To_File(V_charm_OS_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_L.ave(), MD_fit_L.err(), ZV_fit_charm_OS_L.ave(), ZV_fit_charm_OS_L.err(), MV_fit_charm_OS_L.ave(), MV_fit_charm_OS_L.err(), agm2_charm_OS_No_ELM_L.ave(), agm2_charm_OS_No_ELM_L.err()}, "../data/gm2/charm/OS/agm2_fit_L.list", "", "# Ens  L a ml Mpi Mpi_OS MD ZV MV  agm2");
  //M
  Print_To_File(V_charm_OS_1_M.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_M.ave(), MD_fit_M.err(), ZV_fit_charm_OS_M.ave(), ZV_fit_charm_OS_M.err(), MV_fit_charm_OS_M.ave(), MV_fit_charm_OS_M.err(), agm2_charm_OS_No_ELM_M.ave(), agm2_charm_OS_No_ELM_M.err()}, "../data/gm2/charm/OS/agm2_fit_M.list", "", "# Ens  L a ml Mpi Mpi_OS MD ZV MV  agm2");
  //H
  Print_To_File(V_charm_OS_1_H.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_H.ave(), MD_fit_H.err(), ZV_fit_charm_OS_H.ave(), ZV_fit_charm_OS_H.err(), MV_fit_charm_OS_H.ave(), MV_fit_charm_OS_H.err(), agm2_charm_OS_No_ELM_H.ave(), agm2_charm_OS_No_ELM_H.err()}, "../data/gm2/charm/OS/agm2_fit_H.list", "", "# Ens  L a ml Mpi Mpi_OS MD ZV MV  agm2");
  //Extr
  Print_To_File(V_charm_OS_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), agm2_charm_OS_No_ELM_Extr.ave(), agm2_charm_OS_No_ELM_Extr.err()}, "../data/gm2/charm/OS/agm2_fit_Extr.list", "", "# Ens L a ml Mpi Mpi_OS agm2");

  
  //light
  //tm
  Print_To_File(V_light_1.Tag, {ZV_fit_light.ave(), ZV_fit_light.err(), MV_fit_light.ave(), MV_fit_light.err()}, "../data/gm2/light/tm/ZV_MV_fitted.list", "", "#Ens ZV MV");
  //OS
  Print_To_File(V_light_1.Tag, {ZV_fit_light_OS.ave(), ZV_fit_light_OS.err(), MV_fit_light_OS.ave(), MV_fit_light_OS.err()}, "../data/gm2/light/OS/ZV_MV_fitted.list", "", "#Ens ZV MV");



  //print informations on the windows
  //light
  //tm
  Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), agm2_light_W.ave(), agm2_light_W.err(), agm2_light_W_ELM.ave(), agm2_light_W_ELM.err(), agm2_light_SD.ave(), agm2_light_SD.err(), agm2_light_SD_ELM.ave(), agm2_light_SD_ELM.err()}, "../data/gm2/light/tm/windows.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs   W    W(ELM)     SD     SD(ELM) ");
  //OS
  Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), agm2_light_W_OS.ave(), agm2_light_W_OS.err(), agm2_light_W_ELM_OS.ave(), agm2_light_W_ELM_OS.err(), agm2_light_SD_OS.ave(), agm2_light_SD_OS.err(), agm2_light_SD_ELM_OS.ave(), agm2_light_SD_ELM_OS.err()}, "../data/gm2/light/OS/windows.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs   W    W(ELM)     SD     SD(ELM) ");
  //disco
  if(Include_light_disco) {
  Print_To_File(disco_Tags,{L_list_disco, a_list_disco, ml_list_disco, Mpi_fit_disco.ave(), Mpi_fit_disco.err(), Mpi_OS_fit_disco.ave(), Mpi_OS_fit_disco.err(), fp_fit_disco.ave(), fp_fit_disco.err(), Zv_fit_disco.ave(), Zv_fit_disco.err(), agm2_disco_light_W.ave(), agm2_disco_light_W.err(), agm2_disco_light_W_ELM.ave(), agm2_disco_light_W_ELM.err(), agm2_disco_light_SD.ave(), agm2_disco_light_SD.err(), agm2_disco_light_SD_ELM.ave(), agm2_disco_light_SD_ELM.err()} , "../data/gm2/light/disco/windows.list", "", "#ENS L a ml Mpi_tm Mpi_OS fp Zv W W(ELM) SD SD(ELM)");
  }

  
  //strange
  //tm
  //L
  Print_To_File(V_strange_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  MV_fit_strange_L.ave(), MV_fit_strange_L.err(),  agm2_strange_W_L.ave(), agm2_strange_W_L.err(), agm2_strange_W_ELM_L.ave(), agm2_strange_W_ELM_L.err(), agm2_strange_SD_L.ave(), agm2_strange_SD_L.err(), agm2_strange_SD_ELM_L.ave(), agm2_strange_SD_ELM_L.err()}, "../data/gm2/strange/tm/windows_L.list", "", "#ENS L a ml Mpi_tm Mpi_OS MK MV W   W(ELM)   SD    SD(ELM)");
  //M
  Print_To_File(V_strange_1_M.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  MV_fit_strange_M.ave(), MV_fit_strange_M.err(),  agm2_strange_W_M.ave(), agm2_strange_W_M.err(), agm2_strange_W_ELM_M.ave(), agm2_strange_W_ELM_M.err(), agm2_strange_SD_M.ave(), agm2_strange_SD_M.err(), agm2_strange_SD_ELM_M.ave(), agm2_strange_SD_ELM_M.err()}, "../data/gm2/strange/tm/windows_M.list", "", "#ENS L a ml Mpi_tm Mpi_OS MK MV W   W(ELM)   SD    SD(ELM)");
  //H
  /*
  Print_To_File(V_strange_1_H.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(), Mpi_OS_fit_strange.ave(), Mpi_OS_fit_strange.err(), MK_fit_H.ave(), MK_fit_H.err(),  MV_fit_strange_H.ave(), MV_fit_strange_H.err(),  agm2_strange_W_H.ave(), agm2_strange_W_H.err(), agm2_strange_W_ELM_H.ave(), agm2_strange_W_ELM_H.err(), agm2_strange_SD_H.ave(), agm2_strange_SD_H.err(), agm2_strange_SD_ELM_H.ave(), agm2_strange_SD_ELM_H.err()}, "../data/gm2/strange/tm/windows_H.list", "", "#ENS L a ml Mpi_tm Mpi_OS MK MV W   W(ELM)   SD    SD(ELM)");
  */
  //Extr
  Print_To_File(V_strange_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  agm2_strange_W_Extr.ave(), agm2_strange_W_Extr.err(), agm2_strange_W_ELM_Extr.ave(), agm2_strange_W_ELM_Extr.err(), agm2_strange_SD_Extr.ave(), agm2_strange_SD_Extr.err(), agm2_strange_SD_ELM_Extr.ave(), agm2_strange_SD_ELM_Extr.err()}, "../data/gm2/strange/tm/windows_Extr.list", "", "#ENS L a ml Mpi_tm Mpi_OS  W   W(ELM)   SD    SD(ELM)");

  
  //OS
  //L
  Print_To_File(V_strange_OS_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  MV_fit_strange_OS_L.ave(), MV_fit_strange_OS_L.err(), agm2_strange_W_OS_L.ave(), agm2_strange_W_OS_L.err(), agm2_strange_W_ELM_OS_L.ave(), agm2_strange_W_ELM_OS_L.err(), agm2_strange_SD_OS_L.ave(), agm2_strange_SD_OS_L.err(), agm2_strange_SD_ELM_OS_L.ave(), agm2_strange_SD_ELM_OS_L.err()}, "../data/gm2/strange/OS/windows_L.list", "", "#ENS L a ml Mpi_tm Mpi_OS MK MV W   W(ELM)   SD    SD(ELM)");
  //M
  Print_To_File(V_strange_OS_1_M.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(), MV_fit_strange_OS_M.ave(), MV_fit_strange_OS_M.err(), agm2_strange_W_OS_M.ave(), agm2_strange_W_OS_M.err(), agm2_strange_W_ELM_OS_M.ave(), agm2_strange_W_ELM_OS_M.err(), agm2_strange_SD_OS_M.ave(), agm2_strange_SD_OS_M.err(), agm2_strange_SD_ELM_OS_M.ave(), agm2_strange_SD_ELM_OS_M.err()}, "../data/gm2/strange/OS/windows_M.list", "", "#ENS L a ml Mpi_tm Mpi_OS MK MV W   W(ELM)   SD    SD(ELM)");
  //H
  /*
  Print_To_File(V_strange_OS_1_H.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit_strange.ave(), Mpi_fit_strange.err(), Mpi_OS_fit_strange.ave(), Mpi_OS_fit_strange.err(), MK_fit_H.ave(), MK_fit_H.err(),  MV_fit_strange_OS_H.ave(), MV_fit_strange_OS_H.err(), agm2_strange_W_OS_H.ave(), agm2_strange_W_OS_H.err(), agm2_strange_W_ELM_OS_H.ave(), agm2_strange_W_ELM_OS_H.err(), agm2_strange_SD_OS_H.ave(), agm2_strange_SD_OS_H.err(), agm2_strange_SD_ELM_OS_H.ave(), agm2_strange_SD_ELM_OS_H.err()}, "../data/gm2/strange/OS/windows_H.list", "", "#ENS L a ml Mpi_tm Mpi_OS MK MV W   W(ELM)   SD    SD(ELM)");
  */
  //Extr
  Print_To_File(V_strange_OS_1_L.Tag, {L_strange_list, a_strange_list, ml_strange_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  agm2_strange_W_OS_Extr.ave(), agm2_strange_W_OS_Extr.err(), agm2_strange_W_ELM_OS_Extr.ave(), agm2_strange_W_ELM_OS_Extr.err(), agm2_strange_SD_OS_Extr.ave(), agm2_strange_SD_OS_Extr.err(), agm2_strange_SD_ELM_OS_Extr.ave(), agm2_strange_SD_ELM_OS_Extr.err()}, "../data/gm2/strange/OS/windows_Extr.list", "", "#ENS L a ml Mpi_tm Mpi_OS  W   W(ELM)   SD    SD(ELM)");



  
  //charm
  //tm
  //L
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_L.ave(), MD_fit_L.err(),  MV_fit_charm_L.ave(), MV_fit_charm_L.err(),  agm2_charm_W_L.ave(), agm2_charm_W_L.err(), agm2_charm_W_ELM_L.ave(), agm2_charm_W_ELM_L.err(), agm2_charm_SD_L.ave(), agm2_charm_SD_L.err(), agm2_charm_SD_ELM_L.ave(), agm2_charm_SD_ELM_L.err()}, "../data/gm2/charm/tm/windows_L.list", "", "#ENS L a ml Mpi_tm Mpi_OS MD MV W   W(ELM)   SD    SD(ELM)");
  //M
  Print_To_File(V_charm_1_M.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_M.ave(), MD_fit_M.err(),  MV_fit_charm_M.ave(), MV_fit_charm_M.err(),  agm2_charm_W_M.ave(), agm2_charm_W_M.err(), agm2_charm_W_ELM_M.ave(), agm2_charm_W_ELM_M.err(), agm2_charm_SD_M.ave(), agm2_charm_SD_M.err(), agm2_charm_SD_ELM_M.ave(), agm2_charm_SD_ELM_M.err()}, "../data/gm2/charm/tm/windows_M.list", "", "#ENS L a ml Mpi_tm Mpi_OS MD MV W   W(ELM)   SD    SD(ELM)");
  //H
  Print_To_File(V_charm_1_H.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_H.ave(), MD_fit_H.err(),  MV_fit_charm_H.ave(), MV_fit_charm_H.err(),  agm2_charm_W_H.ave(), agm2_charm_W_H.err(), agm2_charm_W_ELM_H.ave(), agm2_charm_W_ELM_H.err(), agm2_charm_SD_H.ave(), agm2_charm_SD_H.err(), agm2_charm_SD_ELM_H.ave(), agm2_charm_SD_ELM_H.err()}, "../data/gm2/charm/tm/windows_H.list", "", "#ENS L a ml Mpi_tm Mpi_OS MD MV W   W(ELM)   SD    SD(ELM)");
  //Extr
  Print_To_File(V_charm_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  agm2_charm_W_Extr.ave(), agm2_charm_W_Extr.err(), agm2_charm_W_ELM_Extr.ave(), agm2_charm_W_ELM_Extr.err(), agm2_charm_SD_Extr.ave(), agm2_charm_SD_Extr.err(), agm2_charm_SD_ELM_Extr.ave(), agm2_charm_SD_ELM_Extr.err()}, "../data/gm2/charm/tm/windows_Extr.list", "", "#ENS L a ml Mpi_tm Mpi_OS  W   W(ELM)   SD    SD(ELM)");
  
  //OS
  //L
  Print_To_File(V_charm_OS_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_L.ave(), MD_fit_L.err(),  MV_fit_charm_OS_L.ave(), MV_fit_charm_OS_L.err(),  agm2_charm_W_OS_L.ave(), agm2_charm_W_OS_L.err(), agm2_charm_W_ELM_OS_L.ave(), agm2_charm_W_ELM_OS_L.err(), agm2_charm_SD_OS_L.ave(), agm2_charm_SD_OS_L.err(), agm2_charm_SD_ELM_OS_L.ave(), agm2_charm_SD_ELM_OS_L.err()}, "../data/gm2/charm/OS/windows_L.list", "", "#ENS L a ml Mpi_tm Mpi_OS MD  MV   W   W(ELM)   SD    SD(ELM)");
  //M
  Print_To_File(V_charm_OS_1_M.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_M.ave(), MD_fit_M.err(),  MV_fit_charm_OS_M.ave(), MV_fit_charm_OS_M.err(),  agm2_charm_W_OS_M.ave(), agm2_charm_W_OS_M.err(), agm2_charm_W_ELM_OS_M.ave(), agm2_charm_W_ELM_OS_M.err(), agm2_charm_SD_OS_M.ave(), agm2_charm_SD_OS_M.err(), agm2_charm_SD_ELM_OS_M.ave(), agm2_charm_SD_ELM_OS_M.err()}, "../data/gm2/charm/OS/windows_M.list", "", "#ENS L a ml Mpi_tm Mpi_OS MD  MV   W   W(ELM)   SD    SD(ELM)");
  //H
  Print_To_File(V_charm_OS_1_H.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(), MD_fit_H.ave(), MD_fit_H.err(),  MV_fit_charm_OS_H.ave(), MV_fit_charm_OS_H.err(),  agm2_charm_W_OS_H.ave(), agm2_charm_W_OS_H.err(), agm2_charm_W_ELM_OS_H.ave(), agm2_charm_W_ELM_OS_H.err(), agm2_charm_SD_OS_H.ave(), agm2_charm_SD_OS_H.err(), agm2_charm_SD_ELM_OS_H.ave(), agm2_charm_SD_ELM_OS_H.err()}, "../data/gm2/charm/OS/windows_H.list", "", "#ENS L a ml Mpi_tm Mpi_OS MD  MV   W   W(ELM)   SD    SD(ELM)");
  //Extr
  Print_To_File(V_charm_OS_1_L.Tag, {L_charm_list, a_charm_list, ml_charm_list, Mpi_fit_charm.ave(), Mpi_fit_charm.err(), Mpi_OS_fit_charm.ave(), Mpi_OS_fit_charm.err(),  agm2_charm_W_OS_Extr.ave(), agm2_charm_W_OS_Extr.err(), agm2_charm_W_ELM_OS_Extr.ave(), agm2_charm_W_ELM_OS_Extr.err(), agm2_charm_SD_OS_Extr.ave(), agm2_charm_SD_OS_Extr.err(), agm2_charm_SD_ELM_OS_Extr.ave(), agm2_charm_SD_ELM_OS_Extr.err()}, "../data/gm2/charm/OS/windows_Extr.list", "", "#ENS L a ml Mpi_tm Mpi_OS  W   W(ELM)   SD    SD(ELM)");
  
  

  if(!Skip_total_light_calc) {
 
  //print fitted pars for all ensembles
  //tm
  Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), Edual.ave(), Edual.err(), Rdual.ave(), Rdual.err(), Mrho.ave(), Mrho.err(), (Mrho*Mpi_fit/a_distr_list).ave(), (Mrho*Mpi_fit/a_distr_list).err(), gpi.ave(), gpi.err()}, "../data/gm2/light/tm/fit_pars.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs  Edual Rdual Mrho g");
  //OS
  Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), Edual_OS.ave(), Edual_OS.err(), Rdual_OS.ave(), Rdual_OS.err(), Mrho_OS.ave(), Mrho_OS.err(), (Mrho_OS*Mpi_fit/a_distr_list).ave(), (Mrho_OS*Mpi_fit/a_distr_list).err(), gpi_OS.ave(), gpi_OS.err()}, "../data/gm2/light/OS/fit_pars.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs  Edual Rdual Mrho g");

  
  
  //tm and OS pion mass, decay constant, RCs and agm2
  //tm
  Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), agm2_light_fit.ave(), agm2_light_fit.err(), agm2_light_2L_fit.ave(), agm2_light_2L_fit.err(), agm2_light_Lprime_fit.ave(), agm2_light_Lprime_fit.err()   , agm2_light_infL_fit.ave(), agm2_light_infL_fit.err()}, "../data/gm2/light/tm/agm2_fit.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs    agm2(L)    agm2(1.5L) agm2(Mpi*L=4.2)   agm2(infL)");
  //OS
    Print_To_File(V_light_1.Tag, {L_list, a_list, ml_list, Mpi_fit.ave(), Mpi_fit.err(), Mpi_OS_fit.ave(), Mpi_OS_fit.err(),  fp_fit.ave(), fp_fit.err(), Zv_fit.ave(), Zv_fit.err(), Za_fit.ave(), Za_fit.err(), Zp_ov_Zs_fit.ave(), Zp_ov_Zs_fit.err(), agm2_light_fit_OS.ave(), agm2_light_fit_OS.err(), agm2_light_2L_fit_OS.ave(), agm2_light_2L_fit_OS.err(), agm2_light_Lprime_fit_OS.ave(), agm2_light_Lprime_fit_OS.err()   , agm2_light_infL_fit_OS.ave(), agm2_light_infL_fit_OS.err()}, "../data/gm2/light/OS/agm2_fit.list", "", "#ENS L a ml  Mpi_tm  Mpi_OS fp  Zv   Za   Zp/Zs    agm2(L)    agm2(1.5L) agm2(Mpi*L=4.2)   agm2(infL)");


  }


  //print Zv and Za (Hadronic and RI-MOM) from strange run
  for(int is=0;is<Nens_strange;is++) {
    cout<<"####################"<<endl;
    cout<<"Ens: "<<V_strange_1_L.Tag[is]<<endl;
    cout<<"Za hadr-RIMOM (s (l), s (h),  u) "<<Za_fit_strange.ave(is)<<" +- "<<Za_fit_strange.err(is)<<" , "<<Za_fit_strange_heavy.ave(is)<<" +- "<<Za_fit_strange_heavy.err(is)<<" , "<<Za_WI.ave(is)<<" +- "<<Za_WI.err(is)<<endl;
    cout<<"Zv hadr-RIMOM (s (l), s (h),  u) "<<Zv_fit_strange.ave(is)<<" +- "<<Zv_fit_strange.err(is)<<" , "<<Zv_fit_strange_heavy.ave(is)<<" +- "<<Zv_fit_strange_heavy.err(is)<<" , "<<Zv_WI.ave(is)<<" +- "<<Zv_WI.err(is)<<endl;
    cout<<"####################"<<endl;
  }
  //print Zv and Za (Hadronic) from charm run
    for(int is=0;is<Nens_charm;is++) {
    cout<<"####################"<<endl;
    cout<<"Ens: "<<V_charm_1_L.Tag[is]<<endl;
    cout<<"Za hadr (c) (L) "<<Za_fit_charm_L.ave(is)<<" +- "<<Za_fit_charm_L.err(is)<<endl;
    cout<<"Za hadr (c) (M) "<<Za_fit_charm_M.ave(is)<<" +- "<<Za_fit_charm_M.err(is)<<endl;
    cout<<"Za hadr (c) (H) "<<Za_fit_charm_H.ave(is)<<" +- "<<Za_fit_charm_H.err(is)<<endl;
    cout<<"Zv hadr (c) (L) "<<Zv_fit_charm_L.ave(is)<<" +- "<<Zv_fit_charm_L.err(is)<<endl;
    cout<<"Zv hadr (c) (M) "<<Zv_fit_charm_M.ave(is)<<" +- "<<Zv_fit_charm_M.err(is)<<endl;
    cout<<"Zv hadr (c) (H) "<<Zv_fit_charm_H.ave(is)<<" +- "<<Zv_fit_charm_H.err(is)<<endl;
    cout<<"####################"<<endl;
    }



  

    //#############################################        CONTINUUM/THERMODYNAMIC/PHYSICAL-POINT EXTRAPOLATION       ###############################################
    vector<string> a2_list({"on", "tm", "OS", "off"});
    vector<string> FSEs_list({"tm", "off", "on"});
    vector<string> a4_list({"off", "tm"});
    vector<string> mass_extr_list({"on", "off"});
    vector<string> single_fit_list({"off", "OS"});
    VPfloat n_m_pair_list({make_pair(0,0), make_pair(1,1), make_pair(2,2), make_pair(3,3)});

    
    Perform_Akaike_fits(agm2_light_W_ELM, agm2_light_W_ELM_OS, a_A, a_B, a_C, a_D, L_list, a_distr_list, Mpi_fit,fp_fit, V_light_1.Tag, UseJack, Njacks, Nboots, "W_win_ELM", "light",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, 200.0);

    
    FSEs_list = {"off"};
    a4_list={"off", "tm", "OS"};
    mass_extr_list={"on", "off"};
    single_fit_list = {"off", "OS", "tm"};


    //W strange
    Perform_Akaike_fits(agm2_strange_W_ELM_Extr, agm2_strange_W_ELM_OS_Extr, a_A, a_B, a_C, a_D, L_strange_list, a_distr_list_strange, Mpi_fit, fp_fit, V_strange_1_L.Tag, UseJack, Njacks, Nboots, "W_win_ELM", "strange",a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, 27.0);


    
    //SD strange
    Perform_Akaike_fits(agm2_strange_SD_Extr, agm2_strange_SD_OS_Extr, a_A, a_B, a_C, a_D, L_strange_list, a_distr_list_strange, Mpi_fit, fp_fit, V_strange_1_L.Tag, UseJack, Njacks, Nboots, "SD_win", "strange", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, 9.0);


    FSEs_list = {"off"};


    //total strange
    Perform_Akaike_fits(agm2_strange_Extr, agm2_strange_OS_Extr, a_A, a_B, a_C, a_D, L_strange_list, a_distr_list_strange, Mpi_fit, fp_fit, V_strange_1_L.Tag, UseJack, Njacks, Nboots, "total", "strange", a2_list, FSEs_list, a4_list, mass_extr_list, single_fit_list, n_m_pair_list, 50.0);




    //################################################################################################################################################################

   


    exit(-1);
   

    



  

  //print kernel
  Plot_kernel_K(2000);
  

  
   



  return; 

}
