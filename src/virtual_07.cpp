#include "../include/virtual_07.h"
#include "Corr_analysis.h"
#include "numerics.h"
using namespace std;

bool verbose_lev_07=1;
Vfloat sigmas_07({1.5,1.25,1.0, 0.75, 0.5,0.4,0.3,0.15}); //sigma in GeV  
int prec_07=128;
const string MODE_FF="TANT";
const bool Skip_spectral_reconstruction_07=false;
const double Mjpsi= 3.0969; //GeV
const double Mphi= 1.019461; //GeV
const double MDs_phys = 1.96847; // GeV
const double E0_fact = 0.90;
const bool CONS_EM_CURRENT = true;
const string SM_TYPE = "FF_Exp";
const double QU = -1.0/3;
const double QD = -1.0/3;



rt_07_Bs Get_virtual_tensor_FF(int n_xg, bool UseJack, int Njacks, string MESON,  string Corr_path,string path_out) {


  rt_07_Bs return_class;

  PrecFloat::setDefaultPrecision(prec_07);
  cout<<"max possible exponent: "<<PrecFloat::getEmax_max()<<endl;
  cout<<"current max exponent: "<<PrecFloat::getEmax()<<endl;
  cout<<"Number of xg to analyze: "<<n_xg<<endl;

  int t_07_s=14;
  int t_07_s_HLT=14;
  int t_07_c=25;

  string TAG_CURR="";
  if(CONS_EM_CURRENT==false) TAG_CURR="LOC_";

 
  int size_mu_nu= 4;

  double sign_kz=-1.0; //correct one is -1

  //BK
  vector<vector<vector<data_t>>> C_B_u_data(size_mu_nu), C_B_d_data(size_mu_nu);
  //T
  vector<vector<vector<data_t>>> C_T_u_data(size_mu_nu), C_T_d_data(size_mu_nu);


  data_t data_2pts_SM, data_2pts_SMSM;

 
  


  
  for(int mu=0;mu<size_mu_nu;mu++) {

    C_B_u_data[mu].resize(size_mu_nu);
    C_B_d_data[mu].resize(size_mu_nu);
    C_T_u_data[mu].resize(size_mu_nu);
    C_T_d_data[mu].resize(size_mu_nu);

   
    for(int nu=0;nu<size_mu_nu;nu++) {

      C_B_u_data[mu][nu].resize(n_xg);
      C_B_d_data[mu][nu].resize(n_xg);
      C_T_u_data[mu][nu].resize(n_xg);
      C_T_d_data[mu][nu].resize(n_xg);

    
    }
  }

  //custom sorting of gauge confs
  auto Sort_confs = [](string A, string B) {

			   

    int conf_length_A= A.length();
    int conf_length_B= B.length();

    int pos_a_slash=-1;
    int pos_b_slash=-1;
    for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
    for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;

    string A_bis= A.substr(pos_a_slash+1);
    string B_bis= B.substr(pos_b_slash+1);

					     
    string conf_num_A = A_bis.substr(0,4);
    string conf_num_B = B_bis.substr(0,4);
							       
		      
    string rA = A_bis.substr(A_bis.length()-2);
    string rB = B_bis.substr(B_bis.length()-2);
    if(rA.substr(0,1) == "r") { 
      int n1 = stoi(A_bis.substr(A_bis.length()-1));
      int n2 = stoi(B_bis.substr(B_bis.length()-1));
      if(rA == rB) {
	if(rA=="r0" || rA=="r2") return conf_num_A > conf_num_B;
	else if(rA=="r1" || rA=="r3") return conf_num_A < conf_num_B;
	else crash("stream not recognized");
      }
      else return n1<n2;
    }
    return A_bis<B_bis;
  };
  
  //read data

  data_2pts_SM.Read(Corr_path+"_mass", "mes_contr_2pts_SM_3", "P5P5", Sort_confs);
  data_2pts_SMSM.Read(Corr_path+"_mass", "mes_contr_2pts_SMSM_3", "P5P5", Sort_confs);

  
 
  //loop over mu and nu axial
  vector<pair<int,int>> mu_nu_pair_B({make_pair(1,1),make_pair(2,2)});
  vector<pair<int,int>> mu_nu_pair_T({make_pair(1,2),make_pair(2,1)});
   
  
   
  for(int ixg=0;ixg<n_xg;ixg++) {


   
    //B
    for(auto &pair_B : mu_nu_pair_B) {
      int mu=pair_B.first;
      int nu=pair_B.second;
  
      string Tag_contr="S0P5";
      if(CONS_EM_CURRENT==false) Tag_contr="V"+to_string(mu)+"P5";
      //B
      //u
      C_B_u_data[mu][nu][ixg].Read(Corr_path, TAG_CURR+"C_u_B_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      //d
      C_B_d_data[mu][nu][ixg].Read(Corr_path, TAG_CURR+"C_d_B_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      
    
    }

    //T
    for(auto &pair_T : mu_nu_pair_T) {
      int mu=pair_T.first;
      int nu=pair_T.second;
    
      string Tag_contr="S0P5";
      if(CONS_EM_CURRENT==false) Tag_contr="V"+to_string(mu)+"P5";
      //T
      //u
      C_T_u_data[mu][nu][ixg].Read(Corr_path, TAG_CURR+"C_u_T_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      //d
      C_T_d_data[mu][nu][ixg].Read(Corr_path, TAG_CURR+"C_d_T_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
      
    
    }
  }



   
  GaussianMersenne GM(652205123);

  //resample RCs
  distr_t ZT_A(UseJack), ZT_B(UseJack), ZT_C(UseJack), ZT_D(UseJack);
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);

  
  double fmTGeV= 1.0/0.197327;

  LatticeInfo L_info_A, L_info_B, L_info_C, L_info_D;
  L_info_A.LatInfo_new_ens("cA211a.12.48");
  L_info_B.LatInfo_new_ens("cB211b.072.96");
  L_info_C.LatInfo_new_ens("cC211a.06.80");
  L_info_D.LatInfo_new_ens("cD211a.054.96");
  

  for(int ijack=0; ijack<Njacks;ijack++) {

    a_A.distr.push_back( L_info_A.a_from_afp*fmTGeV + GM()*L_info_A.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_B.distr.push_back( L_info_B.a_from_afp*fmTGeV + GM()*L_info_B.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_C.distr.push_back( L_info_C.a_from_afp*fmTGeV + GM()*L_info_C.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_D.distr.push_back( L_info_D.a_from_afp*fmTGeV + GM()*L_info_D.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));

        
    ZT_A.distr.push_back( L_info_A.ZT_RI2 + GM()*L_info_A.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_B.distr.push_back( L_info_B.ZT_RI2 + GM()*L_info_B.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_C.distr.push_back( L_info_C.ZT_RI2 + GM()*L_info_C.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    ZT_D.distr.push_back( L_info_D.ZT_RI2 + GM()*L_info_D.ZT_RI2_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    

  }

  int Nens= C_B_u_data[1][1][0].size;
  vector<string> Ens_tags= C_B_u_data[1][1][0].Tag;
  Vint Nts=C_B_u_data[1][1][0].nrows;

  distr_t_list xg_list(UseJack);
  vector<distr_t_list> F_T_u_list;
  vector<distr_t_list> F_T_d_I_list;
  vector<distr_t_list> FV_T_u_real_list;
  vector<distr_t_list> FV_T_d_real_list;
  vector<distr_t_list> FA_T_u_real_list;
  vector<distr_t_list> FA_T_d_real_list;
  vector<vector<distr_t_list>> F_T_d_RE_sm_list(n_xg);
  vector<vector<distr_t_list>> F_T_d_IM_sm_list(n_xg);

  for(int ixg=0; ixg<n_xg;ixg++) {

    F_T_u_list.emplace_back(UseJack);
    F_T_d_I_list.emplace_back(UseJack);
    FV_T_u_real_list.emplace_back(UseJack);
    FV_T_d_real_list.emplace_back(UseJack);
    FA_T_u_real_list.emplace_back(UseJack);
    FA_T_d_real_list.emplace_back(UseJack);
    for(int is=0; is<(signed)sigmas_07.size(); is++) {
      F_T_d_RE_sm_list[ixg].emplace_back(UseJack);
      F_T_d_IM_sm_list[ixg].emplace_back(UseJack);
    }
  }


  auto K_RE= [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {


    if(SM_TYPE=="FF_Gauss") {
      PrecFloat x = (E-m)/(sqrt(PrecFloat(2))*s);
      return sqrt(2)*DawsonF(x)/s;
    }


    if(SM_TYPE=="FF_Gauss_Schwartz") {
      
      PrecFloat x = (E-m);
      if( abs(x)/s > 0.1) return ( 1 - exp( -x*x/(2*s*s)))/x;
      else {
        int n=2;
        bool converged=false;
        PrecFloat fact = -pow(x/(sqrt(PrecFloat(2))*s),2);
        PrecFloat sum= x/(2*s*s);
        PrecFloat M= sum;
        PrecFloat prec_sum;
        while(!converged) {
          M *= fact/n;
          prec_sum=sum;
          sum += M;
          n++;
          if(prec_sum==sum) converged=true;
        }
        return sum;
      }
      
    }
    
    
    if(SM_TYPE=="FF_Cauchy") {
      PrecFloat t = (E-m);
      return t/( t*t + s*s);
    }

    if(SM_TYPE=="FF_Exp") {
      
      PrecFloat x= (E-m);
      
      
      return ( cos(s)*exp(-x)-exp(-2*x))/( 1 + exp(-2*x) -2*cos(s)*exp(-x));
    }
    
    if(SM_TYPE=="FF_Sinh_half") {

      PrecFloat x = (E-m);

      return exp(-x/2)*2*sinh(x/2)/( pow(2*sinh(x/2),2) + s*s);

    }

    PrecFloat norm;
    if( s > 1) norm= PrecFloat(2)*log( s + sqrt( s*s -1))/sqrt(s*s -1);
    else if(s==1) norm=PrecFloat(2);
    else {
      PrecFloat phi= abs(atan( sqrt( 1 - s*s)/s));
      norm= PrecFloat(2)*phi/sqrt( 1 - s*s);
    }
    norm /= precPi();
    norm = 1/norm;
    
       
    PrecFloat t = (E-m);
    PrecFloat x=sinh(t);
    PrecFloat res= (x + (s*s/x));
    res=1/res;

    if( abs(t) >= 1) return norm*res;
    else return norm*x/( s*s + x*x);
      
    exit(-1);
    return 0;
  };

  //smeared kernel of the immaginary part
  auto K_IM = [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

    if(SM_TYPE=="FF_Gauss") return precPi()*Get_exact_gauss(E, m, s, E0);
    if(SM_TYPE=="FF_Gauss_Schwartz") return precPi()*Get_exact_gauss(E, m, s, E0);
    if(SM_TYPE=="FF_Cauchy") {
      PrecFloat t= (E-m);
      return s/( t*t + s*s);
    }

    if(SM_TYPE=="FF_Exp") {

      PrecFloat x= (E-m);

      //PrecFloat phi= abs(atan( sin(s)/(1-cos(s))));
      
      //PrecFloat norm= PrecFloat(2)*phi/precPi() ;

      return (exp(-x)*sin(s))/( 1 + exp(-2*x) -2*cos(s)*exp(-x));
      
    }
    
    if(SM_TYPE=="FF_Sinh_half") {
      
      PrecFloat x = (E-m);
      
      return exp(-x/2)*s/( pow(2*sinh(x/2),2) + s*s);

    }
    
    
    
    PrecFloat norm;
    if( s > 1) norm= PrecFloat(2)*log( s + sqrt( s*s -1))/sqrt(s*s -1);
    else if(s==1) norm=PrecFloat(2);
    else {
      PrecFloat phi= abs(atan( sqrt( 1 - s*s)/s));
      norm= PrecFloat(2)*phi/sqrt( 1 - s*s);
    }
    
    norm /= precPi();
    norm = 1/norm;

    PrecFloat t = E-m;
    PrecFloat x = sinh(t);
    return norm*s/( pow(x,2) + pow(s,2));
    
  };
  
  
  boost::filesystem::create_directory( path_out+"/corr_2pts");
  boost::filesystem::create_directory( path_out+"/corr_3pts");
  boost::filesystem::create_directory( path_out+"/mass");
  boost::filesystem::create_directory( path_out+"/covariance");
  boost::filesystem::create_directory( path_out+"/FF_d_I");
  boost::filesystem::create_directory( path_out+"/FF_d_II");
  boost::filesystem::create_directory( path_out+"/FF_u");
  boost::filesystem::create_directory( path_out+"/FF_d");
  boost::filesystem::create_directory( path_out+"/FF");
  
  
  //loop over ensembles
  for(int iens=0; iens<Nens;iens++) {


    boost::filesystem::create_directory( path_out+"/"+data_2pts_SM.Tag[iens]);
    boost::filesystem::create_directory( path_out+"/corr");
    
    CorrAnalysis Corr(UseJack,Njacks,1000);
    Corr.Nt = Nts[iens];
    CorrAnalysis Corr_boot(false,Njacks, 1000);
    Corr_boot.Nt= Nts[iens];
   

    LatticeInfo L_info;
    L_info.LatInfo_new_ens(data_2pts_SM.Tag[iens]);
   

    //read theta values and loop over them
    Vfloat thetas, masses_u, masses_d;

    thetas= Read_From_File(Corr_path+"/"+Ens_tags[iens]+"/pars_list.dat", 1 , 5);
    masses_u= Read_From_File(Corr_path+"/"+Ens_tags[iens]+"/pars_list.dat", 3 , 5);
    masses_d= Read_From_File( Corr_path+"/"+Ens_tags[iens]+"/pars_list.dat", 4 , 5);

    double mu= masses_u[0];
    double md= masses_d[0];


    //RCs
    distr_t ZT, a_distr;
    if(data_2pts_SM.Tag[iens].substr(1,1)=="A") { ZT=ZT_A;a_distr=a_A;}
    else if(data_2pts_SM.Tag[iens].substr(1,1)=="B") { ZT= ZT_B;a_distr=a_B;}
    else if(data_2pts_SM.Tag[iens].substr(1,1)=="C") { ZT= ZT_C;a_distr=a_C;}
    else if(data_2pts_SM.Tag[iens].substr(1,1)=="D") { ZT= ZT_D;a_distr=a_D;}
    else crash("Ensemble: "+data_2pts_SM.Tag[iens]+" not recognised");


    //2pts plateaux
    if(data_2pts_SM.Tag[iens].substr(1,1) =="A") {
      if(MESON == "B1s") { Corr.Tmin=14; Corr.Tmax= 24;}
      else if(MESON == "B2s" ) { Corr.Tmin=13; Corr.Tmax= 24;}
      else if(MESON=="B3s") { Corr.Tmin=14; Corr.Tmax=20;}
      else if(MESON=="B4s") { Corr.Tmin=14; Corr.Tmax=19;}
      else { Corr.Tmin=14; Corr.Tmax=27;}
      
    }
    else if(data_2pts_SM.Tag[iens] =="cB211b.072.64") {
      if(MESON == "B1s") {Corr.Tmin=24; Corr.Tmax=38;}
      else if(MESON == "B2s") {Corr.Tmin=22; Corr.Tmax=34;}
      else if(MESON=="B3s") {Corr.Tmin=17; Corr.Tmax=23;}
      else if(MESON=="B4s") { Corr.Tmin=16; Corr.Tmax=22;}
      else {Corr.Tmin=23; Corr.Tmax=33;}
    }
    else if(data_2pts_SM.Tag[iens] =="cB211b.072.96") {Corr.Tmin=20; Corr.Tmax=36;}
    
    else if(data_2pts_SM.Tag[iens].substr(1,1) == "C")  {
      if(MESON == "B1s") {Corr.Tmin=27; Corr.Tmax=36;}
      else if(MESON == "B2s") {Corr.Tmin=27; Corr.Tmax=36;}
      else if(MESON=="B3s") {Corr.Tmin=18; Corr.Tmax=31;}
      else if(MESON=="B4s") { Corr.Tmin=18; Corr.Tmax=30;}
      else {Corr.Tmin=29; Corr.Tmax=42;}
    }
    else if(data_2pts_SM.Tag[iens].substr(1,1) == "D")  {
      if(MESON == "B1s") {Corr.Tmin=32; Corr.Tmax=53;}
      else if(MESON == "B2s") { Corr.Tmin=28; Corr.Tmax= 37;}
      else if(MESON=="B3s") {Corr.Tmin=28; Corr.Tmax=37;}
      else if(MESON=="B4s") { Corr.Tmin=27; Corr.Tmax=37;}
      else {Corr.Tmin=38; Corr.Tmax=52;}
    }
    else crash("In fixing [Tmin, Tmax] for smeared MP, Ensemble: "+data_2pts_SM.Tag[iens]+" not recognized");


    //mass and decay constants
    distr_t_list pt2_distr_SM= Corr.corr_t(data_2pts_SM.col(0)[iens], path_out+"/corr_2pts/"+data_2pts_SM.Tag[iens]+"/corr_2pt_SM.dat");
    distr_t_list eff_mass_SM = Corr.effective_mass_t(pt2_distr_SM, path_out+"/masses/"+data_2pts_SM.Tag[iens]+"/eff_mass_SM.dat");
    distr_t_list pt2_distr_SMSM= Corr.corr_t(data_2pts_SMSM.col(0)[iens], path_out+"/corr_2pts/"+data_2pts_SM.Tag[iens]+"/corr_2pt_SMSM.dat");
    distr_t_list eff_mass_SMSM = Corr.effective_mass_t(pt2_distr_SMSM, path_out+"/masses/"+data_2pts_SM.Tag[iens]+"/eff_mass_SMSM.dat");
    distr_t mel_SMSM= Corr.Fit_distr(Corr.mel_ov_mass_t(pt2_distr_SMSM, ""))/2.0;

    distr_t_list pt2_distr_SM_boot= Corr_boot.corr_t(data_2pts_SM.col(0)[iens], "");
    distr_t_list eff_mass_SM_boot = Corr_boot.effective_mass_t(pt2_distr_SM_boot, "");
    distr_t_list pt2_distr_SMSM_boot= Corr_boot.corr_t(data_2pts_SMSM.col(0)[iens], "");
    distr_t_list eff_mass_SMSM_boot = Corr_boot.effective_mass_t(pt2_distr_SMSM_boot,"");


    distr_t M_P_SM = Corr.Fit_distr(eff_mass_SM);
    distr_t M_P_SMSM = Corr.Fit_distr(eff_mass_SMSM);
    distr_t M_P= 0.5*(M_P_SM+M_P_SMSM);
    auto SINH= [](double x) { return sinh(x);};
    distr_t_list FP_SM_distr_list = (mu + md )*Corr.residue_t( pt2_distr_SM, "")/(M_P*distr_t::f_of_distr(SINH, M_P)*Corr.matrix_element_t(pt2_distr_SMSM, ""));
    Print_To_File({}, {FP_SM_distr_list.ave(), FP_SM_distr_list.err()}, path_out+"/decay_const/"+data_2pts_SM.Tag[iens]+"/decay_const_SM.dat.t", "", "");
    distr_t FP_SM= Corr.Fit_distr( FP_SM_distr_list  );

    distr_t M_P_boot = 0.5*( Corr_boot.Fit_distr( eff_mass_SMSM_boot) + Corr_boot.Fit_distr(eff_mass_SM_boot));


    cout<<"MP: "<<(M_P/a_distr).ave()<<" +- "<<(M_P/a_distr).err()<<endl;

    
    distr_t F_P= FP_SM;

    //loop over photon momentum
    for(int ixg=0;ixg<n_xg;ixg++) {

      //get xg, Eg, kz from thetas
      double theta=thetas[ixg];
      pt3_momenta pt3_mom_07(0.0, 0.0, thetas[ixg]/2.0, masses_u[ixg], masses_d[ixg], 0.0, L_info.L, L_info.T);
      double Eg= pt3_mom_07.Egamma();
      distr_t Eg_off = M_P - Eg; 
      double kz = pt3_mom_07.k()[2];
      distr_t xg= pt3_mom_07.x_gamma(M_P);
      distr_t xg_boot= pt3_mom_07.x_gamma(M_P_boot);
      xg_list.distr_list.push_back(xg);

      
    
      cout<<"##### Considering kinematic with..."<<endl;
      cout<<"Eg (on-shell): "<<Eg<<endl;
      cout<<"Eg (virt): "<<Eg_off.ave()<<" +- "<<Eg_off.err()<<endl;
      cout<<"thz: "<<theta<<endl;
      cout<<"kz: "<<kz<<endl;
      cout<<"xg: "<<xg.ave()<<" +- "<<xg.err()<<endl;
      int Im_Re;

      //tensor- electric part
      Corr.Reflection_sign= 1;
      Im_Re=1;
      Corr.Perform_Nt_t_average = 0;
      Corr_boot.Reflection_sign=1;
      Corr_boot.Perform_Nt_t_average=0;
	 
	 
      distr_t_list T_u = 0.5*QU*Corr.corr_t(summ_master(C_T_u_data[1][2][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_T_u_data[2][1][ixg].col(Im_Re)[iens], -1.0)), path_out+"/corr_3pts/"+data_2pts_SM.Tag[iens]+"_T_u_xg_"+to_string(ixg));
      distr_t_list T_d = 0.5*QD*Corr.corr_t(summ_master(C_T_d_data[1][2][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_T_d_data[2][1][ixg].col(Im_Re)[iens], -1.0)), path_out+"/corr_3pts/"+data_2pts_SM.Tag[iens]+"_T_d_xg_"+to_string(ixg));
      distr_t_list T_d_boot= 0.5*QD*Corr_boot.corr_t(summ_master(C_T_d_data[2][1][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_T_d_data[2][1][ixg].col(Im_Re)[iens], -1.0)),"");

      Im_Re=0;

      distr_t_list B_u = 0.5*QU*Corr.corr_t(summ_master(C_B_u_data[1][1][ixg].col(Im_Re)[iens], C_B_u_data[2][2][ixg].col(Im_Re)[iens]), path_out+"/corr_3pts/"+data_2pts_SM.Tag[iens]+"_B_u_xg_"+to_string(ixg));
      distr_t_list B_d = 0.5*QD*Corr.corr_t(summ_master(C_B_d_data[1][1][ixg].col(Im_Re)[iens], C_B_d_data[2][2][ixg].col(Im_Re)[iens]), path_out+"/corr_3pts/"+data_2pts_SM.Tag[iens]+"_B_d_xg_"+to_string(ixg));
      distr_t_list B_d_boot= 0.5*QD*Corr_boot.corr_t(summ_master(C_B_d_data[1][1][ixg].col(Im_Re)[iens], C_B_d_data[2][2][ixg].col(Im_Re)[iens]),"");


     
      //determine the form factors that do not need spectral-reconstruction techniques, i.e. up-type quark (heavy) contribution and down-type quark (s-quark) contribution in I-TO


      
      
      distr_t F_T_u(UseJack, Njacks);
      distr_t F_T_d_I(UseJack, Njacks);
      distr_t FV_T_u_real(UseJack, Njacks);
      distr_t FV_T_d_real(UseJack, Njacks);
      distr_t FA_T_u_real(UseJack, Njacks);
      distr_t FA_T_d_real(UseJack, Njacks);

      
      auto HeavyTheta=[](const int x) {  return ((x>=0)+(x>0))/2.0;   };
      auto Exp= [&Njacks, &UseJack](const distr_t &A) -> distr_t { distr_t ret(UseJack); for(int ijack=0;ijack<Njacks;ijack++) ret.distr.push_back( exp(A.distr[ijack])); return ret;};
      double T=Corr.Nt;
      
      for(int ty=0; ty < Corr.Nt; ty++) {
	
	const distr_t f1=  Exp(-(T/2-ty)*Eg_off);

	const distr_t f2= Exp(-((3*T/2)-ty)*Eg_off);

	const double f1_real= exp(-(T/2-ty)*Eg);
	const double f2_real= exp(-((3*T/2)-ty)*Eg);
	
	const double h1=HeavyTheta((T/2)-ty);
	const double h2=HeavyTheta(ty-(T/2));



		
	F_T_u = F_T_u +  (sign_kz*T_u.distr_list[ty]*(xg/2.0)   + B_u.distr_list[ty]*xg/2.0 )*(h1*f1+ h2*f2);
	FV_T_u_real = FV_T_u_real + (sign_kz*T_u.distr_list[ty]*(1.0- xg/2.0)   + B_u.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FV_T_d_real = FV_T_d_real + (sign_kz*T_d.distr_list[ty]*(1.0- xg/2.0)   + B_d.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FA_T_u_real = FA_T_u_real + ( B_u.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_u.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	FA_T_d_real = FA_T_d_real + ( B_d.distr_list[ty]*(1.0- xg/2.0)   + sign_kz*T_d.distr_list[ty]*xg/2.0 )*(h1*f1_real+ h2*f2_real);
	if(ty<= t_07_s) {
	  F_T_d_I = F_T_d_I + (sign_kz*T_d.distr_list[ty]*(xg/2.0)   + B_d.distr_list[ty]*xg/2.0 )*(h1*f1+ h2*f2);
	}
      }

      
         
      //normalize
      F_T_u = F_T_u*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_c))*Exp( M_P*t_07_c) ;
      F_T_d_I = F_T_d_I*ZT*(1.0/(mel_SMSM*Eg))*Exp( Eg_off*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      FV_T_u_real= FV_T_u_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c) ;
      FV_T_d_real= FV_T_d_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      FA_T_u_real= FA_T_u_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_c))*Exp( M_P*t_07_c) ;
      FA_T_d_real= FA_T_d_real*ZT*(1.0/(mel_SMSM*Eg))*exp( Eg*abs(T/2 - t_07_s))*Exp( M_P*t_07_s) ;
      //push_back
      F_T_u_list[ixg].distr_list.push_back( F_T_u);
      F_T_d_I_list[ixg].distr_list.push_back( F_T_d_I);
      FV_T_u_real_list[ixg].distr_list.push_back( FV_T_u_real);
      FV_T_d_real_list[ixg].distr_list.push_back( FV_T_d_real);
      FA_T_u_real_list[ixg].distr_list.push_back( FA_T_u_real);
      FA_T_d_real_list[ixg].distr_list.push_back( FA_T_d_real);
      

      //define corrs to be used in spec dens reconstruction
      distr_t_list Corr_T(UseJack);
      distr_t_list Corr_T_boot(UseJack);

      distr_t_list Corr_TV_real(UseJack), Corr_TV_real_boot(UseJack);

      for(int t=t_07_s_HLT; t <= Corr.Nt/2; t++) {
	
	Corr_T.distr_list.push_back( sign_kz*T_d.distr_list[t] + B_d.distr_list[t]);
	Corr_T_boot.distr_list.push_back( sign_kz*T_d_boot.distr_list[t] + B_d_boot.distr_list[t]);

	Corr_TV_real.distr_list.push_back(   sign_kz*T_d.distr_list[t]*(1 - xg/2.0) + B_d.distr_list[t]*(xg/2.0));
	Corr_TV_real_boot.distr_list.push_back(   sign_kz*T_d_boot.distr_list[t]*(1 - xg_boot/2.0) + B_d_boot.distr_list[t]*(xg_boot/2.0));
	
      }

      int tmax= Corr_T.size();

      distr_t FACT= ZT*(1.0/(mel_SMSM*Eg) )*Exp( M_P*t_07_s_HLT)*xg/2.0;
      distr_t FACT_real = ZT*(1.0/(mel_SMSM*Eg) )*Exp( M_P*t_07_s_HLT);

      
      CorrAnalysis Corr_HLT(UseJack, Njacks,1000);
      Corr_HLT.Tmin= 20;
      Corr_HLT.Tmax= 30;
      Corr_HLT.Nt=2*(Corr_T.size()-1);
      distr_t_list Corr_T_symm= Corr_T;
      for(int t=Corr_HLT.Nt/2 +1;t<Corr_HLT.Nt;t++) Corr_T_symm.distr_list.push_back( Corr_T.distr_list[Corr_HLT.Nt -t]);
      distr_t eff_M_T= Corr_HLT.Fit_distr(Corr_HLT.effective_mass_t(Corr_T_symm, "" ))/a_distr;
      double Mphi_motion= sqrt( Mphi*Mphi + pow(kz/a_distr.ave(),2));
     
    

      cout<<"eff_M_T: "<<eff_M_T.ave()<<" +- "<<eff_M_T.err()<<" expected: "<<Mphi_motion<<endl;
 
      //generate covariance matrix
      Vfloat cov_T, corr_T;
      Vfloat cov_TV_real, corr_TV_real;
      Vfloat TT, RR;
      for(int tt=0;tt< tmax;tt++)
	for(int rr=0;rr< tmax;rr++) {
	  TT.push_back(tt);
	  RR.push_back(rr);
	  cov_T.push_back( Corr_T_boot.distr_list[tt]%Corr_T_boot.distr_list[rr]);
	  corr_T.push_back( ( Corr_T_boot.distr_list[tt]%Corr_T_boot.distr_list[rr])/( Corr_T_boot.err(tt)*Corr_T_boot.err(rr)));

	  cov_TV_real.push_back( Corr_TV_real_boot.distr_list[tt]%Corr_TV_real_boot.distr_list[rr]);
	  corr_TV_real.push_back( ( Corr_TV_real_boot.distr_list[tt]%Corr_TV_real_boot.distr_list[rr])/( Corr_TV_real_boot.err(tt)*Corr_TV_real_boot.err(rr)));
	  
	}
      

      //print covariance matrix
      Print_To_File({},{TT,RR, cov_T, corr_T}, path_out+"/covariance/"+Ens_tags[iens]+"/cov_T_xg_"+to_string_with_precision(xg.ave(),2)+".cov", "" , "");
  

      cout<<"Starting spectral reconstruction:..."<<endl;
      

      if(!Skip_spectral_reconstruction_07) {

	#pragma omp parallel for schedule(dynamic)
	//spectral reconstruction for second time ordering
	for(int isg=0;isg<(signed)sigmas_07.size();isg++) {

	  double s= sigmas_07[isg]*a_distr.ave();
	  
	  double syst_T;
	  double mult_T= 0.5;
	  double Ag_target=1e-4;
	  if(sigmas_07[isg] < 0.5) Ag_target=5e-2;
	  else if(sigmas_07[isg] < 0.7) Ag_target=5e-3;
	  double th= E0_fact*Mphi_motion;
	  double l_re_T;
	  
	  
	  F_T_d_RE_sm_list[ixg][isg].distr_list.push_back( Get_Laplace_transfo ( Eg_off.ave(),  s, th*a_distr.ave(),  Nts[iens], tmax-1, prec_07, SM_TYPE+"_RE",K_RE, Corr_T, syst_T, mult_T ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"T_"+Ens_tags[iens], Ag_target,0, FACT , "07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0,1));
	  	  
	  F_T_d_IM_sm_list[ixg][isg].distr_list.push_back( Get_Laplace_transfo ( Eg_off.ave(),  s, th*a_distr.ave(),  Nts[iens], tmax-1, prec_07, SM_TYPE+"_IM",K_IM, Corr_T, syst_T, mult_T ,  l_re_T, MODE_FF, "E0_"+to_string_with_precision(E0_fact,1), TAG_CURR+"T_"+Ens_tags[iens], Ag_target,0, FACT , "07_FF_Tw_"+to_string(t_07_s_HLT), cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0,1));
	  
	 	  
	  
	}
      }
            
      
    }

  }

  //print results
  //per kinematic
  for(int ixg=0;ixg<n_xg;ixg++) {
    Print_To_File( Ens_tags, { F_T_u_list[ixg].ave(), F_T_u_list[ixg].err()     } , path_out+"/FF_u/F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { F_T_d_I_list[ixg].ave(), F_T_d_I_list[ixg].err()     } , path_out+"/FF_d_I/F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FV_T_u_real_list[ixg].ave(), FV_T_u_real_list[ixg].err()     } , path_out+"/FF_u/FV_T_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FV_T_d_real_list[ixg].ave(), FV_T_d_real_list[ixg].err()     } , path_out+"/FF_d/FV_T_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FA_T_u_real_list[ixg].ave(), FA_T_u_real_list[ixg].err()     } , path_out+"/FF_u/FA_T_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    Print_To_File( Ens_tags, { FA_T_d_real_list[ixg].ave(), FA_T_d_real_list[ixg].err()     } , path_out+"/FF_d/FA_T_real_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");
    
    for(int isg=0;isg<(signed)sigmas_07.size();isg++) {

      Print_To_File( Ens_tags, { F_T_d_RE_sm_list[ixg][isg].ave(), F_T_d_RE_sm_list[ixg][isg].err(), F_T_d_IM_sm_list[ixg][isg].ave(), F_T_d_IM_sm_list[ixg][isg].err()     } , path_out+"/FF_d_II/F_T_sm_"+to_string_with_precision(sigmas_07[isg],3)+"_xg_"+to_string_with_precision(xg_list.ave(ixg),2), "", "");


    }
  }
  //per ensemble
  for(int iens=0; iens<Nens;iens++) {
    distr_t_list F_T_u_per_ens(UseJack),  F_T_d_I_per_ens(UseJack);
    distr_t_list FV_T_u_real_per_ens(UseJack),  FV_T_d_real_per_ens(UseJack);
    distr_t_list FA_T_u_real_per_ens(UseJack),  FA_T_d_real_per_ens(UseJack);
   

    for(int ixg=0;ixg<n_xg;ixg++) {
      F_T_u_per_ens.distr_list.push_back( F_T_u_list[ixg].distr_list[iens]);
      F_T_d_I_per_ens.distr_list.push_back( F_T_d_I_list[ixg].distr_list[iens]);
      FV_T_u_real_per_ens.distr_list.push_back( FV_T_u_real_list[ixg].distr_list[iens]);
      FV_T_d_real_per_ens.distr_list.push_back( FV_T_d_real_list[ixg].distr_list[iens]);
      FA_T_u_real_per_ens.distr_list.push_back( FA_T_u_real_list[ixg].distr_list[iens]);
      FA_T_d_real_per_ens.distr_list.push_back( FA_T_d_real_list[ixg].distr_list[iens]);

      distr_t_list F_T_d_RE_sm_per_ens_per_kin, F_T_d_IM_sm_per_ens_per_kin;
      
      for(int isg=0;isg<(signed)sigmas_07.size(); isg++) {
	F_T_d_RE_sm_per_ens_per_kin.distr_list.push_back( F_T_d_RE_sm_list[ixg][isg].distr_list[iens]);
	F_T_d_IM_sm_per_ens_per_kin.distr_list.push_back( F_T_d_IM_sm_list[ixg][isg].distr_list[iens]);
      }
      
      Print_To_File( {}, {sigmas_07, F_T_d_RE_sm_per_ens_per_kin.ave(), F_T_d_RE_sm_per_ens_per_kin.err(), F_T_d_IM_sm_per_ens_per_kin.ave(), F_T_d_IM_sm_per_ens_per_kin.err()}, path_out+"/FF_d_II/F_T_xg_"+to_string_with_precision(xg_list.ave(ixg),2)+"_"+Ens_tags[iens], "", "");
      
    }

    Print_To_File( {}, {xg_list.ave(),  F_T_u_per_ens.ave(), F_T_u_per_ens.err()     } , path_out+"/FF_u/F_T_"+Ens_tags[iens], "", "");
    Print_To_File( {}, {xg_list.ave(), F_T_d_I_per_ens.ave(), F_T_d_I_per_ens.err()     } , path_out+"/FF_d_I/F_T_"+Ens_tags[iens], "", "");

    Print_To_File( {}, {xg_list.ave(),  FV_T_u_real_per_ens.ave(), FV_T_u_real_per_ens.err()     } , path_out+"/FF_u/FV_T_real_"+Ens_tags[iens], "", "");
    Print_To_File( {}, {xg_list.ave(), FV_T_d_real_per_ens.ave(), FV_T_d_real_per_ens.err()     } , path_out+"/FF_d/FV_T_real_"+Ens_tags[iens], "", "");

    Print_To_File( {}, {xg_list.ave(),  FA_T_u_real_per_ens.ave(), FA_T_u_real_per_ens.err()     } , path_out+"/FF_u/FA_T_real_"+Ens_tags[iens], "", "");
    Print_To_File( {}, {xg_list.ave(), FA_T_d_real_per_ens.ave(), FA_T_d_real_per_ens.err()     } , path_out+"/FF_d/FA_T_real_"+Ens_tags[iens], "", "");
    
    for(int isg=0;isg<(signed)sigmas_07.size(); isg++) {

      distr_t_list F_T_d_RE_sm_per_ens_per_sigma, F_T_d_IM_sm_per_ens_per_sigma;
      
       for(int ixg=0;ixg<n_xg;ixg++) {
	 F_T_d_RE_sm_per_ens_per_sigma.distr_list.push_back( F_T_d_RE_sm_list[ixg][isg].distr_list[iens]);
	 F_T_d_IM_sm_per_ens_per_sigma.distr_list.push_back( F_T_d_IM_sm_list[ixg][isg].distr_list[iens]);
       }

       Print_To_File( {}, {xg_list.ave(), F_T_d_RE_sm_per_ens_per_sigma.ave(), F_T_d_RE_sm_per_ens_per_sigma.err(), F_T_d_IM_sm_per_ens_per_sigma.ave(), F_T_d_IM_sm_per_ens_per_sigma.err()}, path_out+"/FF_d_II/F_T_sm_"+to_string_with_precision(sigmas_07[isg],3)+"_"+Ens_tags[iens], "", "");
       Print_To_File( {}, {xg_list.ave(), (F_T_d_RE_sm_per_ens_per_sigma + F_T_u_per_ens+ F_T_d_I_per_ens).ave(), (F_T_d_RE_sm_per_ens_per_sigma + F_T_u_per_ens+ F_T_d_I_per_ens).err(),  F_T_d_IM_sm_per_ens_per_sigma.ave(), F_T_d_IM_sm_per_ens_per_sigma.err()}, path_out+"/FF/F_T_sm_"+to_string_with_precision(sigmas_07[isg],3)+"_"+Ens_tags[iens], "", "");
       
       
    }
  }

    
  


  

 





  return return_class;

}
