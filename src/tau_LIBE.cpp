#include "../include/tau_LIBE.h"
#include "Corr_analysis.h"
#include "numerics.h"
using namespace std;


const int Nboots = 800;
const bool UseJack=true;
//const int Njacks = 606;
const vector<vector<string>> muval({ {"7.2000e-04","1.8250e-02"},  { "6.0000e-04", "1.6044e-02"},    {"5.4000e-04", "1.3576e-02"} } );
const double alpha = 1.0 / 137.035999;
const vector<int> Nconfs_per_ens({605,321,397});
const double qu = 2.0 / 3;
const double qd = -qu / 2;
const double qs= qd;
const double MK_plus = 0.493677;
const double MK_zero = 0.497611;
const double DMK = -0.003934;
const double kappa = 2.837297;
const double m_tau = 1.77686;
const double C_V = 2*M_PI/(pow(m_tau,3));
const double GAMMA_FACT= 12*M_PI; 
const vector<string> Ens_list({"cB211b.072.64", "cC211a.06.80", "cD211a.054.96"});
enum { V0P5, P5P5, VkVk, AkAk, A0A0, V0V0 };
enum { ISO_Q, EXCH, SELF_U, SELF_D };
enum { ISO_S, MCR_U, MCR_D, M_U, M_D};
const vector<string> Corr_Tags({"V0P5", "P5P5", "VkVk", "AkAk", "A0A0", "V0V0"}); //ordering should be the same as in enum
const vector<string> QED_Tags({"ISO_Q", "EXCH", "SELF_U", "SELF_D"});  //ordering here should be the same as in ENUM
const vector<string> SIB_Tags({"ISO_S", "MCR_U", "MCR_D" , "M_U" , "M_D"});  //ordering here should be the same as in ENUM


void Read_file_QED(string path, VVVfloat &QED, int reim) {

  FILE *stream= fopen(path.c_str(), "rb");
  if(stream == NULL ) crash("File "+path+" not open correctly");
  int TT ,Nconfs,  Nsubs;
  bin_read(TT, stream);
  bin_read(Nconfs, stream);
  bin_read(Nsubs, stream);

  assert(QED.size() == 4);
  for(int k=0;k<(signed)QED.size();k++) {
    QED[k].resize(TT);
    for(int t=0;t<TT;t++) {
      QED[k][t].resize(Nconfs,0.0);
    }
  }


  cout<<"Reading: "<<path<<endl;
  cout<<"(TT,Nconfs,Nsubs): "<<"("<<TT<<","<<Nconfs<<","<<Nsubs<<")"<<endl;

  
  for(int i=0;i<5;i++) {
    double fact=1.0;
    if(i==1 || i==2) fact=0.5;
    int what;
    bin_read(what,stream);
    //cout<<"what?: "<<what<<endl;
    for(int ic=0;ic<Nconfs;ic++) {
      for(int isub=0;isub<Nsubs;isub++) {
	for(int t=0;t<TT/2+1;t++) {
	  double re,im;
	  bin_read(re,stream); bin_read(im,stream);
	  QED[(i>=2)?(i-1):i][t][ic] += fact*(reim?re:im)/Nsubs;
	  if(t!= 0) QED[(i>=2)?(i-1):i][TT-t][ic] = QED[(i>=2)?(i-1):i][t][ic];
	}
      }
    }
  }

  cout<<"done!"<<endl;
   
  return;
}



void Read_file_SIB(string path, VVVfloat &SIB,  int r1, int r2, int reim) {


  FILE *stream= fopen(path.c_str(), "rb");
  if(stream == NULL ) crash("File "+path+" not open correctly");
  int TT ,Nconfs,  Nsubs;
  bin_read(TT, stream);
  bin_read(Nconfs, stream);
  bin_read(Nsubs, stream);

  assert(SIB.size() == 5);
  for(int k=0;k<(signed)SIB.size();k++) {
    SIB[k].resize(TT);
    for(int t=0;t<TT;t++) {
      SIB[k][t].resize(Nconfs,0.0);
    }
  }


  cout<<"Reading: "<<path<<endl;
  cout<<"(TT,Nconfs,Nsubs): "<<"("<<TT<<","<<Nconfs<<","<<Nsubs<<")"<<endl;


  //for mass insertion I still have to multiply by i*r , where r is the twisted Wilson term of the quark propagator with the scalar insertion


  for(int i=0;i<5;i++) {
    double fact=1.0;
    //double fact=(i==0)?1.0:-1.0;
    int what;
    bin_read(what,stream);
    //cout<<"what?: "<<what<<endl;
    for(int ic=0;ic<Nconfs;ic++) {
      for(int isub=0;isub<Nsubs;isub++) {
	for(int t=0;t<TT/2+1;t++) {
	  double re,im;
	  bin_read(re,stream); bin_read(im,stream);
	  
	  SIB[i][t][ic] += fact*(reim?re:im)/Nsubs;

	  //if(i<=2) SIB[i][t][ic] += fact*(reim?re:im)/Nsubs;
	  //else if(i==3) SIB[i][t][ic] += (reim?-1:1)*r2*fact*(reim?im:re)/Nsubs;
	  //else if(i==4) SIB[i][t][ic] += (reim?-1:1)*r1*fact*(reim?im:re)/Nsubs;
	    
	  if(t!= 0) SIB[i][TT-t][ic] = SIB[i][t][ic];
	}
      }
    }
  }

   
  return;
}


distr_t_list Get_LIBE_correlator(const vector<distr_t_list>& C_QED, const vector<distr_t_list>& C_SIB, const FLAV& U, const FLAV& D, string path) {

  distr_t_list Corr_QED =   4.0*M_PI*alpha*( C_QED[SELF_U]*U.q*U.q + C_QED[SELF_D]*D.q*D.q + C_QED[EXCH]*U.q*D.q );
  distr_t_list Corr_SIB =   C_SIB[MCR_U]*U.dmcr + C_SIB[MCR_D]*D.dmcr + C_SIB[M_U]*U.dm + C_SIB[M_D]*D.dm;

  if(path != "") {
    distr_t_list C= Corr_QED+Corr_SIB;
    Print_To_File({}, {C.ave(), C.err()}, path, "","");
  }
  
  return Corr_QED + Corr_SIB;
}

distr_t_list Get_LIBE_correlator(const vector<distr_t_list>& C_QED, const vector<distr_t_list>& C_SIB, const FLAV& U, const FLAV& D, const distr_t& F,  string path) {


  distr_t_list C = Get_LIBE_correlator(C_QED, C_SIB, U,D, "");

  if(path != "") {
    Print_To_File({}, {(F*C).ave(), (F*C).err()}, path, "","");
  }
  
  return F*C;
}




void Compute_tau_LIBE() {


  
  double s;
  
  double fmTGeV= 1.0/0.197327;

  LatticeInfo L_info_A, L_info_B, L_info_C, L_info_D;
  L_info_A.LatInfo_new_ens("cA211a.12.48");
  L_info_B.LatInfo_new_ens("cB211b.072.96");
  L_info_C.LatInfo_new_ens("cC211a.06.80");
  L_info_D.LatInfo_new_ens("cD211a.054.96");

  distr_t_list A0_distr_TM(UseJack);
  distr_t_list A0_distr_OS(UseJack);

  distr_t_list V0_distr_TM(UseJack);
  distr_t_list V0_distr_OS(UseJack);

  distr_t_list Ak_distr_TM(UseJack);
  distr_t_list Ak_distr_OS(UseJack);
  
  distr_t_list Vk_distr_TM(UseJack);
  distr_t_list Vk_distr_OS(UseJack);

  distr_t_list a_distr_list(UseJack);

  
  int Nens= Ens_list.size();

  boost::filesystem::create_directory("../data/tau_LIBE");
  

  //loop over ensembles
  for(int iens=0;iens<Nens;iens++) {



    GaussianMersenne GM(78821);

    int Njacks= Nconfs_per_ens[iens];
    
    //resample RCs
    distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
    distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
    distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);
    

    for(int ijack=0; ijack<Njacks;ijack++) {
      
      ZA_A.distr.push_back( L_info_A.Za_WI_strange + GM()*L_info_A.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
      ZV_A.distr.push_back( L_info_A.Zv_WI_strange + GM()*L_info_A.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
      
      ZA_B.distr.push_back( L_info_B.Za_WI_strange + GM()*L_info_B.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
      ZV_B.distr.push_back( L_info_B.Zv_WI_strange + GM()*L_info_B.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
      
      ZA_C.distr.push_back( L_info_C.Za_WI_strange + GM()*L_info_C.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
      ZV_C.distr.push_back( L_info_C.Zv_WI_strange + GM()*L_info_C.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
      
      ZA_D.distr.push_back( L_info_D.Za_WI_strange + GM()*L_info_D.Za_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
      ZV_D.distr.push_back( L_info_D.Zv_WI_strange + GM()*L_info_D.Zv_WI_strange_err/((UseJack==true)?sqrt(Njacks -1.0):1.0));
      
      a_A.distr.push_back( L_info_A.a_from_afp*fmTGeV + GM()*L_info_A.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
      a_B.distr.push_back( L_info_B.a_from_afp*fmTGeV + GM()*L_info_B.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
      a_C.distr.push_back( L_info_C.a_from_afp*fmTGeV + GM()*L_info_C.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
      a_D.distr.push_back( L_info_D.a_from_afp*fmTGeV + GM()*L_info_D.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    }
    
    cout<<"RC generated!"<<endl;
    
    string Ens_tag= Ens_list[iens];
    cout<<"Analyzing ensemble: "<<Ens_tag<<endl;
    boost::filesystem::create_directory("../data/tau_LIBE/"+Ens_tag);
    boost::filesystem::create_directory("../data/tau_LIBE/"+Ens_tag+"/corr");
    boost::filesystem::create_directory("../data/tau_LIBE/"+Ens_tag+"/eff_slope");
    boost::filesystem::create_directory("../data/tau_LIBE/"+Ens_tag+"/mass");
    boost::filesystem::create_directory("../data/tau_LIBE/"+Ens_tag+"/REN");
      
    
    //Lattice info
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(Ens_tag);
    CorrAnalysis Corr(UseJack, Njacks,Nboots, iens);
    Corr.Nt = L_info.T;
    int L = L_info.L;
    int T= 2*L;
    double V= pow(L_info.L,3);

    double ams=stod(muval[iens][1]);

    auto mass_shift = [&Corr](const distr_t_list &A, const distr_t_list &B, string path) -> distr_t_list { return -1.0*Corr.effective_slope_t(A,B, path);};
    

    
    //find ensemble info
    distr_t Za, Zv, a_distr;
    if(Ens_tag.substr(1,1)=="A") { Za= ZA_A; Zv=ZV_A; a_distr=a_A;}
    else if(Ens_tag.substr(1,1)=="B") { Za= ZA_B; Zv=ZV_B; a_distr=a_B;}
    else if(Ens_tag.substr(1,1)=="C") { Za= ZA_C; Zv=ZV_C; a_distr=a_C;}
    else if(Ens_tag.substr(1,1)=="D") { Za= ZA_D; Zv=ZV_D; a_distr=a_D;}
    else crash("Ensemble: "+Ens_tag+" not recognised");
    
    
    vector<vector<distr_t_list>> C_TM_ll_QED(Corr_Tags.size()), C_TM_ll_SIB(Corr_Tags.size());
    vector<vector<distr_t_list>> C_OS_ll_QED(Corr_Tags.size()), C_OS_ll_SIB(Corr_Tags.size());
    vector<vector<distr_t_list>> C_TM_ls_QED(Corr_Tags.size()), C_TM_ls_SIB(Corr_Tags.size());
    vector<vector<distr_t_list>> C_OS_ls_QED(Corr_Tags.size()), C_OS_ls_SIB(Corr_Tags.size());
    vector<vector<distr_t_list>> C_TM_ss_QED(Corr_Tags.size()), C_TM_ss_SIB(Corr_Tags.size());
    vector<vector<distr_t_list>> C_OS_ss_QED(Corr_Tags.size()), C_OS_ss_SIB(Corr_Tags.size());
    

        
    //read data from binary files
    for(int i=0; i<(signed)Corr_Tags.size();i++) {
      int k_Q=4;
      int k_SIB=5;

      VVVfloat ll_QED(4), ll_SIB(5), OS_ll_QED(4), OS_ll_SIB(5);
      VVVfloat ls_QED(4), ls_SIB(5), OS_ls_QED(4) , OS_ls_SIB(5);
      VVVfloat ss_QED(4), ss_SIB(5), OS_ss_QED(4) , OS_ss_SIB(5);

      int reim=1;
      if(i==V0P5) reim=0;

          
      //TM_
      Read_file_QED("../tau_LIBE/"+Ens_tag+"/"+Corr_Tags[i]+"_tm/"+muval[iens][0]+"_"+muval[iens][0]+"_QED.bin", ll_QED, reim);
      Read_file_SIB("../tau_LIBE/"+Ens_tag+"/"+Corr_Tags[i]+"_tm/"+muval[iens][0]+"_"+muval[iens][0]+"_SIB.bin", ll_SIB, 1,-1, reim);
      Read_file_QED("../tau_LIBE/"+Ens_tag+"/"+Corr_Tags[i]+"_tm/"+muval[iens][0]+"_"+muval[iens][1]+"_QED.bin", ls_QED, reim);
      Read_file_SIB("../tau_LIBE/"+Ens_tag+"/"+Corr_Tags[i]+"_tm/"+muval[iens][0]+"_"+muval[iens][1]+"_SIB.bin", ls_SIB,1, -1,reim);
      Read_file_QED("../tau_LIBE/"+Ens_tag+"/"+Corr_Tags[i]+"_tm/"+muval[iens][1]+"_"+muval[iens][1]+"_QED.bin", ss_QED, reim);
      Read_file_SIB("../tau_LIBE/"+Ens_tag+"/"+Corr_Tags[i]+"_tm/"+muval[iens][1]+"_"+muval[iens][1]+"_SIB.bin", ss_SIB, 1,-1, reim);
      //OS
      Read_file_QED("../tau_LIBE/"+Ens_tag+"/"+Corr_Tags[i]+"_OS/"+muval[iens][0]+"_-"+muval[iens][0]+"_QED.bin", OS_ll_QED, reim);
      Read_file_SIB("../tau_LIBE/"+Ens_tag+"/"+Corr_Tags[i]+"_OS/"+muval[iens][0]+"_-"+muval[iens][0]+"_SIB.bin", OS_ll_SIB, 1, 1, reim);
      Read_file_QED("../tau_LIBE/"+Ens_tag+"/"+Corr_Tags[i]+"_OS/"+muval[iens][0]+"_-"+muval[iens][1]+"_QED.bin", OS_ls_QED, reim);
      Read_file_SIB("../tau_LIBE/"+Ens_tag+"/"+Corr_Tags[i]+"_OS/"+muval[iens][0]+"_-"+muval[iens][1]+"_SIB.bin", OS_ls_SIB, 1,1,reim);
      Read_file_QED("../tau_LIBE/"+Ens_tag+"/"+Corr_Tags[i]+"_OS/"+muval[iens][1]+"_-"+muval[iens][1]+"_QED.bin", OS_ss_QED, reim);
      Read_file_SIB("../tau_LIBE/"+Ens_tag+"/"+Corr_Tags[i]+"_OS/"+muval[iens][1]+"_-"+muval[iens][1]+"_SIB.bin", OS_ss_SIB, 1, 1, reim);


      for(int Q=0;Q<k_Q;Q++) {
	C_TM_ll_QED[i].push_back( Corr.corr_t( ll_QED[Q], "../data/tau_LIBE/"+Ens_tag+"/corr/ll_"+Corr_Tags[i]+"_TM_"+QED_Tags[Q]));
	C_OS_ll_QED[i].push_back( Corr.corr_t( OS_ll_QED[Q],  "../data/tau_LIBE/"+Ens_tag+"/corr/ll_"+Corr_Tags[i]+"_OS_"+QED_Tags[Q]));
	C_TM_ls_QED[i].push_back( Corr.corr_t( ls_QED[Q],  "../data/tau_LIBE/"+Ens_tag+"/corr/ls_"+Corr_Tags[i]+"_TM_"+QED_Tags[Q]));
	C_OS_ls_QED[i].push_back( Corr.corr_t( OS_ls_QED[Q],  "../data/tau_LIBE/"+Ens_tag+"/corr/ls_"+Corr_Tags[i]+"_OS_"+QED_Tags[Q]));
	C_TM_ss_QED[i].push_back( Corr.corr_t( ss_QED[Q], "../data/tau_LIBE/"+Ens_tag+"/corr/ss_"+Corr_Tags[i]+"_TM_"+QED_Tags[Q]));
	C_OS_ss_QED[i].push_back( Corr.corr_t( OS_ss_QED[Q],  "../data/tau_LIBE/"+Ens_tag+"/corr/ss_"+Corr_Tags[i]+"_OS_"+QED_Tags[Q]));
      }
      for(int S=0;S<k_SIB;S++) {
	C_TM_ll_SIB[i].push_back( Corr.corr_t( ll_SIB[S],  "../data/tau_LIBE/"+Ens_tag+"/corr/ll_"+Corr_Tags[i]+"_TM_"+SIB_Tags[S]));
	C_OS_ll_SIB[i].push_back( Corr.corr_t( OS_ll_SIB[S],  "../data/tau_LIBE/"+Ens_tag+"/corr/ll_"+Corr_Tags[i]+"_OS_"+SIB_Tags[S]));
	C_TM_ls_SIB[i].push_back( Corr.corr_t( ls_SIB[S],  "../data/tau_LIBE/"+Ens_tag+"/corr/ls_"+Corr_Tags[i]+"_TM_"+SIB_Tags[S]));
	C_OS_ls_SIB[i].push_back( Corr.corr_t( OS_ls_SIB[S],  "../data/tau_LIBE/"+Ens_tag+"/corr/ls_"+Corr_Tags[i]+"_OS_"+SIB_Tags[S]));
	C_TM_ss_SIB[i].push_back( Corr.corr_t( ss_SIB[S],  "../data/tau_LIBE/"+Ens_tag+"/corr/ss_"+Corr_Tags[i]+"_TM_"+SIB_Tags[S]));
	C_OS_ss_SIB[i].push_back( Corr.corr_t( OS_ss_SIB[S],  "../data/tau_LIBE/"+Ens_tag+"/corr/ss_"+Corr_Tags[i]+"_OS_"+SIB_Tags[S]));
      }
    }


       
    
    distr_t_list pion_distr= Corr.effective_mass_t( C_TM_ll_QED[P5P5][ISO_Q], "../data/tau_LIBE/"+Ens_tag+"/mass/pion");
    distr_t_list kaon_distr= Corr.effective_mass_t( C_TM_ls_QED[P5P5][ISO_Q], "../data/tau_LIBE/"+Ens_tag+"/mass/kaon");

      if(Ens_tag=="cB211b.072.48") {
	Corr.Tmin=20;
	Corr.Tmax=35;
      }
      else {
	Corr.Tmin=30;
	Corr.Tmax=45;
      }

      distr_t Mpi=Corr.Fit_distr(pion_distr);
      distr_t M_K= Corr.Fit_distr(kaon_distr);
      distr_t aDMpi_plus= 0.13957039*a_distr - 0.135*a_distr ;
      distr_t Delta_MK_iso=  (MK_plus*a_distr - 0.4946*a_distr);
      distr_t aDMK= DMK*a_distr; 

      distr_t FSE_charged_pion= alpha*(2*kappa/(2.0*Mpi*pow(L,2)))*( 1.0 + 0.5*L*Mpi);
      distr_t FSE_charged_kaon= alpha*(2*kappa/(2.0*M_K*pow(L,2)))*( 1.0 + 0.5*L*M_K);  
      
            
      //#################################################
      //compute pion mass splitting
      distr_t_list Delta_pion_distr= Corr.effective_slope_t( (4*M_PI)*0.5*alpha*C_TM_ll_QED[P5P5][EXCH], C_TM_ll_QED[P5P5][ISO_Q] , "../data/tau_LIBE/"+Ens_tag+"/eff_slope/pion_splitting" );
      distr_t_list Delta_pion_subL_distr= Delta_pion_distr + alpha*(2*kappa/(2.0*Mpi*pow(L,2)))*( 1.0 + 0.5*L*Mpi);
     
      Print_To_File({}, {Delta_pion_distr.ave(), Delta_pion_distr.err(), Delta_pion_subL_distr.ave(), Delta_pion_subL_distr.err()}, "../data/tau_LIBE/"+Ens_tag+"/eff_slope/pion_splitting_subL", "", "");
      //#################################################
      

      cout<<"### PRINTING COUNTERTERMS ###"<<endl;

      //#################################################
      //compute critical mass shift
      distr_t_list VP_QED= 4.0*M_PI*alpha*( distr_t_list::derivative((qu*qu+qd*qd)*C_TM_ll_QED[V0P5][SELF_U],0)/distr_t_list::derivative(C_TM_ll_QED[V0P5][ISO_Q],0) -  ((qu*qu+qd*qd)*C_TM_ll_QED[P5P5][SELF_U])/C_TM_ll_QED[P5P5][ISO_Q]);
      distr_t_list VP_DM = distr_t_list::derivative(C_TM_ll_SIB[V0P5][MCR_U],0)/distr_t_list::derivative(C_TM_ll_SIB[V0P5][ISO_S],0) - C_TM_ll_SIB[P5P5][MCR_U]/C_TM_ll_SIB[P5P5][ISO_S];
      distr_t_list dmcr_distr = -1.0*VP_QED/VP_DM;
      Print_To_File( {}, {dmcr_distr.ave(), dmcr_distr.err()}, "../data/tau_LIBE/"+Ens_tag+"/REN/dmcr.eff", "", "");
      if(Ens_tag=="cB211b.072.48") {Corr.Tmin = 12; Corr.Tmax = 25;}
      else {Corr.Tmin = 18; Corr.Tmax = 30; }
      distr_t dmcr = Corr.Fit_distr( dmcr_distr);
      distr_t dmcr_u = (4.0/5.0)*dmcr;
      distr_t dmcr_d = (1.0/5.0)*dmcr;
      distr_t dmcr_s = dmcr_d;
      cout<<"dmcr:        "<<dmcr.ave()/2<<" "<<dmcr.err()/2<<endl;
      //#################################################


      
      //#################################################
      //compute average ud mass shift from Mpi^+    DM = (dmu+dmud)/2
      distr_t_list DM_distr = 4.0*M_PI*alpha*(   (qu*qu+qd*qd)*mass_shift( C_TM_ll_QED[P5P5][SELF_U],C_TM_ll_QED[P5P5][ISO_Q], "")); //self
      DM_distr = DM_distr + 4.0*M_PI*alpha*qu*qd*mass_shift( C_TM_ll_QED[P5P5][EXCH], C_TM_ll_QED[P5P5][ISO_Q] , ""); //exchange
      DM_distr = DM_distr + dmcr*mass_shift( C_TM_ll_SIB[P5P5][MCR_U], C_TM_ll_SIB[P5P5][ISO_S],  ""); //critical mass
      DM_distr = aDMpi_plus -FSE_charged_pion - DM_distr; 
      DM_distr = 0.5*DM_distr/mass_shift( C_TM_ll_SIB[P5P5][M_U], C_TM_ll_SIB[P5P5][ISO_S], ""); //solve equation
      Print_To_File( {}, {DM_distr.ave(), DM_distr.err()}, "../data/tau_LIBE/"+Ens_tag+"/REN/DM.eff", "", "");
      if(Ens_tag=="cB211b.072.48") { Corr.Tmin=20;	Corr.Tmax=35;  }
      else {	Corr.Tmin=30; 	Corr.Tmax=40;     }
      distr_t DM= Corr.Fit_distr(DM_distr);
      cout<<"(dmu+dmd)/2: "<<DM.ave()<<" "<<DM.err()<<endl;
      //#################################################



      //#################################################
      //compute  ud mass difference from MK^+ -MK^0    dmud = (dmd-dmu)/2
      distr_t_list dmud_distr = 4.0*M_PI*alpha*(qu*qu -qd*qd)*mass_shift( C_TM_ls_QED[P5P5][SELF_U], C_TM_ls_QED[P5P5][ISO_Q], ""); //self
      dmud_distr = dmud_distr + 4.0*M_PI*alpha*(qu*qs - qd*qs)*mass_shift( C_TM_ls_QED[P5P5][EXCH], C_TM_ls_QED[P5P5][ISO_Q], ""); //exchange
      dmud_distr = dmud_distr + (dmcr_u -dmcr_d)*mass_shift( C_TM_ls_SIB[P5P5][MCR_U], C_TM_ls_SIB[P5P5][ISO_S], ""); //critical mass
      dmud_distr = aDMK - FSE_charged_kaon - dmud_distr;
      dmud_distr= -0.5*dmud_distr/mass_shift( C_TM_ls_SIB[P5P5][M_U], C_TM_ls_SIB[P5P5][ISO_S], "");   // (dmd-dmu)/2
      if(Ens_tag=="cB211b.072.48") { Corr.Tmin=20;Corr.Tmax=35;}
      else {Corr.Tmin=27;Corr.Tmax=40;}
      Print_To_File( {}, {dmud_distr.ave(), dmud_distr.err()}, "../data/tau_LIBE/"+Ens_tag+"/REN/dmud.eff", "", "");
      distr_t dmud = Corr.Fit_distr(dmud_distr);
      distr_t dmu= (DM-dmud);
      distr_t dmd= (DM+dmud);
      cout<<"(md-mu)/2:    "<<dmud.ave()<<" "<<dmud.err()<<endl;
      //#################################################


      //#################################################
      //compute  dms mass difference from MK^+
      distr_t_list dms_distr = 4.0*M_PI*alpha*(  qu*qu*mass_shift( C_TM_ls_QED[P5P5][SELF_U], C_TM_ls_QED[P5P5][ISO_Q],"") + qs*qs*mass_shift(C_TM_ls_QED[P5P5][SELF_D],C_TM_ls_QED[P5P5][ISO_Q],"")); //self
      dms_distr = dms_distr + 4.0*M_PI*alpha*qu*qs*mass_shift( C_TM_ls_QED[P5P5][EXCH], C_TM_ls_QED[P5P5][ISO_Q], ""); //exch
      dms_distr = dms_distr + dmcr_u*mass_shift( C_TM_ls_SIB[P5P5][MCR_U], C_TM_ls_SIB[P5P5][ISO_S],"") + dmcr_s*mass_shift( C_TM_ls_SIB[P5P5][MCR_D], C_TM_ls_SIB[P5P5][ISO_S], ""); //critical mass
      dms_distr = dms_distr + dmu*mass_shift( C_TM_ls_SIB[P5P5][M_U], C_TM_ls_SIB[P5P5][ISO_S], ""); //mass shift of u-quark
      dms_distr = Delta_MK_iso -FSE_charged_kaon -dms_distr;
      dms_distr = dms_distr/mass_shift( C_TM_ls_SIB[P5P5][M_D], C_TM_ls_SIB[P5P5][ISO_S],"");
      if(Ens_tag=="cB211b.072.48") { Corr.Tmin=20;Corr.Tmax=35;}
      else {Corr.Tmin=31;Corr.Tmax=40;}
      Print_To_File( {}, {dms_distr.ave(), dms_distr.err()}, "../data/tau_LIBE/"+Ens_tag+"/REN/dms.eff", "", "");
      distr_t dms = Corr.Fit_distr(dms_distr);
      cout<<"dms:         "<<dms.ave()<<" "<<dms.err()<<endl;
      //#################################################


      cout<<"###   Counterterms for ensemble "<<Ens_tag<<" computed!  ###"<<endl;


      
      FLAV FLAV_U(qu, dmu, dmcr_u);
      FLAV FLAV_D(qd, dmd, dmcr_d);
      FLAV FLAV_S(qs, dms, dmcr_s);

      distr_t id= Get_id_jack_distr(Njacks);
      
      

      //determine VkVk and AkAk ls correlators

      distr_t_list dC_test_OS= -1.0*2*Zv*Zv*(qu*qu*dmu + qd*qd*dmd)*C_OS_ll_SIB[VkVk][M_D]  ;

      Print_To_File({ }, { dC_test_OS.ave(), dC_test_OS.err()} , "../data/tau_LIBE/"+Ens_tag+"/corr/SIB_test", "", "");


      distr_t_list dC_VkVk_TM = Get_LIBE_correlator( C_TM_ls_QED[VkVk], C_TM_ls_SIB[VkVk], FLAV_U, FLAV_S , Za*Za, "../data/tau_LIBE/"+Ens_tag+"/corr/dC_VkVk_TM");
      distr_t_list dC_AkAk_TM = Get_LIBE_correlator( C_TM_ls_QED[AkAk], C_TM_ls_SIB[AkAk], FLAV_U, FLAV_S , Zv*Zv,  "../data/tau_LIBE/"+Ens_tag+"/corr/dC_AkAk_TM");
      distr_t_list dC_VkVk_OS = Get_LIBE_correlator( C_OS_ls_QED[VkVk], C_OS_ls_SIB[VkVk], FLAV_U, FLAV_S , Zv*Zv, "../data/tau_LIBE/"+Ens_tag+"/corr/dC_VkVk_OS");
      distr_t_list dC_AkAk_OS = Get_LIBE_correlator( C_OS_ls_QED[AkAk], C_OS_ls_SIB[AkAk], FLAV_U, FLAV_S , Za*Za, "../data/tau_LIBE/"+Ens_tag+"/corr/dC_AkAk_OS");

     
      distr_t_list dC_V0V0_TM = Get_LIBE_correlator( C_TM_ls_QED[V0V0], C_TM_ls_SIB[V0V0], FLAV_U, FLAV_S , Za*Za, "../data/tau_LIBE/"+Ens_tag+"/corr/dC_V0V0_TM");
      distr_t_list dC_A0A0_TM = Get_LIBE_correlator( C_TM_ls_QED[A0A0], C_TM_ls_SIB[A0A0], FLAV_U, FLAV_S , Zv*Zv,  "../data/tau_LIBE/"+Ens_tag+"/corr/dC_A0A0_TM");
      distr_t_list dC_V0V0_OS = Get_LIBE_correlator( C_OS_ls_QED[V0V0], C_OS_ls_SIB[V0V0], FLAV_U, FLAV_S , Zv*Zv, "../data/tau_LIBE/"+Ens_tag+"/corr/dC_V0V0_OS");
      distr_t_list dC_A0A0_OS = Get_LIBE_correlator( C_OS_ls_QED[A0A0], C_OS_ls_SIB[A0A0], FLAV_U, FLAV_S , Za*Za, "../data/tau_LIBE/"+Ens_tag+"/corr/dC_A0A0_OS");


      //get delta correlator to be used for Z factors

      //distr_t_list dC_P5P5_TM_ss =  Get_LIBE_correlator( C_TM_ss_QED[P5P5], C_TM_ss_SIB[P5P5], FLAV_S, FLAV_S , id, "../data/tau_LIBE/"+Ens_tag+"/corr/dC_P5P5_TM_ss");
      //distr_t_list dC_A0P5_TM_ss =  Get_LIBE_correlator( C_TM_ss_QED[A0P5], C_TM_ss_SIB[A0P5], FLAV_S, FLAV_S , id, "../data/tau_LIBE/"+Ens_tag+"/corr/dC_A0P5_TM_ss");
      //distr_t_list dC_P5P5_OS_ss =  Get_LIBE_correlator( C_OS_ss_QED[P5P5], C_OS_ss_SIB[P5P5], FLAV_S, FLAV_S , id, "../data/tau_LIBE/"+Ens_tag+"/corr/dC_P5P5_OS_ss");
      //distr_t_list dC_A0P5_OS_ss =  Get_LIBE_correlator( C_OS_ss_QED[A0P5], C_OS_ss_SIB[A0P5], FLAV_S, FLAV_S , id, "../data/tau_LIBE/"+Ens_tag+"/corr/dC_A0P5_OS_ss");

      //distr_t_list dt_dC_A0P5_TM_ss = distr_t_list::derivative(dC_A0P5_TM_ss, 0);
      //distr_t_list dt_A0P5_TM_ss = distr_t_list::derivative(C_TM_ss_QED[A0P5][ISO_Q], 0);


      // get dZV and dZA

      //distr_t_list RV_QCD_QED = 2*ams*dC_P5P5_TM_ss/dt_A0P5_TM_ss  - 2*ams*C_TM_ss_QED[P5P5][ISO_Q]*dt_dC_A0P5_TM_ss/(dt_A0P5_TM_ss*dt_A0P5_TM_ss) + 2*(DM-dmud)*C_TM_ss_QED[P5P5][ISO_Q]/dt_A0P5_TM_ss ;
   


      //DO HLT

      //set kernel function in transverse and longitudinal channel (so far only transverse needed)

      const auto K1 = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

	PrecFloat X;
	PrecFloat X_ave = E/(m_tau*a_distr.ave());
	PrecFloat XX= E/(m_tau*a_distr.ave());
	if( X_ave < E0) return 0.0;	
	//if(X_ave > 1e4) return 0.0;
	if(ijack==-1) {
	  X=X_ave;
	}
	else X= E/(m_tau*a_distr.distr[ijack]);
	PrecFloat  sm_theta= 1/(1+ exp(-(1-X)/s));
	return (1 + 2*pow(X,2))*(1/(XX))*pow(( 1 -pow(X,2)),2)*sm_theta;
      };


      const auto K0 = [&a_distr](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

	PrecFloat X;
	PrecFloat X_ave = E/(m_tau*a_distr.ave());
	PrecFloat XX= E/(m_tau*a_distr.ave());
	if( X_ave < E0) return 0.0;	
	//if(X_ave > 1e4) return 0.0;
	if(ijack==-1) {
	  X=X_ave;
	}
	else X= E/(m_tau*a_distr.distr[ijack]);
	PrecFloat  sm_theta= 1/(1+ exp(-(1-X)/s));
	return (1/XX)*pow(( 1 -pow(X,2)),2)*sm_theta*sm_theta;
      };

      
  
     

      //compute covariance matrix
      Vfloat cov_Vk_TM, cov_Vk_OS, cov_Ak_TM, cov_Ak_OS;
      Vfloat corr_Vk_TM, corr_Vk_OS, corr_Ak_TM, corr_Ak_OS;
      
      Vfloat cov_V0_TM, cov_V0_OS, cov_A0_TM, cov_A0_OS;
      Vfloat corr_V0_TM, corr_V0_OS, corr_A0_TM, corr_A0_OS;
      
      Vfloat TT, RR;
      for(int tt=0;tt<Corr.Nt;tt++)
	for(int rr=0;rr<Corr.Nt;rr++) {
	  TT.push_back(tt);
	  RR.push_back(rr);
	  //cov
	  cov_Vk_TM.push_back( dC_VkVk_TM.distr_list[tt]%dC_VkVk_TM.distr_list[rr] );
	  cov_Vk_OS.push_back( dC_VkVk_OS.distr_list[tt]%dC_VkVk_OS.distr_list[rr] );
	  cov_Ak_TM.push_back( dC_AkAk_TM.distr_list[tt]%dC_AkAk_TM.distr_list[rr] );
	  cov_Ak_OS.push_back( dC_AkAk_OS.distr_list[tt]%dC_AkAk_OS.distr_list[rr] );
	  //corr
	  corr_Vk_TM.push_back( cov_Vk_TM[rr+ tt*Corr.Nt]/(dC_VkVk_TM.err(tt)*dC_VkVk_TM.err(rr)));
	  corr_Vk_OS.push_back( cov_Vk_OS[rr+ tt*Corr.Nt]/(dC_VkVk_OS.err(tt)*dC_VkVk_OS.err(rr))); 
	  corr_Ak_TM.push_back( cov_Ak_TM[rr+ tt*Corr.Nt]/(dC_AkAk_TM.err(tt)*dC_AkAk_TM.err(rr)));
	  corr_Ak_OS.push_back( cov_Ak_OS[rr+ tt*Corr.Nt]/(dC_AkAk_OS.err(tt)*dC_AkAk_OS.err(rr)));

	  //cov
	  cov_V0_TM.push_back( dC_V0V0_TM.distr_list[tt]%dC_V0V0_TM.distr_list[rr] );
	  cov_V0_OS.push_back( dC_V0V0_OS.distr_list[tt]%dC_V0V0_OS.distr_list[rr] );
	  cov_A0_TM.push_back( dC_A0A0_TM.distr_list[tt]%dC_A0A0_TM.distr_list[rr] );
	  cov_A0_OS.push_back( dC_A0A0_OS.distr_list[tt]%dC_A0A0_OS.distr_list[rr] );
	  //corr
	  corr_V0_TM.push_back( cov_V0_TM[rr+ tt*Corr.Nt]/(dC_V0V0_TM.err(tt)*dC_V0V0_TM.err(rr)));
	  corr_V0_OS.push_back( cov_V0_OS[rr+ tt*Corr.Nt]/(dC_V0V0_OS.err(tt)*dC_V0V0_OS.err(rr))); 
	  corr_A0_TM.push_back( cov_A0_TM[rr+ tt*Corr.Nt]/(dC_A0A0_TM.err(tt)*dC_A0A0_TM.err(rr)));
	  corr_A0_OS.push_back( cov_A0_OS[rr+ tt*Corr.Nt]/(dC_A0A0_OS.err(tt)*dC_A0A0_OS.err(rr)));


	  
	}
      cout<<"Covariance matrix computed"<<endl;

      //print
      boost::filesystem::create_directory("../data/tau_LIBE/"+Ens_tag+"/covariance");
      Print_To_File({}, {TT,RR, cov_Vk_TM, corr_Vk_TM}, "../data/tau_LIBE/"+Ens_tag+"/covariance/VkVk_TM.cov","","");
      Print_To_File({}, {TT,RR, cov_Vk_OS, corr_Vk_OS}, "../data/tau_LIBE/"+Ens_tag+"/covariance/VkVk_OS.cov","","");
      Print_To_File({}, {TT,RR, cov_Ak_TM, corr_Ak_TM}, "../data/tau_LIBE/"+Ens_tag+"/covariance/AkAk_TM.cov","","");
      Print_To_File({}, {TT,RR, cov_Ak_OS, corr_Ak_OS}, "../data/tau_LIBE/"+Ens_tag+"/covariance/AkAk_OS.cov","","");
      
      Print_To_File({}, {TT,RR, cov_V0_TM, corr_V0_TM}, "../data/tau_LIBE/"+Ens_tag+"/covariance/V0V0_TM.cov","","");
      Print_To_File({}, {TT,RR, cov_V0_OS, corr_V0_OS}, "../data/tau_LIBE/"+Ens_tag+"/covariance/V0V0_OS.cov","","");
      Print_To_File({}, {TT,RR, cov_A0_TM, corr_A0_TM}, "../data/tau_LIBE/"+Ens_tag+"/covariance/A0A0_TM.cov","","");
      Print_To_File({}, {TT,RR, cov_A0_OS, corr_A0_OS}, "../data/tau_LIBE/"+Ens_tag+"/covariance/A0A0_OS.cov","","");


      //######  HLT PARAMETERS #######
      int tmax=(int)( 2.5/(a_distr.ave()/fmTGeV));
      double mult_A = 7e4;
      double mult_V = 5e4;
      double mult_A0= 2e4;
      double mult_V0 = 2e4;
      double E0 = 0.9*MK_plus;
      double alpha= 3.99;
      double sigma=0.02;
      double Emax= 4.0;
      int prec=128; // ~quadruple prec
      int Is_Emax_Finite=1;
      double syst;
      double l;
      distr_t resc_GeV = C_V*GAMMA_FACT*Get_id_jack_distr(Njacks)/(pow(a_distr.ave(),3));
      double Ag_ov_A0_target=3e-3;
      //#############################
      
      distr_t dR_Ak_TM = Get_Laplace_transfo(  0.0, sigma, E0*a_distr.ave(),  T, tmax, prec, "KT", K1, dC_AkAk_TM, syst, mult_A, l, "TANT", "TM", "dAk_"+Ens_tag, Ag_ov_A0_target,0, resc_GeV, 0.0, "LIBE_tau_decay", cov_Ak_TM, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, alpha);
      dR_Ak_TM = dR_Ak_TM.ave() + (dR_Ak_TM-dR_Ak_TM.ave())*sqrt( 1.0 + pow(syst/dR_Ak_TM.err(),2));
      
      
      distr_t dR_Ak_OS = Get_Laplace_transfo(  0.0, sigma, E0*a_distr.ave(),  T, tmax, prec, "KT", K1, dC_AkAk_OS, syst, mult_A, l, "TANT", "OS", "dAk_"+Ens_tag, Ag_ov_A0_target,0, resc_GeV, 0.0, "LIBE_tau_decay", cov_Ak_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, alpha);
      dR_Ak_OS = dR_Ak_OS.ave() + (dR_Ak_OS-dR_Ak_OS.ave())*sqrt( 1.0 + pow(syst/dR_Ak_OS.err(),2));

      distr_t dR_Vk_TM = Get_Laplace_transfo(  0.0, sigma, E0*a_distr.ave(),  T, tmax, prec, "KT", K1, dC_VkVk_TM, syst, mult_V, l, "TANT", "TM", "dVk_"+Ens_tag, Ag_ov_A0_target,0, resc_GeV, 0.0, "LIBE_tau_decay", cov_Vk_TM, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, alpha);
      dR_Vk_TM = dR_Vk_TM.ave() + (dR_Vk_TM-dR_Vk_TM.ave())*sqrt( 1.0 + pow(syst/dR_Vk_TM.err(),2));
      
      distr_t dR_Vk_OS = Get_Laplace_transfo(  0.0, sigma, E0*a_distr.ave(),  T, tmax, prec, "KT", K1, dC_VkVk_OS, syst, mult_V, l, "TANT", "OS", "dVk_"+Ens_tag, Ag_ov_A0_target,0, resc_GeV, 0.0, "LIBE_tau_decay", cov_Vk_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, alpha);
      dR_Vk_OS = dR_Vk_OS.ave() + (dR_Vk_OS-dR_Vk_OS.ave())*sqrt( 1.0 + pow(syst/dR_Vk_OS.err(),2));


      //############################

      distr_t dR_A0_TM = Get_Laplace_transfo(  0.0, sigma, E0*a_distr.ave(),  T, tmax, prec, "KL", K0, dC_A0A0_TM, syst, mult_A0, l, "TANT", "TM", "dA0_"+Ens_tag, Ag_ov_A0_target,0, resc_GeV, 0.0, "LIBE_tau_decay", cov_A0_TM, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, alpha);
      dR_A0_TM = dR_A0_TM.ave() + (dR_A0_TM-dR_A0_TM.ave())*sqrt( 1.0 + pow(syst/dR_A0_TM.err(),2));
      
      distr_t dR_A0_OS = Get_Laplace_transfo(  0.0, sigma, E0*a_distr.ave(),  T, tmax, prec, "KL", K0, dC_A0A0_OS, syst, mult_A0, l, "TANT", "OS", "dA0_"+Ens_tag, Ag_ov_A0_target,0, resc_GeV, 0.0, "LIBE_tau_decay", cov_A0_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, alpha);
      dR_A0_OS = dR_A0_OS.ave() + (dR_A0_OS-dR_A0_OS.ave())*sqrt( 1.0 + pow(syst/dR_A0_OS.err(),2));

      distr_t dR_V0_TM = Get_Laplace_transfo(  0.0, sigma, E0*a_distr.ave(),  T, tmax, prec, "KL", K0, dC_V0V0_TM, syst, mult_V0, l, "TANT", "TM", "dV0_"+Ens_tag, Ag_ov_A0_target,0, resc_GeV, 0.0, "LIBE_tau_decay", cov_V0_TM, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, alpha);
      dR_V0_TM = dR_V0_TM.ave() + (dR_V0_TM-dR_V0_TM.ave())*sqrt( 1.0 + pow(syst/dR_V0_TM.err(),2));
      
      distr_t dR_V0_OS = Get_Laplace_transfo(  0.0, sigma, E0*a_distr.ave(),  T, tmax, prec, "KL", K0, dC_V0V0_OS, syst, mult_V0, l, "TANT", "OS", "dV0_"+Ens_tag, Ag_ov_A0_target,0, resc_GeV, 0.0, "LIBE_tau_decay", cov_V0_OS, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, alpha);
      dR_V0_OS = dR_V0_OS.ave() + (dR_V0_OS-dR_V0_OS.ave())*sqrt( 1.0 + pow(syst/dR_V0_OS.err(),2));


      //########  GET FACTORIZABLE CORRECTIONS TO RENORMALIZATION CONSTANTS  #####################

      
      





      //##########################################################################################

      


      //##########################
      cout<<"#####  HLT SUMMARY #####" <<endl;
      cout<<"dR(Ak,TM): "<<dR_Ak_TM.ave()<<" +- "<<dR_Ak_TM.err()<<endl;
      cout<<"dR(Ak,OS): "<<dR_Ak_OS.ave()<<" +- "<<dR_Ak_OS.err()<<endl;
      cout<<"dR(Vk,TM): "<<dR_Vk_TM.ave()<<" +- "<<dR_Vk_TM.err()<<endl;
      cout<<"dR(Vk,OS): "<<dR_Vk_OS.ave()<<" +- "<<dR_Vk_OS.err()<<endl;

      cout<<"dR(Ak,TM): "<<dR_A0_TM.ave()<<" +- "<<dR_A0_TM.err()<<endl;
      cout<<"dR(Ak,OS): "<<dR_A0_OS.ave()<<" +- "<<dR_A0_OS.err()<<endl;
      cout<<"dR(Vk,TM): "<<dR_V0_TM.ave()<<" +- "<<dR_V0_TM.err()<<endl;
      cout<<"dR(Vk,OS): "<<dR_V0_OS.ave()<<" +- "<<dR_V0_OS.err()<<endl;
      
      cout<<"########################"<<endl;

      Ak_distr_TM.distr_list.push_back( dR_Ak_TM);
      Vk_distr_TM.distr_list.push_back( dR_Vk_TM);
      A0_distr_TM.distr_list.push_back( dR_A0_TM);
      V0_distr_TM.distr_list.push_back( dR_V0_TM);

      Ak_distr_OS.distr_list.push_back( dR_Ak_OS);
      Vk_distr_OS.distr_list.push_back( dR_Vk_OS);
      A0_distr_OS.distr_list.push_back( dR_A0_OS);
      V0_distr_OS.distr_list.push_back( dR_V0_OS);

      

      a_distr_list.distr_list.push_back(a_distr/fmTGeV);

      s=sigma;
	
      
  }

  
  boost::filesystem::create_directory("../data/tau_LIBE/cont");
  Print_To_File({}, {a_distr_list.ave(), Ak_distr_TM.ave(), Ak_distr_TM.err()}, "../data/tau_LIBE/cont/Ak_TM_s_"+to_string_with_precision(s,3), "","");
  Print_To_File({}, {a_distr_list.ave(), Vk_distr_TM.ave(), Vk_distr_TM.err()}, "../data/tau_LIBE/cont/Vk_TM_s_"+to_string_with_precision(s,3), "","");
  Print_To_File({}, {a_distr_list.ave(), A0_distr_TM.ave(), A0_distr_TM.err()}, "../data/tau_LIBE/cont/A0_TM_s_"+to_string_with_precision(s,3), "","");
  Print_To_File({}, {a_distr_list.ave(), V0_distr_TM.ave(), V0_distr_TM.err()}, "../data/tau_LIBE/cont/V0_TM_s_"+to_string_with_precision(s,3), "","");
  
  Print_To_File({}, {a_distr_list.ave(), Ak_distr_OS.ave(), Ak_distr_OS.err()}, "../data/tau_LIBE/cont/Ak_OS_s_"+to_string_with_precision(s,3), "","");
  Print_To_File({}, {a_distr_list.ave(), Vk_distr_OS.ave(), Vk_distr_OS.err()}, "../data/tau_LIBE/cont/Vk_OS_s_"+to_string_with_precision(s,3), "","");
  Print_To_File({}, {a_distr_list.ave(), A0_distr_OS.ave(), A0_distr_OS.err()}, "../data/tau_LIBE/cont/A0_OS_s_"+to_string_with_precision(s,3), "","");
  Print_To_File({}, {a_distr_list.ave(), V0_distr_OS.ave(), V0_distr_OS.err()}, "../data/tau_LIBE/cont/V0_OS_s_"+to_string_with_precision(s,3), "","");

  return; 

}
