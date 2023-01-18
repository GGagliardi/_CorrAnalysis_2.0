#include "../include/vph_Nissa_3d.h"

using namespace std;


const double M2PiPhys=pow(0.135,2);
const double alpha= 1/137.04;
const double e2 = alpha*4.0*M_PI;
const int Nboots= 800;
const bool UseJack=1;
const int nboots=150;
const int Njacks=30;
const double qu = 2.0/3.0; //electric charge of u-type quark
const double qd = -1.0/3.0; //electric charge of d-type quark
const string Meson="Ds";
double L_QCD= 0.3; //300 MeV
int n_xg=4;
int n_xg_rev=3;
Vfloat virt_list;
bool verbose_lev=0;
bool P5_ON_SOURCE=true;
bool Is_rep=false;
Vfloat sigmas({0.6, 0.5, 0.4, 0.3, 0.2}); //sigma in GeV
int prec=256;
const string MODE_FF="TANT";
bool CONS_EM_CURR=false;
const bool Skip_spectral_reconstruction=false;
const bool Reconstruct_axial_part=false;
const bool Reconstruct_vector_part=true;
const double Mjpsi= 3.0969;
const double Mphi= 1.019461;

void Get_virt_list() {

  int Nvirts=40;
  for(int nv=0;nv<Nvirts;nv++) virt_list.push_back( nv/(Nvirts-1.0));

  return;
}


void Integrate_over_photon_insertion(const distr_t_list &W, vector<distr_t_list> &H, double Eg, int t_weak, double MP, int TO) {

  int Nt= W.size();
  if(verbose_lev==1) {
    cout<<"Call to: Integrate_over_photon_insertion: "<<endl;
    cout<<"T: "<<Nt<<", t_weak: "<<t_weak<<", Eg: "<<Eg<<endl;
    cout<<"######################"<<endl;
  }

  auto HeavyTheta=[](const int x)
  {
    return ((x>=0)+(x>0))/2.0;
  };


  for(int nv=0;nv<(signed)virt_list.size();nv++) {

   
    
    double off= virt_list[nv];
    double Eg_n= sqrt( Eg*Eg + pow(MP*off,2));
    bool BW=(t_weak>(Nt/2));

    int tmin=0;
    int tmax=Nt;

    if(TO==0) {tmin=0; tmax=Nt;}
    else if(TO==1) {tmin=0; tmax=t_weak+1;}
    else if(TO==2) {tmin=t_weak+1; tmax=Nt;}
    else crash("Time ordering: "+to_string(TO)+" not defined!");

    for(int tcut=0;tcut<Nt;tcut++) {

      distr_t H_virt(UseJack, UseJack?Njacks:Nboots); //constructor sets to zero

      H_virt= 0.0*H_virt;

      int tmax_cut= min(tmax,tcut);
      
      for(int t=tmin;t<tmax_cut;t++) {

	int ty= t;

	if(tmax_cut==0) crash("tmax_cut should not be zero inside the t-loop");
      
	const double f1=                                                                                                                                                        
	  (BW==true)?exp(-(Nt/2+ty)*Eg_n):exp(-(Nt/2-ty)*Eg_n);

	const double f2=
	  (BW==true)?exp((Nt/2-ty)*Eg_n):exp(-((3*Nt/2)-ty)*Eg_n);

	const double h1=HeavyTheta((Nt/2)-ty);
	const double h2=HeavyTheta(ty-(Nt/2));

	H_virt= H_virt + W.distr_list[t]*(h1*f1+h2*f2);

      }

      H[nv].distr_list.push_back(exp(Eg_n*abs(Nt/2 - t_weak))*H_virt);
    }
  }


  return;

}


void Integrate_over_photon_insertion_w_subtraction(const distr_t_list &W, vector<distr_t_list> &H, double Eg, int t_weak, double MP, int ixg, int Tmin_mass, int Tmax_mass,  string out_path, string obs) {

  int Nt= W.size();
  if(verbose_lev==1) {
    cout<<"Call to: Integrate_over_photon_insertion and perform subtraction of first exponential: "<<endl;
    cout<<"T: "<<Nt<<", t_weak: "<<t_weak<<", Eg: "<<Eg<<endl;
    cout<<"######################"<<endl;
  }



  for(int nv=0;nv<(signed)virt_list.size();nv++) {

   
    
    double off= virt_list[nv];
    double Eg_n= sqrt( Eg*Eg + pow(MP*off,2));
    bool BW=(t_weak>(Nt/2));
    if(BW) crash("Integrate_over_photon_insertion_w_subtraction assumes BW=false, but BW=true");

    int tmin=t_weak+1;
    int tmax=Nt;

    distr_t_list V_tcut(UseJack);

      
    for(int tcut=0;tcut<=Nt/2;tcut++) {

      distr_t H_virt(UseJack, UseJack?Njacks:Nboots); //constructor sets to zero

     
      int tmax_cut= min(tmax,tcut);
      
      for(int t=tmin;t<tmax_cut;t++) {



	if(tmax_cut==0) crash("tmax_cut should not be zero inside the t-loop");
      

	H_virt= H_virt + W.distr_list[t]*exp(Eg_n*(t-t_weak));

      }

      V_tcut.distr_list.push_back(H_virt);
      
    }

    //remove leading exponential
    auto LOG= [](double m) { return log(m);};
    distr_t_list eff_mass_distr(UseJack);
    for(int tcut=tmin;tcut<Nt/2;tcut++) {
      eff_mass_distr.distr_list.push_back(Eg_n - distr_t::f_of_distr(LOG, (W.distr_list[tcut]/W.distr_list[tcut+1])));
    }
    //fit effective mass
    CorrAnalysis Corr_tcut(UseJack, Njacks,Nboots);
    Corr_tcut.Nt = Nt/2-tmin;
    Corr_tcut.Reflection_sign=0;
    Corr_tcut.Perform_Nt_t_average=0;
    //time interval fot effective mass fit
    Corr_tcut.Tmin=Tmin_mass;
    Corr_tcut.Tmax=Tmax_mass;
    distr_t eff_mass= Corr_tcut.Fit_distr(eff_mass_distr);
    cout<<"ixk: "<<nv<<" Eg-En: "<<eff_mass.ave()<<" +- "<<eff_mass.err()<<endl;
    distr_t_list overlap_distr(UseJack);
    auto EXP= [](double m, double t, double nt) { return exp(m*(t+1));};
    distr_t_list exp_eff_mass= distr_t_list::f_of_distr(EXP, eff_mass, Corr_tcut.Nt+1);
    for(int tcut=tmin;tcut<Nt/2;tcut++) overlap_distr.distr_list.push_back( V_tcut.distr_list[tcut]/exp_eff_mass.distr_list[tcut-tmin]);
    //time interval for effective residue fit
    if(ixg >= 1) {
      if(nv <=10) { Corr_tcut.Tmin=25; Corr_tcut.Tmax=35;}
      else if(nv==11) {Corr_tcut.Tmin=28; Corr_tcut.Tmax=35;}
      else if(nv==12) {Corr_tcut.Tmin=25; Corr_tcut.Tmax=35;}
      else if(nv==13) {Corr_tcut.Tmin=24; Corr_tcut.Tmax=35;}
      else if(nv==14) {Corr_tcut.Tmin=24; Corr_tcut.Tmax=33;}
      else if(nv==15) {Corr_tcut.Tmin=23; Corr_tcut.Tmax=30;}
      else if(nv==16) {Corr_tcut.Tmin=16; Corr_tcut.Tmax=30;}
      else if(nv==17) {Corr_tcut.Tmin=15; Corr_tcut.Tmax=30;}
      else if(nv==18) {Corr_tcut.Tmin=14; Corr_tcut.Tmax=30;}
      else if(nv >= 19) {Corr_tcut.Tmin=12; Corr_tcut.Tmax=30;}
    }
    
    distr_t eff_overlap= Corr_tcut.Fit_distr(overlap_distr);

  
    //print eff_mass_distr and overlap_distr
    Print_To_File({}, { eff_mass_distr.ave(), eff_mass_distr.err()}, out_path+"/eff_vector_mass_"+obs+"_ixg_"+to_string(ixg)+"_ixk_"+to_string(nv)+".dat", "", "");
    Print_To_File({}, { overlap_distr.ave(),  overlap_distr.err()}, out_path+"/eff_vector_residue_"+obs+"_ixg_"+to_string(ixg)+"_ixk_"+to_string(nv)+".dat", "", "");

    for(int tcut=0;tcut<Nt/2+1;tcut++) {
      distr_t sub_exp(UseJack, UseJack?Njacks:Nboots);
      if(tcut >= tmin) {
	for(int ijack=0;ijack<Njacks;ijack++) sub_exp.distr[ijack] =  ((eff_mass.distr[ijack] > 0)?eff_overlap.distr[ijack]*exp_eff_mass.distr_list[tcut-tmin].distr[ijack]:0.0);
      }
      H[nv].distr_list.push_back( V_tcut.distr_list[tcut] - sub_exp);
    }

  }

  cout<<"subtraction performed"<<endl;
  return;

}



void Get_radiative_form_factors_3d() {


  Get_virt_list();
  
  int rank, size;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /*
  Vfloat beta_List({0.0, 0.0, 0.0, 0.99, 0.99, 0.99, -0.99, -0.99});
  vector<bool> Integrate_Up_To_Emax_List({0,0,0,0,0,0,0,0});
  Vfloat Emax_List({10,10,10,10,10,10,10,10});
  vector<bool> Perform_theta_average_List({1,1,1,1,1,1,1,1});
  vector<string> SM_TYPE_List({"FF_Sinh", "FF_Sinh", "FF_Sinh", "FF_Sinh", "FF_Sinh", "FF_Sinh", "FF_Sinh", "FF_Sinh"});
  Vfloat E0_List({0.6,0.8,0.9,0.6,0.8,0.9,0.8,0.9});
  */

  Vfloat beta_List({0.0, 0.0});
  vector<bool> Integrate_Up_To_Emax_List({0,0});
  Vfloat Emax_List({10,10});
  vector<bool> Perform_theta_average_List({1,1});
  vector<string> SM_TYPE_List({"FF_Sinh", "FF_Sinh"});
  vector<bool> CONS_EM_CURR_LIST({true, false});
  Vfloat E0_List({0.9,0.8});
  
  int N= beta_List.size();

  
  if(N%size != 0) crash("MPI called with -np= "+to_string(size)+". np does not divide vector size N="+to_string(N));
  
  for(int i=rank*N/size;i<(rank+1)*N/size;i++) {
    cout<<"################# DETERMINATION OF VIRTUAL-RADIATIVE FF USING HLT & EXP-SUB METHODS #################"<<endl;
    cout<<"Rank: "<<rank<<" pid: "<<getpid()<<" core id: "<<"("<<sched_getcpu()<<")"<<endl;
    cout<<"RECONSTRUCTION CALLED FOR:"<<endl;
    cout<<"alpha: "<<beta_List[i]<<", Use_Emax: "<<Integrate_Up_To_Emax_List[i]<<", SM_TYPE: "<<SM_TYPE_List[i]<<", theta average: "<<Perform_theta_average_List[i]<<endl;
    CONS_EM_CURR= CONS_EM_CURR_LIST[i];
    cout<<"electromagnetic current: "<<((CONS_EM_CURR)?"conserved":"local")<<endl;
    Compute_form_factors_Nissa_3d(beta_List[i], Integrate_Up_To_Emax_List[i], Emax_List[i], SM_TYPE_List[i], Perform_theta_average_List[i], E0_List[i]);
  }
  cout<<"##########################################"<<endl;

 
  return ;
}

void Compute_form_factors_Nissa_3d(double beta, bool Integrate_Up_To_Emax, double Emax, string SM_TYPE, bool Perform_theta_average, double E0_fact) {


  PrecFloat::setDefaultPrecision(prec);
  cout<<"max possible exponent: "<<PrecFloat::getEmax_max()<<endl;
  cout<<"current max exponent: "<<PrecFloat::getEmax()<<endl;


  string TAG_CURR="";
  if(CONS_EM_CURR==false) TAG_CURR="LOC_";

  double E0_fact_u=E0_fact;
  double E0_fact_d=E0_fact;

  int size_mu_nu= 4;
  string ph_type= Is_rep?"rph":"vph";

  string ph_type_mes=ph_type+"/"+Meson;
 
  //create directories
  boost::filesystem::create_directory("../data/ph_emission_3d");
  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type);
  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type+"/"+Meson);
  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type+"/"+Meson+"/C");
  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type+"/"+Meson+"/H");
  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type+"/"+Meson+"/mass");
  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type+"/"+Meson+"/covariance");
  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type+"/"+Meson+"/decay_const");
  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type+"/"+Meson+"/FF");
  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type+"/"+Meson+"/FF/continuum");
  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type+"/"+Meson+"/FF/per_kin");
  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type+"/"+Meson+"/FF_u");
  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type+"/"+Meson+"/FF_d");
   
  //axial
  vector<vector<vector<data_t>>> C_A_u_data(size_mu_nu), C_A_d_data(size_mu_nu);
  //vector
  vector<vector<vector<data_t>>> C_V_u_data(size_mu_nu), C_V_d_data(size_mu_nu);

  //axial rev theta
  vector<vector<vector<data_t>>> C_A_u_data_rev(size_mu_nu), C_A_d_data_rev(size_mu_nu);
  //vector rew theta
  vector<vector<vector<data_t>>> C_V_u_data_rev(size_mu_nu), C_V_d_data_rev(size_mu_nu);
  


  //2pts
  data_t data_2pts;
  data_t data_2pts_SM;
  
  for(int mu=0;mu<size_mu_nu;mu++) {

    C_A_u_data[mu].resize(size_mu_nu);
    C_A_d_data[mu].resize(size_mu_nu);
    C_V_u_data[mu].resize(size_mu_nu);
    C_V_d_data[mu].resize(size_mu_nu);

    C_A_u_data_rev[mu].resize(size_mu_nu);
    C_A_d_data_rev[mu].resize(size_mu_nu);
    C_V_u_data_rev[mu].resize(size_mu_nu);
    C_V_d_data_rev[mu].resize(size_mu_nu);
   
    for(int nu=0;nu<size_mu_nu;nu++) {

      C_A_u_data[mu][nu].resize(n_xg);
      C_A_d_data[mu][nu].resize(n_xg);
      C_V_u_data[mu][nu].resize(n_xg);
      C_V_d_data[mu][nu].resize(n_xg);

      C_A_u_data_rev[mu][nu].resize(n_xg_rev);
      C_A_d_data_rev[mu][nu].resize(n_xg_rev);
      C_V_u_data_rev[mu][nu].resize(n_xg_rev);
      C_V_d_data_rev[mu][nu].resize(n_xg_rev);
    
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

  data_2pts.Read("../new_vph_3d_gpu_data", "mes_contr_2pts", "P5P5", Sort_confs);
  //data_2pts_SM.Read("../new_vph_3d_gpu_data", "mes_contr_2pts_SM", "P5P5");


  //loop over mu and nu axial
  vector<pair<int,int>> mu_nu_pair_A({make_pair(1,1),make_pair(2,2),make_pair(3,3),make_pair(0,3),make_pair(3,0)});
  vector<pair<int,int>> mu_nu_pair_V({make_pair(1,2),make_pair(2,1)});
  vector<pair<int,int>> red_mu_nu_pair_V({make_pair(1,2)});
  vector<pair<int,int>> red_mu_nu_pair_A({make_pair(1,1), make_pair(3,3), make_pair(3,0), make_pair(0,3) });
  
  
  if(Is_rep) mu_nu_pair_A = { make_pair(1,1), make_pair(2,2) };
  
  for(int ixg=0;ixg<n_xg;ixg++) {


   
    //axial
    for(auto &pair_A : mu_nu_pair_A) {
      int mu=pair_A.first;
      int nu=pair_A.second;
      if(!P5_ON_SOURCE) {

	//axial
	//u
	C_A_u_data[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"C_mu_"+to_string(mu)+"_u_ixg_"+to_string(ixg), "S0A"+to_string(nu), Sort_confs);
	//d
	C_A_d_data[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"C_mu_"+to_string(mu)+"_d_ixg_"+to_string(ixg), "S0A"+to_string(nu), Sort_confs);
      }

      else {

	string Tag_contr="S0P5";
	if(CONS_EM_CURR==false) Tag_contr="V"+to_string(mu)+"P5";
	//axial
	//u
	C_A_u_data[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"C_u_A_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
	//d
	C_A_d_data[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"C_d_A_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);

      }
    }

    //vector
    for(auto &pair_V : mu_nu_pair_V) {
      int mu=pair_V.first;
      int nu=pair_V.second;
      if(!P5_ON_SOURCE) {
	
	//vector
	//u
	C_V_u_data[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"C_mu_"+to_string(mu)+"_u_ixg_"+to_string(ixg), "S0V"+to_string(nu), Sort_confs);
	//d
	C_V_d_data[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"C_mu_"+to_string(mu)+"_d_ixg_"+to_string(ixg), "S0V"+to_string(nu), Sort_confs);
      }

      else {

	string Tag_contr="S0P5";
	if(CONS_EM_CURR==false) Tag_contr="V"+to_string(mu)+"P5";
	//vector
	//u
	C_V_u_data[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"C_u_V_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
	//d
	C_V_d_data[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"C_d_V_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
	
      }
    }
  }



  for(int ixg=0;ixg<n_xg_rev;ixg++) {


   
    //axial
    for(auto &pair_A : mu_nu_pair_A) {
      int mu=pair_A.first;
      int nu=pair_A.second;
      if(!P5_ON_SOURCE) {

	//axial
	//u
	C_A_u_data_rev[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"REV_C_mu_"+to_string(mu)+"_u_ixg_"+to_string(ixg), "S0A"+to_string(nu), Sort_confs);
	//d
	C_A_d_data_rev[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"REV_C_mu_"+to_string(mu)+"_d_ixg_"+to_string(ixg), "S0A"+to_string(nu), Sort_confs);
      }

      else {

	string Tag_contr="S0P5";
	if(CONS_EM_CURR==false) Tag_contr="V"+to_string(mu)+"P5";
	//axial
	//u
	C_A_u_data_rev[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"REV_C_u_A_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
	//d
	C_A_d_data_rev[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"REV_C_d_A_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);

      }
    }

    //vector
    for(auto &pair_V : mu_nu_pair_V) {
      int mu=pair_V.first;
      int nu=pair_V.second;
      if(!P5_ON_SOURCE) {
	
	//vector
	//u
	C_V_u_data_rev[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"REV_C_mu_"+to_string(mu)+"_u_ixg_"+to_string(ixg), "S0V"+to_string(nu), Sort_confs);
	//d
	C_V_d_data_rev[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"REV_C_mu_"+to_string(mu)+"_d_ixg_"+to_string(ixg), "S0V"+to_string(nu), Sort_confs);
      }

      else {

	string Tag_contr="S0P5";
	if(CONS_EM_CURR==false) Tag_contr="V"+to_string(mu)+"P5";
	//vector
	//u
	C_V_u_data_rev[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"REV_C_u_V_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
	//d
	C_V_d_data_rev[mu][nu][ixg].Read("../new_vph_3d_gpu_data", TAG_CURR+"REV_C_d_V_nu_"+to_string(nu)+"_mu_"+to_string(mu)+"_ixg_"+to_string(ixg), Tag_contr, Sort_confs);
	
      }
    }
  }
  
  GaussianMersenne GM(652205123);
 



  //resample RCs
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack);

  double fmTGeV= 1.0/0.197327;

  LatticeInfo L_info_A, L_info_B, L_info_C, L_info_D;
  L_info_A.LatInfo_new_ens("cA211a.12.48");
  L_info_B.LatInfo_new_ens("cB211b.072.96");
  L_info_C.LatInfo_new_ens("cC211a.06.80");
  L_info_D.LatInfo_new_ens("cD211a.054.96");
  

  for(int ijack=0; ijack<Njacks;ijack++) {

    ZA_A.distr.push_back( L_info_A.Za_WI_strange + GM()*L_info_A.Za_WI_strange_err/sqrt( Njacks -1.0));
    ZV_A.distr.push_back( L_info_A.Zv_WI_strange + GM()*L_info_A.Zv_WI_strange_err/sqrt( Njacks -1.0));

    ZA_B.distr.push_back( L_info_B.Za_WI_strange + GM()*L_info_B.Za_WI_strange_err/sqrt( Njacks -1.0));
    ZV_B.distr.push_back( L_info_B.Zv_WI_strange + GM()*L_info_B.Zv_WI_strange_err/sqrt( Njacks -1.0));

    ZA_C.distr.push_back( L_info_C.Za_WI_strange + GM()*L_info_C.Za_WI_strange_err/sqrt( Njacks -1.0));
    ZV_C.distr.push_back( L_info_C.Zv_WI_strange + GM()*L_info_C.Zv_WI_strange_err/sqrt( Njacks -1.0));

    ZA_D.distr.push_back( L_info_D.Za_WI_strange + GM()*L_info_D.Za_WI_strange_err/sqrt( Njacks -1.0));
    ZV_D.distr.push_back( L_info_D.Zv_WI_strange + GM()*L_info_D.Zv_WI_strange_err/sqrt( Njacks -1.0));

    a_A.distr.push_back( L_info_A.a_from_afp*fmTGeV + GM()*L_info_A.a_from_afp_err*fmTGeV/sqrt(Njacks-1.0));
    a_B.distr.push_back( L_info_B.a_from_afp*fmTGeV + GM()*L_info_B.a_from_afp_err*fmTGeV/sqrt(Njacks-1.0));
    a_C.distr.push_back( L_info_C.a_from_afp*fmTGeV + GM()*L_info_C.a_from_afp_err*fmTGeV/sqrt(Njacks-1.0));
    a_D.distr.push_back( L_info_D.a_from_afp*fmTGeV + GM()*L_info_D.a_from_afp_err*fmTGeV/sqrt(Njacks-1.0));

  }


  int Nens= C_A_u_data[1][1][0].size;
  vector<string> Ens_tags= C_A_u_data[1][1][0].Tag;
  Vint Nts=C_A_u_data[1][1][0].nrows;


  //vector where to store smeared form factors
  vector<vector<vector<vector<vector<distr_t_list>>>>>  RE_HA_sm_u(Nens);
  vector<vector<vector<vector<vector<distr_t_list>>>>>  RE_HA_sm_d(Nens);
  vector<vector<vector<vector<vector<distr_t_list>>>>>  RE_HV_sm_u(Nens);
  vector<vector<vector<vector<vector<distr_t_list>>>>>  RE_HV_sm_d(Nens);
  vector<vector<vector<vector<vector<distr_t_list>>>>>  IM_HA_sm_u(Nens);
  vector<vector<vector<vector<vector<distr_t_list>>>>>  IM_HA_sm_d(Nens);
  vector<vector<vector<vector<vector<distr_t_list>>>>>  IM_HV_sm_u(Nens);
  vector<vector<vector<vector<vector<distr_t_list>>>>>  IM_HV_sm_d(Nens);

  //cascade resize
  for(int i=0;i<Nens;i++) {

    RE_HA_sm_u[i].resize(size_mu_nu);
    RE_HA_sm_d[i].resize(size_mu_nu);
    RE_HV_sm_u[i].resize(size_mu_nu);
    RE_HV_sm_d[i].resize(size_mu_nu);
    IM_HA_sm_u[i].resize(size_mu_nu);
    IM_HA_sm_d[i].resize(size_mu_nu);
    IM_HV_sm_u[i].resize(size_mu_nu);
    IM_HV_sm_d[i].resize(size_mu_nu);

    for(int mu=0;mu<size_mu_nu;mu++) {

      RE_HA_sm_u[i][mu].resize(size_mu_nu);
      RE_HA_sm_d[i][mu].resize(size_mu_nu);
      RE_HV_sm_u[i][mu].resize(size_mu_nu);
      RE_HV_sm_d[i][mu].resize(size_mu_nu);
      IM_HA_sm_u[i][mu].resize(size_mu_nu);
      IM_HA_sm_d[i][mu].resize(size_mu_nu);
      IM_HV_sm_u[i][mu].resize(size_mu_nu);
      IM_HV_sm_d[i][mu].resize(size_mu_nu);

      for(int nu=0;nu<size_mu_nu;nu++) {

	RE_HA_sm_u[i][mu][nu].resize(n_xg);
	RE_HA_sm_d[i][mu][nu].resize(n_xg);
	RE_HV_sm_u[i][mu][nu].resize(n_xg);
	RE_HV_sm_d[i][mu][nu].resize(n_xg);
	IM_HA_sm_u[i][mu][nu].resize(n_xg);
	IM_HA_sm_d[i][mu][nu].resize(n_xg);
	IM_HV_sm_u[i][mu][nu].resize(n_xg);
	IM_HV_sm_d[i][mu][nu].resize(n_xg);

	for(int ixg=0;ixg<n_xg;ixg++) {

	  for(int is=0;is<(signed)sigmas.size();is++) {

	    RE_HA_sm_u[i][mu][nu][ixg].emplace_back(UseJack, virt_list.size());
	    RE_HA_sm_d[i][mu][nu][ixg].emplace_back(UseJack, virt_list.size());
	    RE_HV_sm_u[i][mu][nu][ixg].emplace_back(UseJack, virt_list.size());
	    RE_HV_sm_d[i][mu][nu][ixg].emplace_back(UseJack, virt_list.size());
	    IM_HA_sm_u[i][mu][nu][ixg].emplace_back(UseJack, virt_list.size());
	    IM_HA_sm_d[i][mu][nu][ixg].emplace_back(UseJack, virt_list.size());
	    IM_HV_sm_u[i][mu][nu][ixg].emplace_back(UseJack, virt_list.size());
	    IM_HV_sm_d[i][mu][nu][ixg].emplace_back(UseJack, virt_list.size());

	  }
	}
      }
    }
  }

  
  
  int t_weak=22;




  //smeared kernel of the real part
  auto K_RE= [&SM_TYPE](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {


    if(SM_TYPE=="FF_Gauss") {
      PrecFloat x = (E-m)/(sqrt(PrecFloat(2))*s);
      return sqrt(2)*DawsonF(x)/s;
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
  auto K_IM = [&SM_TYPE](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

    if(SM_TYPE=="FF_Gauss") return precPi()*Get_exact_gauss(E, m, s, E0);

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

  for(int iens=0; iens<Nens;iens++) {

    boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]);
    boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/covariance/"+Ens_tags[iens]);
    boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]);
    boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]);
    boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/fit_results");
    boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/FF_u/"+Ens_tags[iens]);
    boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/FF_u/"+Ens_tags[iens]+"/fit_results");
    boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/FF_d/"+Ens_tags[iens]);
    boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/FF_d/"+Ens_tags[iens]+"/fit_results");
    boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/mass/"+Ens_tags[iens]);
    boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/decay_const/"+Ens_tags[iens]);


    cout<<"Analyzing ensemble: "<<Ens_tags[iens]<<endl;

    //Lattice info
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(Ens_tags[iens]);
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    CorrAnalysis Corr_boot(0, Njacks,Nboots, iens);
    Corr.Nt = Nts[iens];
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;
    Corr_boot.Nt= Nts[iens];
    Corr_boot.Reflection_sign=1;
    Corr_boot.Perform_Nt_t_average=1;


  


    //read theta values and loop over them
    Vfloat thetas, masses_u, masses_d;

    thetas= Read_From_File("../new_vph_3d_gpu_data/"+Ens_tags[iens]+"/pars_list.dat", 1 , 5);
    masses_u= Read_From_File("../new_vph_3d_gpu_data/"+Ens_tags[iens]+"/pars_list.dat", 3 , 5);
    masses_d= Read_From_File("../new_vph_3d_gpu_data/"+Ens_tags[iens]+"/pars_list.dat", 4 , 5);

    //read opposite theta values
    Vfloat thetas_rev;
    thetas_rev= Read_From_File("../new_vph_3d_gpu_data/"+Ens_tags[iens]+"/pars_list_rev.dat", 1 , 5);

    cout<<"pars_list.dat: Read!"<<endl;

    if((signed)thetas.size() != n_xg) crash("Number of rows in pars_list.dat does not match n_xg");
    if((signed)thetas_rev.size() != n_xg_rev) crash("Number of rows in pars_list_rev.dat does not match n_xg_rev");

    //load smeared 2-pt function

    //RCs
    distr_t Za, Zv, a_distr;
    if(data_2pts.Tag[iens].substr(1,1)=="A") { Za= ZA_A; Zv=ZV_A;a_distr=a_A;}
    else if(data_2pts.Tag[iens].substr(1,1)=="B") { Za= ZA_B; Zv=ZV_B;a_distr=a_B;}
    else if(data_2pts.Tag[iens].substr(1,1)=="C") { Za= ZA_C; Zv=ZV_C; a_distr=a_C;}
    else if(data_2pts.Tag[iens].substr(1,1)=="D") { Za= ZA_D; Zv=ZV_D;a_distr=a_D;}
    else crash("Ensemble: "+data_2pts.Tag[iens]+" not recognised");


   

    //read masses
    double mu= masses_u[0];
    double md= masses_d[0];


    cout<<"ZA: "<<Za.ave()<<" +- "<<Za.err()<<endl;
    cout<<"ZV: "<<Zv.ave()<<" +- "<<Zv.err()<<endl;
    cout<<"mu: "<<mu<<endl;
    cout<<"md: "<<md<<endl;
    

    //set time interval for eff_mass_fit
    if(data_2pts.Tag[iens].substr(1,1) =="A") {Corr.Tmin=19; Corr.Tmax=35;}
    else if(data_2pts.Tag[iens].substr(1,1) =="B") {Corr.Tmin=25; Corr.Tmax=40;}
    else if(data_2pts.Tag[iens].substr(1,1) == "C") {Corr.Tmin=33;Corr.Tmax=51;}
    else if(data_2pts.Tag[iens].substr(1,1) == "D") {Corr.Tmin=35;Corr.Tmax=54;}
    else crash("In fixing [Tmin, Tmax] for MP, Ensemble: "+data_2pts.Tag[iens]+" not recognized");

    distr_t_list pt2_distr= Corr.corr_t(data_2pts.col(0)[iens], "../data/ph_emission_3d/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt.dat");
    distr_t_list eff_mass = Corr.effective_mass_t(pt2_distr, "../data/ph_emission_3d/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass.dat");
    distr_t_list fp_distr= Corr.decay_constant_t( pow( mu+md,2)*pt2_distr, "../data/ph_emission_3d/"+ph_type_mes+"/"+"decay_const/"+data_2pts.Tag[iens]+"/decay_const.dat");
    distr_t MP= Corr.Fit_distr(eff_mass);
    distr_t FP= Corr.Fit_distr(fp_distr);

    cout<<"MP: "<<MP.ave()<<" +- "<<MP.err()<<endl;
    cout<<"FP: "<<FP.ave()<<" +- "<<FP.err()<<endl;

    //smeared
    //set time interval for eff_mass_fit SM
    /*
    if(data_2pts.Tag[iens].substr(1,1) =="A") {Corr.Tmin=24; Corr.Tmax=35;}
    else if(data_2pts.Tag[iens].substr(1,1) =="B") {Corr.Tmin=20; Corr.Tmax=36;}
    else if(data_2pts.Tag[iens].substr(1,1) == "C")  {Corr.Tmin=33;Corr.Tmax=51;}
    else if(data_2pts.Tag[iens].substr(1,1) == "D")  {Corr.Tmin=41;Corr.Tmax=54;}
    else crash("In fixing [Tmin, Tmax] for smeared MP, Ensemble: "+data_2pts.Tag[iens]+" not recognized");
    distr_t_list pt2_distr_SM= Corr.corr_t(data_2pts_SM.col(0)[iens], "../data/ph_emission_3d/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_SM.dat");
    distr_t_list eff_mass_SM = Corr.effective_mass_t(pt2_distr_SM, "../data/ph_emission_3d/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_SM.dat");
    distr_t M_P_SM = Corr.Fit_distr(eff_mass_SM);
    */

    //define renormalization factor for axial and vector currents in terms of axial 3pt at k=0
    //#################################################################################
    Corr.Perform_Nt_t_average=0;
    //sum 3pt axial at k=0 over ty
    distr_t_list ax_0 = qu*( Corr.corr_t(C_A_u_data[1][1][0].col(0)[iens],"") + Corr.corr_t(C_A_u_data[2][2][0].col(0)[iens],"") + Corr.corr_t(C_A_u_data[3][3][0].col(0)[iens],"") )  - qd*( Corr.corr_t(C_A_d_data[1][1][0].col(0)[iens],"") + Corr.corr_t(C_A_d_data[2][2][0].col(0)[iens],"") + Corr.corr_t(C_A_d_data[3][3][0].col(0)[iens],"")      );
    distr_t_list ax_0_u = qu*( Corr.corr_t(C_A_u_data[1][1][0].col(0)[iens],"") + Corr.corr_t(C_A_u_data[2][2][0].col(0)[iens],"") + Corr.corr_t(C_A_u_data[3][3][0].col(0)[iens],"") );
    distr_t_list ax_0_d = qd*( Corr.corr_t(C_A_d_data[1][1][0].col(0)[iens],"") + Corr.corr_t(C_A_d_data[2][2][0].col(0)[iens],"") + Corr.corr_t(C_A_d_data[3][3][0].col(0)[iens],"")      );
    distr_t ax_0_sum(UseJack, UseJack?Njacks:Nboots);
    distr_t ax_0_u_sum(UseJack, UseJack?Njacks:Nboots);
    distr_t ax_0_d_sum(UseJack, UseJack?Njacks:Nboots);
    for(int ty=0;ty<Corr.Nt;ty++) ax_0_sum = ax_0_sum + ax_0.distr_list[ty];
    for(int ty=0;ty<Corr.Nt;ty++) ax_0_u_sum = ax_0_u_sum + ax_0_u.distr_list[ty];
    for(int ty=0;ty<Corr.Nt;ty++) ax_0_d_sum = ax_0_d_sum + ax_0_d.distr_list[ty];
    distr_t renorm_A = FP/(ax_0_sum/3.0);
    distr_t FP_bare_3pt= (ax_0_sum/3.0);
    distr_t FP_bare_3pt_u= (ax_0_u_sum/3.0);
    distr_t FP_bare_3pt_d= (ax_0_d_sum/3.0);
    distr_t renorm_V = renorm_A*(Za/Zv);
    Corr.Perform_Nt_t_average=1;
    //#################################################################################


			  

    for(int ixg=1;ixg<n_xg;ixg++) {

      //get xg, Eg, kz from thetas
      double theta=thetas[ixg];
      pt3_momenta pt3_mom(0.0, 0.0, thetas[ixg]/2.0, masses_u[ixg], masses_d[ixg], 0.0, L_info.L, L_info.T);
      double Eg= pt3_mom.Egamma();
      double kz = pt3_mom.k()[2];


      //check if opposite theta is present
      bool theta_rev_present=false;
      int ixg_rev=-1;
      for(int loop_rev=0;loop_rev<n_xg_rev;loop_rev++)
	if( thetas_rev[loop_rev] == -1*theta ) {
	  theta_rev_present=true; ixg_rev=loop_rev; break;
	}
       
       

      cout<<"##### Considering kinematic with......"<<endl;
      cout<<"Eg: "<<Eg<<endl;
      cout<<"thz: "<<theta<<endl;
      cout<<"kz: "<<kz<<endl;
      cout<<"Opposite kz present: "<<(theta_rev_present?"true":"false")<<endl;

      distr_t renorm_V_w_kz= renorm_V/kz;

      PrecFloat precK=PrecFloat(kz);
      PrecFloat Ep_u= sqrt( pow(PrecFloat(1.249),2) + pow(precK,2));
      PrecFloat Lp_u= PrecFloat(3*0.421)/PrecFloat(10);

      PrecFloat Ep_d= sqrt( pow(PrecFloat(0.421),2) + pow(precK,2));
      PrecFloat Lp_d= PrecFloat(3*0.421)/PrecFloat(10);



      //########################################################################
      //########################################################################
      //########################################################################
      //@@@@@@@@ DEFINE CUSTOMIZED NORM TO BE USED IN HLT RECONSTRUCTION @@@@@@@
      
      auto F_NORM_u= [Ep=Ep_u,Lp=Lp_u ](PrecFloat E, PrecFloat m, PrecFloat s, PrecFloat E0, int jack_id) {

	PrecFloat res= pow(Lp*Lp/( pow(E-Ep,2)+ Lp*Lp),2);

	return res;

      };
      auto F_NORM_d= [Ep=Ep_d,Lp=Lp_d](PrecFloat E, PrecFloat m, PrecFloat s, PrecFloat E0, int jack_id) {

	
	return pow(Lp*Lp/( pow(E-Ep,2)+ Lp*Lp),2);

      };
      
      auto Atr_GEN_u= [Ep=Ep_u, Lp=Lp_u](PrecFloat TT) {
	
	
	PrecFloat MOD= sqrt( Lp*Lp + Ep*Ep)*TT;
	PrecFloat PH= Lp/Ep;
	PrecFloat s= Lp*TT;
	
	PrecFloat A1= 0.5*(Ep*Lp*Lp)/(Ep*Ep + Lp*Lp);
	PrecFloat A2= 0.5*Lp*exp(-Ep*TT)*precPi()*( cos(Lp*TT) + TT*Lp*sin(Lp*TT));
	PrecFloat A3= 0.5*Lp*exp(-Ep*TT)*( ExpEiComplexSum(MOD,PH,s,0) -Lp*TT*ExpEiComplexSum(MOD,PH,s,1));
	
	return A1+A2+A3;
	    

      };
      auto Atr_GEN_d= [Ep=Ep_d, Lp=Lp_d](PrecFloat TT) {
	    

	PrecFloat MOD= sqrt( Lp*Lp + Ep*Ep)*TT;
	PrecFloat PH= Lp/Ep;
	PrecFloat s= Lp*TT;

	PrecFloat A1= 0.5*(Ep*Lp*Lp)/(Ep*Ep + Lp*Lp);
	PrecFloat A2= 0.5*Lp*exp(-Ep*TT)*precPi()*( cos(Lp*TT) + TT*Lp*sin(Lp*TT));
	PrecFloat A3= 0.5*Lp*exp(-Ep*TT)*( ExpEiComplexSum(MOD,PH,s,0) -Lp*TT*ExpEiComplexSum(MOD,PH,s,1));

	return A1+A2+A3;

      };
      //########################################################################
      //########################################################################
      //########################################################################


      //vector
      for(auto &pair_V:red_mu_nu_pair_V) {

	int mu=pair_V.first;
	int nu=pair_V.second;

	double rev_theta_sign=-1;

	//vector
	int Im_Re;
	Corr.Reflection_sign = -1;
	Im_Re=1;
	Corr.Perform_Nt_t_average = 0;
	Corr_boot.Reflection_sign=-1;
	Corr_boot.Perform_Nt_t_average=0;
	 
	 
	distr_t_list vec_u = 0.5*qu*Corr.corr_t(summ_master(C_V_u_data[mu][nu][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_V_u_data[nu][mu][ixg].col(Im_Re)[iens], -1.0)),"../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"V_u_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg));
	distr_t_list vec_d = 0.5*qd*Corr.corr_t(summ_master(C_V_d_data[mu][nu][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_V_d_data[nu][mu][ixg].col(Im_Re)[iens], -1.0)) ,"../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"V_d_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg));

	

	distr_t_list vec_u_boot= 0.5*qu*Corr_boot.corr_t(summ_master(C_V_u_data[mu][nu][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_V_u_data[nu][mu][ixg].col(Im_Re)[iens], -1.0)),"");
	distr_t_list vec_d_boot= 0.5*qd*Corr_boot.corr_t(summ_master(C_V_d_data[mu][nu][ixg].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_V_d_data[nu][mu][ixg].col(Im_Re)[iens], -1.0)),"");


	if(theta_rev_present && Perform_theta_average) {

	  //jackknife
	  distr_t_list vec_u_rev=  rev_theta_sign*0.5*qu*Corr.corr_t(summ_master(C_V_u_data_rev[mu][nu][ixg_rev].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_V_u_data_rev[nu][mu][ixg_rev].col(Im_Re)[iens], -1.0)),"../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"V_u_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_rev");
	  distr_t_list vec_d_rev = rev_theta_sign*0.5*qd*Corr.corr_t(summ_master(C_V_d_data_rev[mu][nu][ixg_rev].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_V_d_data_rev[nu][mu][ixg_rev].col(Im_Re)[iens], -1.0)) ,"../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"V_d_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_rev");

	  //bootstrap
	  distr_t_list vec_u_boot_rev= rev_theta_sign*0.5*qu*Corr_boot.corr_t(summ_master(C_V_u_data_rev[mu][nu][ixg_rev].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_V_u_data[nu][mu][ixg_rev].col(Im_Re)[iens], -1.0)),"");
	  distr_t_list vec_d_boot_rev= rev_theta_sign*0.5*qd*Corr_boot.corr_t(summ_master(C_V_d_data_rev[mu][nu][ixg_rev].col(Im_Re)[iens], Multiply_Vvector_by_scalar(C_V_d_data[nu][mu][ixg_rev].col(Im_Re)[iens], -1.0)),"");

	  //compute difference between kz and -kz

	  distr_t_list vec_u_diff= vec_u - vec_u_rev;
	  distr_t_list vec_d_diff= vec_d - vec_d_rev;

	  
	  //average kz and -kz contributions

	  //jackknife
	  vec_u= 0.5*(vec_u + vec_u_rev);
	  vec_d= 0.5*(vec_d + vec_d_rev);

	  //bootstrap
	  vec_u_boot= 0.5*(vec_u_boot + vec_u_boot_rev);
	  vec_d_boot= 0.5*(vec_d_boot + vec_d_boot_rev);

	  //print averaged kz -kz
	  Print_To_File({}, { (vec_u/(0.5*qu)).ave(), (vec_u/(0.5*qu)).err(), (vec_u_diff/(0.5*qu)).ave(), (vec_u_diff/(0.5*qu)).err()}, "../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"V_u_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_ave.t", "", "");
	  Print_To_File({}, { (vec_d/(0.5*qd)).ave(), (vec_d/(0.5*qd)).err(), (vec_d_diff/(0.5*qd)).ave(), (vec_d_diff/(0.5*qd)).err()}, "../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"V_d_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_ave.t", "", "");
	    
	}
	

	vector<distr_t_list>  HV_u, HV_d, HV_tot;
	vector<distr_t_list>  HV_u_1_TO, HV_d_1_TO, HV_tot_1_TO;
	vector<distr_t_list>  HV_u_2_TO, HV_d_2_TO, HV_tot_2_TO;
	vector<distr_t_list>  HV_u_2_TO_w_sub, HV_d_2_TO_w_sub, HV_tot_2_TO_w_sub;


	for(int iv=0;iv<(signed)virt_list.size();iv++) {
	  HV_u.emplace_back( UseJack);
	  HV_d.emplace_back( UseJack);
	  HV_tot.emplace_back(UseJack);
	  HV_u_1_TO.emplace_back(UseJack);
	  HV_d_1_TO.emplace_back(UseJack);
	  HV_tot_1_TO.emplace_back(UseJack);
	  HV_u_2_TO.emplace_back(UseJack);
	  HV_d_2_TO.emplace_back(UseJack);
	  HV_tot_2_TO.emplace_back(UseJack);
	  HV_u_2_TO_w_sub.emplace_back(UseJack);
	  HV_d_2_TO_w_sub.emplace_back(UseJack);
	  HV_tot_2_TO_w_sub.emplace_back(UseJack);
	}
	 
	//standard integration
	Integrate_over_photon_insertion(vec_u, HV_u, Eg, t_weak,MP.ave(),0); 
	Integrate_over_photon_insertion(vec_d, HV_d, Eg, t_weak,MP.ave(),0);
	//first time ordering
	Integrate_over_photon_insertion(vec_u, HV_u_1_TO, Eg, t_weak,MP.ave(),1); 
	Integrate_over_photon_insertion(vec_d, HV_d_1_TO, Eg, t_weak,MP.ave(),1);
	//second time ordering
	Integrate_over_photon_insertion(vec_u, HV_u_2_TO, Eg, t_weak,MP.ave(),2); 
	Integrate_over_photon_insertion(vec_d, HV_d_2_TO, Eg, t_weak,MP.ave(),2);
	//second time ordering w substraction of exponential
	string obs_u="u_V_mu_"+to_string(mu)+"_nu_"+to_string(nu);
	string obs_d="d_V_mu_"+to_string(mu)+"_nu_"+to_string(nu);
	string out_path="../data/ph_emission_3d/"+ph_type_mes+"/mass/"+Ens_tags[iens];
	//set time intervals for masses and residue
	int Tmin_mass, Tmax_mass;
	if(ixg==1) { Tmin_mass=7; Tmax_mass=16; }
	else if(ixg==2) { Tmin_mass=7; Tmax_mass=16; }
	else if(ixg==3) { Tmin_mass=7; Tmax_mass=16; }
	else if(ixg==4) { Tmin_mass=7; Tmax_mass=16; }
	else crash("ixg: "+to_string(ixg)+" has no associate Tmin-Tmax interval for masses and residue in 2nd TO");
	Integrate_over_photon_insertion_w_subtraction(vec_u, HV_u_2_TO_w_sub, Eg, t_weak, MP.ave(), ixg, Tmin_mass, Tmax_mass,  out_path, obs_u);
	Integrate_over_photon_insertion_w_subtraction(vec_d, HV_d_2_TO_w_sub, Eg, t_weak, MP.ave(), ixg, Tmin_mass, Tmax_mass,  out_path, obs_d);
	
	
	
	 

	//loop over virtualities and renormalize contributions
	for(int iv=0;iv<(signed)virt_list.size();iv++) {
	  HV_u[iv] = renorm_V_w_kz*HV_u[iv];
	  HV_d[iv] = renorm_V_w_kz*HV_d[iv];
	  HV_u_1_TO[iv] = renorm_V_w_kz*HV_u_1_TO[iv];
	  HV_d_1_TO[iv] = renorm_V_w_kz*HV_d_1_TO[iv];
	  HV_u_2_TO[iv]= HV_u_2_TO[iv]*renorm_V_w_kz;
	  HV_d_2_TO[iv]= HV_d_2_TO[iv]*renorm_V_w_kz;
	  HV_u_2_TO_w_sub[iv] = HV_u_2_TO_w_sub[iv]*renorm_V_w_kz;
	  HV_d_2_TO_w_sub[iv] = HV_d_2_TO_w_sub[iv]*renorm_V_w_kz;
	  //sum ud contributions
	  HV_tot[iv]= HV_u[iv] +HV_d[iv];
	  HV_tot_1_TO[iv]= HV_u_1_TO[iv] +  HV_d_1_TO[iv];
	  HV_tot_2_TO[iv]= HV_u_2_TO[iv] +  HV_d_2_TO[iv];
	  HV_tot_2_TO_w_sub[iv] = HV_u_2_TO_w_sub[iv] + HV_d_2_TO_w_sub[iv];
	  
	}


	//print as a function of tcut for fixed virtuality
	for(int iv=0;iv<(signed)virt_list.size();iv++) {
	  //1+2 time orderings
	  Print_To_File({}, { HV_u[iv].ave(), HV_u[iv].err(), HV_d[iv].ave(), HV_d[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"V_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin   Vu  Vd");
	  Print_To_File({}, { HV_tot[iv].ave(), HV_tot[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"V_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin V");
	  //1 time ordering
	  Print_To_File({}, { HV_u_1_TO[iv].ave(), HV_u_1_TO[iv].err(), HV_d_1_TO[iv].ave(), HV_d_1_TO[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_1_V_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin   Vu  Vd");
	  Print_To_File({}, {   HV_tot_1_TO[iv].ave(), HV_tot_1_TO[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_1_V_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin V");
	  //2 time ordering
	  Print_To_File({}, {  HV_u_2_TO[iv].ave(), HV_u_2_TO[iv].err(), HV_d_2_TO[iv].ave(), HV_d_2_TO[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_2_V_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin   Vu  Vd");
	  Print_To_File({}, {   HV_tot_2_TO[iv].ave(), HV_tot_2_TO[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_2_V_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin V");
	  //2 time ordering w sub
	  Print_To_File({}, {  HV_u_2_TO_w_sub[iv].ave(), HV_u_2_TO_w_sub[iv].err(), HV_d_2_TO_w_sub[iv].ave(), HV_d_2_TO_w_sub[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"sub_V_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin   Vu  Vd");
	  Print_To_File({}, {   HV_tot_2_TO_w_sub[iv].ave(), HV_tot_2_TO_w_sub[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"sub_V_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin V");
	}

	//print as a function of virtuality for fixed tcut
	for(int tcut=0;tcut<Nts[iens];tcut++) {
	   
	  distr_t_list HV_u_tcut(UseJack), HV_d_tcut(UseJack), HV_tot_tcut(UseJack);
	  distr_t_list HV_u_1_TO_tcut(UseJack), HV_d_1_TO_tcut(UseJack), HV_tot_1_TO_tcut(UseJack);
	  distr_t_list HV_u_2_TO_tcut(UseJack), HV_d_2_TO_tcut(UseJack), HV_tot_2_TO_tcut(UseJack);
	  distr_t_list HV_u_2_TO_tcut_w_sub(UseJack), HV_d_2_TO_tcut_w_sub(UseJack), HV_tot_2_TO_tcut_w_sub(UseJack);
	   
	  for(int iv=0;iv<(signed)virt_list.size();iv++) {
	    HV_u_tcut.distr_list.push_back( HV_u[iv].distr_list[tcut]);
	    HV_d_tcut.distr_list.push_back( HV_d[iv].distr_list[tcut]);
	    HV_u_1_TO_tcut.distr_list.push_back( HV_u_1_TO[iv].distr_list[tcut]);
	    HV_d_1_TO_tcut.distr_list.push_back( HV_d_1_TO[iv].distr_list[tcut]);
	    HV_u_2_TO_tcut.distr_list.push_back( HV_u_2_TO[iv].distr_list[tcut]);
	    HV_d_2_TO_tcut.distr_list.push_back( HV_d_2_TO[iv].distr_list[tcut]);
	    if(tcut <= Nts[iens]/2) {
	    HV_u_2_TO_tcut_w_sub.distr_list.push_back( HV_u_2_TO_w_sub[iv].distr_list[tcut]);
	    HV_d_2_TO_tcut_w_sub.distr_list.push_back( HV_d_2_TO_w_sub[iv].distr_list[tcut]);
	    }
	  }
	  HV_tot_tcut= HV_u_tcut + HV_d_tcut;
	  HV_tot_1_TO_tcut= HV_u_1_TO_tcut + HV_d_1_TO_tcut;
	  HV_tot_2_TO_tcut= HV_u_2_TO_tcut + HV_d_2_TO_tcut;
	  if(tcut<=Nts[iens]/2) HV_tot_2_TO_tcut_w_sub = HV_u_2_TO_tcut_w_sub + HV_d_2_TO_tcut_w_sub;
	  //1+2 time orderings
	  Print_To_File({}, { virt_list, HV_u_tcut.ave(), HV_u_tcut.err(), HV_d_tcut.ave(), HV_d_tcut.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"V_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off   Vu  Vd");
	  Print_To_File({}, { virt_list, HV_tot_tcut.ave(), HV_tot_tcut.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"V_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off V");
	  //1 time ordering
	  Print_To_File({}, { virt_list, HV_u_1_TO_tcut.ave(), HV_u_1_TO_tcut.err(), HV_d_1_TO_tcut.ave(), HV_d_1_TO_tcut.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_1_V_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off   Vu  Vd");
	  Print_To_File({}, { virt_list,  HV_tot_1_TO_tcut.ave(), HV_tot_1_TO_tcut.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_1_V_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off V");
	  //2 time ordering
	  Print_To_File({}, { virt_list, HV_u_2_TO_tcut.ave(), HV_u_2_TO_tcut.err(), HV_d_2_TO_tcut.ave(), HV_d_2_TO_tcut.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_2_V_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off   Vu  Vd");
	  Print_To_File({}, { virt_list,  HV_tot_2_TO_tcut.ave(), HV_tot_2_TO_tcut.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_2_V_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off V");
	  //2 time ordering w sub
	  if(tcut<=Nts[iens]/2) {
	  Print_To_File({}, { virt_list, HV_u_2_TO_tcut_w_sub.ave(), HV_u_2_TO_tcut_w_sub.err(), HV_d_2_TO_tcut_w_sub.ave(), HV_d_2_TO_tcut_w_sub.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"sub_V_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off   Vu  Vd");
	  Print_To_File({}, { virt_list,  HV_tot_2_TO_tcut_w_sub.ave(), HV_tot_2_TO_tcut_w_sub.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"sub_V_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off V");
	  }
	   
	}
	 

	//define vec_u and vec_d in 2nd time ordering
	distr_t_list vec_u_TO_2(UseJack);
	distr_t_list vec_d_TO_2(UseJack);
	//boostrap
	distr_t_list vec_u_TO_2_boot(0);
	distr_t_list vec_d_TO_2_boot(0);

	for(int t=t_weak;t<=Nts[iens]/2;t++) {
	  vec_u_TO_2.distr_list.push_back( vec_u.distr_list[t]);
	  vec_d_TO_2.distr_list.push_back( vec_d.distr_list[t]);
	  vec_u_TO_2_boot.distr_list.push_back( vec_u_boot.distr_list[t]);
	  vec_d_TO_2_boot.distr_list.push_back( vec_d_boot.distr_list[t]);
	}

	int tmax_reco_u= vec_u_TO_2.size();
	int tmax_reco_d= vec_d_TO_2.size();
	 
	//generate covariance matrix
	Vfloat cov_vec_u, cov_vec_d, corr_vec_u, corr_vec_d;
	Vfloat TT, RR;
	for(int tt=0;tt< tmax_reco_u;tt++)
	  for(int rr=0;rr< tmax_reco_u;rr++) {
	    TT.push_back(tt);
	    RR.push_back(rr);
	    cov_vec_u.push_back( vec_u_TO_2_boot.distr_list[tt]%vec_u_TO_2_boot.distr_list[rr]);
	    cov_vec_d.push_back( vec_d_TO_2_boot.distr_list[tt]%vec_d_TO_2_boot.distr_list[rr]);
	    corr_vec_u.push_back( (vec_u_TO_2_boot.distr_list[tt]%vec_u_TO_2_boot.distr_list[rr])/(vec_u_TO_2_boot.err(tt)*vec_u_TO_2_boot.err(rr)));
	    corr_vec_d.push_back( (vec_d_TO_2_boot.distr_list[tt]%vec_d_TO_2_boot.distr_list[rr])/(vec_d_TO_2_boot.err(tt)*vec_d_TO_2_boot.err(rr)));
	  }


	//print covariance matrix
	Print_To_File({},{TT,RR, cov_vec_u, corr_vec_u}, "../data/ph_emission_3d/"+ph_type_mes+"/covariance/"+Ens_tags[iens]+"/"+TAG_CURR+"cov_Vu_mu_"+to_string(mu)+"_nu_"+to_string(nu)+".cov", "" , "");
	Print_To_File({},{TT,RR, cov_vec_d, corr_vec_d}, "../data/ph_emission_3d/"+ph_type_mes+"/covariance/"+Ens_tags[iens]+"/"+TAG_CURR+"cov_Vd_mu_"+to_string(mu)+"_nu_"+to_string(nu)+".cov", "" , "");

	if(!Skip_spectral_reconstruction && Reconstruct_vector_part) {

	
	 
	  //spectral reconstruction for second time ordering
	  for(int isg=0;isg<(signed)sigmas.size();isg++) {
	    cout<<"Calling spectral reconstruction 2nd-TO with sigma= "<<sigmas[isg]<<" GeV, vector channel, (mu,nu) : ("<<mu<<", "<<nu<<")"<<endl<<flush;
	    Vfloat syst_re_u(virt_list.size()), syst_re_d(virt_list.size()), syst_im_u(virt_list.size()), syst_im_d(virt_list.size());
	   
#pragma omp parallel for schedule(dynamic)
	    for(int ie=0;ie<(signed)virt_list.size();ie++) {
	   
	      double mult_re_u=1e-2;
	      double mult_re_d=1e-2;
	      double mult_im_u=1e-2;
	      double mult_im_d=1e-2;
	      double s= sigmas[isg]*a_distr.ave();
	      double Eg_virt= sqrt( Eg*Eg + pow(MP.ave()*virt_list[ie],2));
	      double E0_u_RE= E0_fact_u*sqrt( pow(Mjpsi*a_distr.ave(),2) + pow(kz,2));
	      double E0_d_RE= E0_fact_d*sqrt( pow(Mphi*a_distr.ave(),2) + pow(kz,2));
	      double E0_u_IM= min( Eg+3*s, E0_u_RE);
	      double E0_d_IM= min( Eg+3*s, E0_d_RE);
	      double l_re_u, l_re_d;
	      double l_im_u, l_im_d;

	      cout<<"Computing ixg: "<<ixg<<" xk: "<<virt_list[ie]<<" Eg: "<<Eg_virt<<" E0_fact_u: "<<E0_fact_u<<", E0_fact_d: "<<E0_fact_d<<" sigma: "<<sigmas[isg]<<" SM_TYPE: "<<SM_TYPE<<" CONS CURRENT: "<<CONS_EM_CURR<<endl<<flush;
	      cout<<"MV_u*a "<<Mjpsi*a_distr.ave()<<" MV_d*a "<<Mphi*a_distr.ave()<<endl<<flush;
	     
	      //Real part
	      syst_re_u[ie]=0.0;
	      RE_HV_sm_u[iens][mu][nu][ixg][isg].distr_list[ie] = 0.0*Get_id_jack_distr(Njacks);
	      //RE_HV_sm_u[iens][mu][nu][ixg][isg].distr_list[ie] = Get_Laplace_transfo(  Eg_virt,  s, E0_u_RE,  Nts[iens], tmax_reco_u-1, prec, SM_TYPE+"_RE",K_RE, vec_u_TO_2, syst_re_u[ie], mult_re_u, l_re_u, MODE_FF, "Ef_"+to_string_with_precision(E0_fact_u,1), TAG_CURR+"Vu_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_"+Ens_tags[iens], -1,0, renorm_V_w_kz, "virtual_FF", cov_vec_u, fake_func,0, fake_func_d ,  Integrate_Up_To_Emax, Emax, beta, 1,0, F_NORM_u, Atr_GEN_u );
	      RE_HV_sm_d[iens][mu][nu][ixg][isg].distr_list[ie] = Get_Laplace_transfo(  Eg_virt,  s, E0_d_RE,  Nts[iens], tmax_reco_d-1, prec, SM_TYPE+"_RE",K_RE, vec_d_TO_2, syst_re_d[ie], mult_re_d, l_re_d, MODE_FF, "Ef_"+to_string_with_precision(E0_fact_d,1), TAG_CURR+"Vd_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_"+Ens_tags[iens], -1,0, renorm_V_w_kz, "virtual_FF", cov_vec_d, fake_func,0, fake_func_d ,  Integrate_Up_To_Emax, Emax, beta, 1,0, F_NORM_d, Atr_GEN_d);

	      //Imag part
	      syst_im_u[ie]= 0.0;
	      IM_HV_sm_u[iens][mu][nu][ixg][isg].distr_list[ie] = 0.0*Get_id_jack_distr(Njacks) ;
	      //IM_HV_sm_u[iens][mu][nu][ixg][isg].distr_list[ie] = Get_Laplace_transfo(  Eg_virt,  s, E0_u,  Nts[iens], tmax_reco_u-1, prec, SM_TYPE+"_IM",K_IM, vec_u_TO_2, syst_im_u[ie], mult_im_u, l_im_u, MODE_FF, "Ef_"+to_string_with_precision(E0_fact_u,1), TAG_CURR+"Vu_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_"+Ens_tags[iens], -1,0, renorm_V_w_kz, "virtual_FF", cov_vec_u, fake_func,0, fake_func_d ,  Integrate_Up_To_Emax, Emax, beta, 1,0, F_NORM_u, Atr_GEN_u);
	      IM_HV_sm_d[iens][mu][nu][ixg][isg].distr_list[ie] = Get_Laplace_transfo(  Eg_virt,  s, E0_d_IM,  Nts[iens], tmax_reco_d-1, prec, SM_TYPE+"_IM",K_IM, vec_d_TO_2, syst_im_d[ie], mult_im_d, l_im_d, MODE_FF, "Ef_"+to_string_with_precision(E0_fact_d,1), TAG_CURR+"Vd_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_"+Ens_tags[iens], -1,0, renorm_V_w_kz, "virtual_FF", cov_vec_d, fake_func,0, fake_func_d ,  Integrate_Up_To_Emax, Emax, beta, 1,0, F_NORM_d, Atr_GEN_d);


	      cout<<"Computed ixg: "<<ixg<<" xk: "<<virt_list[ie]<<" Eg: "<<Eg_virt<<" E0_fact_u: "<<E0_fact_u<<", E0_fact_d: "<<E0_fact_d<<" sigma: "<<sigmas[isg]<<" SM_TYPE: "<<SM_TYPE<<" CONS CURRENT: "<<CONS_EM_CURR<<endl<<flush;


	    }

	    //print to file
	    //Real part
	    Print_To_File({}, {virt_list, RE_HV_sm_u[iens][mu][nu][ixg][isg].ave(), RE_HV_sm_u[iens][mu][nu][ixg][isg].err(), syst_re_u, RE_HV_sm_d[iens][mu][nu][ixg][isg].ave(), RE_HV_sm_d[iens][mu][nu][ixg][isg].err(), syst_re_d, (RE_HV_sm_u[iens][mu][nu][ixg][isg] + RE_HV_sm_d[iens][mu][nu][ixg][isg]).ave(), (RE_HV_sm_u[iens][mu][nu][ixg][isg] + RE_HV_sm_d[iens][mu][nu][ixg][isg]).err()}, "../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+"RE_V_alpha_"+to_string_with_precision(beta,2)+"_E0_"+to_string_with_precision(E0_fact,2)+"_SM_"+SM_TYPE+"_ixg_"+to_string(ixg)+"_sigma_"+to_string_with_precision(sigmas[isg],3)+"_mu_"+to_string(mu)+"_nu_"+to_string(nu)+".dat", "", "xk u d  u+d");
	    //Imag part
	    Print_To_File({}, {virt_list, IM_HV_sm_u[iens][mu][nu][ixg][isg].ave(), IM_HV_sm_u[iens][mu][nu][ixg][isg].err(), syst_im_u, IM_HV_sm_d[iens][mu][nu][ixg][isg].ave(), IM_HV_sm_d[iens][mu][nu][ixg][isg].err(), syst_im_d, (IM_HV_sm_u[iens][mu][nu][ixg][isg]+ IM_HV_sm_d[iens][mu][nu][ixg][isg]).ave(), (IM_HV_sm_u[iens][mu][nu][ixg][isg]+ IM_HV_sm_d[iens][mu][nu][ixg][isg]).err() }, "../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+"IM_V_alpha_"+to_string_with_precision(beta,2)+"_E0_"+to_string_with_precision(E0_fact,2)+"_SM_"+SM_TYPE+"_ixg_"+to_string(ixg)+"_sigma_"+to_string_with_precision(sigmas[isg],3)+"_mu_"+to_string(mu)+"_nu_"+to_string(nu)+".dat", "", "xk u d u+d");

	    cout<<"sigma: "<<sigmas[isg]<<" computed!"<<endl;
	   
	  }
	
	cout<<"done!"<<endl<<flush;
	}
      }

      
      
      //axial
      for(auto &pair_A:red_mu_nu_pair_A) {

	int mu=pair_A.first;
	int nu=pair_A.second;
	double rev_theta_sign= ( (mu==0 || nu==0) && (mu != 0 || nu != 0))?-1.0:1.0;

	//axial
	int Im_Re;
	double parity;
	Corr.Perform_Nt_t_average=0;
	Corr_boot.Perform_Nt_t_average=0;
	if( (mu==0 || nu==0) && (mu != 0 || nu != 0)) {Im_Re=1; Corr.Reflection_sign=-1; Corr_boot.Reflection_sign=-1;  parity=1.0;}
	else { Im_Re=0; Corr.Reflection_sign=1; Corr_boot.Reflection_sign=1; parity=1.0;}

	distr_t_list ax_u(UseJack), ax_d(UseJack), ax_u_boot(0), ax_d_boot(0);

	if((mu != 1) || (nu != 1)) {

	  ax_u = parity*qu*Corr.corr_t(C_A_u_data[mu][nu][ixg].col(Im_Re)[iens],"../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"A_u_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg));
	  ax_d = parity*qd*Corr.corr_t(C_A_d_data[mu][nu][ixg].col(Im_Re)[iens],"../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"A_d_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg));

	  ax_u_boot= parity*qu*Corr_boot.corr_t(C_A_u_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  ax_d_boot= parity*qd*Corr_boot.corr_t(C_A_d_data[mu][nu][ixg].col(Im_Re)[iens],"");

	}
	else {
	  ax_u = 0.5*parity*qu*Corr.corr_t(summ_master(C_A_u_data[mu][nu][ixg].col(Im_Re)[iens], C_A_u_data[mu+1][nu+1][ixg].col(Im_Re)[iens]),"../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"A_u_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg));
	  ax_d = 0.5*parity*qd*Corr.corr_t(summ_master(C_A_d_data[mu][nu][ixg].col(Im_Re)[iens], C_A_d_data[mu+1][nu+1][ixg].col(Im_Re)[iens]) ,"../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"A_d_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg));

	  ax_u_boot= 0.5*parity*qu*Corr_boot.corr_t(summ_master(C_A_u_data[mu][nu][ixg].col(Im_Re)[iens], C_A_u_data[mu+1][nu+1][ixg].col(Im_Re)[iens]),"");
	  ax_d_boot= 0.5*parity*qd*Corr_boot.corr_t(summ_master(C_A_d_data[mu][nu][ixg].col(Im_Re)[iens], C_A_d_data[mu+1][nu+1][ixg].col(Im_Re)[iens]),"");

	}


	if(theta_rev_present && Perform_theta_average) {

	  distr_t_list ax_u_rev(UseJack), ax_d_rev(UseJack);
	  distr_t_list ax_u_boot_rev(0), ax_d_boot_rev(0);

	  if((mu != 1) || (nu != 1)) {

	    ax_u_rev = rev_theta_sign*parity*qu*Corr.corr_t(C_A_u_data_rev[mu][nu][ixg_rev].col(Im_Re)[iens],"../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"A_u_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_rev");
	    ax_d_rev = rev_theta_sign*parity*qd*Corr.corr_t(C_A_d_data_rev[mu][nu][ixg_rev].col(Im_Re)[iens],"../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"A_d_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_rev");

	    ax_u_boot_rev= rev_theta_sign*parity*qu*Corr_boot.corr_t(C_A_u_data_rev[mu][nu][ixg_rev].col(Im_Re)[iens],"");
	    ax_d_boot_rev= rev_theta_sign*parity*qd*Corr_boot.corr_t(C_A_d_data_rev[mu][nu][ixg_rev].col(Im_Re)[iens],"");

	  }
	  else {
	    ax_u_rev = rev_theta_sign*0.5*parity*qu*Corr.corr_t(summ_master(C_A_u_data_rev[mu][nu][ixg_rev].col(Im_Re)[iens], C_A_u_data_rev[mu+1][nu+1][ixg_rev].col(Im_Re)[iens]),"../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"A_u_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_rev");
	    ax_d_rev = rev_theta_sign*0.5*parity*qd*Corr.corr_t(summ_master(C_A_d_data_rev[mu][nu][ixg_rev].col(Im_Re)[iens], C_A_d_data_rev[mu+1][nu+1][ixg_rev].col(Im_Re)[iens]) ,"../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"A_d_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_rev");
	     
	    ax_u_boot_rev= rev_theta_sign*0.5*parity*qu*Corr_boot.corr_t(summ_master(C_A_u_data_rev[mu][nu][ixg_rev].col(Im_Re)[iens], C_A_u_data_rev[mu+1][nu+1][ixg_rev].col(Im_Re)[iens]),"");
	    ax_d_boot_rev= rev_theta_sign*0.5*parity*qd*Corr_boot.corr_t(summ_master(C_A_d_data_rev[mu][nu][ixg_rev].col(Im_Re)[iens], C_A_d_data_rev[mu+1][nu+1][ixg_rev].col(Im_Re)[iens]),"");

	  }
	    
	  //compute difference between kz and -kz

	  distr_t_list ax_u_diff= ax_u - ax_u_rev;
	  distr_t_list ax_d_diff= ax_d - ax_d_rev;
	  
	  //average kz and -kz contributions

	  //jackknife
	  ax_u= 0.5*(ax_u + ax_u_rev);
	  ax_d= 0.5*(ax_d + ax_d_rev);

	  

	  //bootstrap
	  ax_u_boot= 0.5*(ax_u_boot + ax_u_boot_rev);
	  ax_d_boot= 0.5*(ax_d_boot + ax_d_boot_rev);

	  double resc_fact= (mu==1 && nu==1)?0.5:1.0;

	  //print averaged kz -kz
	  Print_To_File({}, { (ax_u/(parity*resc_fact*qu)).ave(), (ax_u/(parity*resc_fact*qu)).err(), (ax_u_diff/(parity*resc_fact*qu)).ave(), (ax_u_diff/(parity*resc_fact*qu)).err()}, "../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"A_u_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_ave.t", "", "");
	  Print_To_File({}, { (ax_d/(parity*resc_fact*qd)).ave(), (ax_d/(parity*resc_fact*qd)).err(), (ax_d_diff/(parity*resc_fact*qd)).ave(), (ax_d_diff/(parity*resc_fact*qd)).err()}, "../data/ph_emission_3d/"+ph_type_mes+"/C/"+Ens_tags[iens]+"/"+TAG_CURR+"A_d_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_ave.t", "", "");
	    
	}

	
	vector<distr_t_list>  HA_u, HA_d, HA_tot;
	vector<distr_t_list>  HA_u_1_TO, HA_d_1_TO, HA_tot_1_TO;
	vector<distr_t_list>  HA_u_2_TO, HA_d_2_TO, HA_tot_2_TO;
	vector<distr_t_list>  HA_u_2_TO_w_sub, HA_d_2_TO_w_sub, HA_tot_2_TO_w_sub;

	for(int iv=0;iv<(signed)virt_list.size();iv++) {
	  HA_u.emplace_back( UseJack);
	  HA_d.emplace_back( UseJack);
	  HA_tot.emplace_back(UseJack);
	  HA_u_1_TO.emplace_back(UseJack);
	  HA_d_1_TO.emplace_back(UseJack);
	  HA_tot_1_TO.emplace_back(UseJack);
	  HA_u_2_TO.emplace_back(UseJack);
	  HA_d_2_TO.emplace_back(UseJack);
	  HA_tot_2_TO.emplace_back(UseJack);
	  HA_u_2_TO_w_sub.emplace_back(UseJack);
	  HA_d_2_TO_w_sub.emplace_back(UseJack);
	  HA_tot_2_TO_w_sub.emplace_back(UseJack);
	}

	distr_t_list kin_fact_point_sub(UseJack);

	if( (pair_A == make_pair(1,1)) || (pair_A == make_pair(2,2))) {
	  for(auto &virt: virt_list)  kin_fact_point_sub.distr_list.push_back(FP/FP);}
	else if( pair_A == make_pair(3,0)) {
	  for(auto &virt: virt_list) {
	    double off2= pow(MP.ave()*virt,2);
	    double Eg_virt= sqrt( Eg*Eg+ off2);	     
	    kin_fact_point_sub.distr_list.push_back( 0.0*kz*(MP-Eg_virt)/(2*MP*Eg_virt -off2) );
	  }
	}
	else if( pair_A == make_pair(0,3)) {
	  for(auto &virt: virt_list) {
	    double off2= pow(MP.ave()*virt,2);
	    double Eg_virt= sqrt( Eg*Eg+ off2);
	    
	    kin_fact_point_sub.distr_list.push_back( 0.0*kz*(2*MP - Eg_virt)/(2*MP*Eg_virt - off2) );
	  }
	}
	else if( pair_A == make_pair(3,3)) {
	  for(auto &virt: virt_list) {
	    double off2= pow(MP.ave()*virt,2);
	    double Eg_virt= sqrt( Eg*Eg+ off2);
	    kin_fact_point_sub.distr_list.push_back( Eg_virt*( 2*MP-Eg_virt)/(2*MP*Eg_virt - off2));
	  }
	}
	else crash("(mu,nu) = ("+to_string(mu)+", "+to_string(nu)+") is not a valid pair in axial channel");


		 
	//standard integration
	Integrate_over_photon_insertion(ax_u, HA_u, Eg, t_weak, MP.ave(), 0); 
	Integrate_over_photon_insertion(ax_d, HA_d, Eg, t_weak,MP.ave(),0);
	//first time ordering
	Integrate_over_photon_insertion(ax_u, HA_u_1_TO, Eg, t_weak, MP.ave(), 1); 
	Integrate_over_photon_insertion(ax_d, HA_d_1_TO, Eg, t_weak,MP.ave(),1);
	//second time ordering
	Integrate_over_photon_insertion(ax_u, HA_u_2_TO, Eg, t_weak,MP.ave(),2); 
	Integrate_over_photon_insertion(ax_d, HA_d_2_TO, Eg, t_weak,MP.ave(),2);
	//second time ordering w substraction of exponential
	string obs_u="u_A_mu_"+to_string(mu)+"_nu_"+to_string(nu);
	string obs_d="d_A_mu_"+to_string(mu)+"_nu_"+to_string(nu);
	string out_path="../data/ph_emission_3d/"+ph_type_mes+"/mass/"+Ens_tags[iens];
	//set time intervals for masses and residue
	int Tmin_mass, Tmax_mass;
	if(ixg==1) { Tmin_mass=7; Tmax_mass=16; }
	else if(ixg==2) { Tmin_mass=7; Tmax_mass=16;}
	else if(ixg==3) { Tmin_mass=7; Tmax_mass=16;}
	else if(ixg==4) { Tmin_mass=7; Tmax_mass=16;}
	else crash("ixg: "+to_string(ixg)+" has no associate Tmin-Tmax interval for masses and residue in 2nd TO");
	Integrate_over_photon_insertion_w_subtraction(ax_u, HA_u_2_TO_w_sub, Eg, t_weak, MP.ave(), ixg, Tmin_mass, Tmax_mass, out_path, obs_u);
	Integrate_over_photon_insertion_w_subtraction(ax_d, HA_d_2_TO_w_sub, Eg, t_weak, MP.ave(), ixg, Tmin_mass, Tmax_mass, out_path, obs_d);

       	 
	//loop over virtualities and renormalize contributions
	for(int iv=0;iv<(signed)virt_list.size();iv++) {
	  HA_u[iv] = renorm_A*(HA_u[iv] - FP_bare_3pt_u*kin_fact_point_sub.distr_list[iv]);
	  HA_d[iv] = renorm_A*(HA_d[iv] - FP_bare_3pt_d*kin_fact_point_sub.distr_list[iv]);
	  HA_u_1_TO[iv] = renorm_A*( HA_u_1_TO[iv] - FP_bare_3pt_u*kin_fact_point_sub.distr_list[iv]);
	  HA_d_1_TO[iv] = renorm_A*( HA_d_1_TO[iv] - FP_bare_3pt_d*kin_fact_point_sub.distr_list[iv]);
	  HA_u_2_TO[iv]= HA_u_2_TO[iv]*renorm_A;
	  HA_d_2_TO[iv]= HA_d_2_TO[iv]*renorm_A;
	  HA_u_2_TO_w_sub[iv] = HA_u_2_TO_w_sub[iv]*renorm_A;
	  HA_d_2_TO_w_sub[iv] = HA_d_2_TO_w_sub[iv]*renorm_A;
	  //sum ud contributions
	  HA_tot[iv]= HA_u[iv] -HA_d[iv];
	  HA_tot_1_TO[iv]= HA_u_1_TO[iv] - HA_d_1_TO[iv];
	  HA_tot_2_TO[iv]= HA_u_2_TO[iv] - HA_d_2_TO[iv];
	  HA_tot_2_TO_w_sub[iv] = HA_u_2_TO_w_sub[iv] - HA_d_2_TO_w_sub[iv];

	}



	//print as a function of tcut for fixed virtuality
	for(int iv=0;iv<(signed)virt_list.size();iv++) {
	  //1+2 time orderings
	  Print_To_File({}, { HA_u[iv].ave(), HA_u[iv].err(), HA_d[iv].ave(), HA_d[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"A_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin   Au  Ad");
	  Print_To_File({}, { HA_tot[iv].ave(), HA_tot[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"A_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin A");
	  //1 time ordering
	  Print_To_File({}, { HA_u_1_TO[iv].ave(), HA_u_1_TO[iv].err(), HA_d_1_TO[iv].ave(), HA_d_1_TO[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_1_A_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin   Au  Ad");
	  Print_To_File({}, {   HA_tot_1_TO[iv].ave(), HA_tot_1_TO[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_1_A_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin A");
	  //2 time ordering
	  Print_To_File({}, {  HA_u_2_TO[iv].ave(), HA_u_2_TO[iv].err(), HA_d_2_TO[iv].ave(), HA_d_2_TO[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_2_A_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin   Au  Ad");
	  Print_To_File({}, {   HA_tot_2_TO[iv].ave(), HA_tot_2_TO[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_2_A_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin A");
	  //2 time ordering w sub
	  Print_To_File({}, {  HA_u_2_TO_w_sub[iv].ave(), HA_u_2_TO_w_sub[iv].err(), HA_d_2_TO_w_sub[iv].ave(), HA_d_2_TO_w_sub[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"sub_A_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin   Au  Ad");
	  Print_To_File({}, {   HA_tot_2_TO_w_sub[iv].ave(), HA_tot_2_TO_w_sub[iv].err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"sub_A_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_ixk_"+to_string(iv), "", "#tmin A");
	  
	}

	//print as a function of virtuality for fixed tcut
	for(int tcut=0;tcut<Nts[iens];tcut++) {
	   
	  distr_t_list HA_u_tcut(UseJack), HA_d_tcut(UseJack), HA_tot_tcut(UseJack);
	  distr_t_list HA_u_1_TO_tcut(UseJack), HA_d_1_TO_tcut(UseJack), HA_tot_1_TO_tcut(UseJack);
	  distr_t_list HA_u_2_TO_tcut(UseJack), HA_d_2_TO_tcut(UseJack), HA_tot_2_TO_tcut(UseJack);
	  distr_t_list HA_u_2_TO_tcut_w_sub(UseJack), HA_d_2_TO_tcut_w_sub(UseJack), HA_tot_2_TO_tcut_w_sub(UseJack);

	  for(int iv=0;iv<(signed)virt_list.size();iv++) {
	    HA_u_tcut.distr_list.push_back( HA_u[iv].distr_list[tcut]);
	    HA_d_tcut.distr_list.push_back( HA_d[iv].distr_list[tcut]);
	    HA_u_1_TO_tcut.distr_list.push_back( HA_u_1_TO[iv].distr_list[tcut]);
	    HA_d_1_TO_tcut.distr_list.push_back( HA_d_1_TO[iv].distr_list[tcut]);
	    HA_u_2_TO_tcut.distr_list.push_back( HA_u_2_TO[iv].distr_list[tcut]);
	    HA_d_2_TO_tcut.distr_list.push_back( HA_d_2_TO[iv].distr_list[tcut]);
	    if(tcut<= Nts[iens]/2) {
	    HA_u_2_TO_tcut_w_sub.distr_list.push_back( HA_u_2_TO_w_sub[iv].distr_list[tcut]);
	    HA_d_2_TO_tcut_w_sub.distr_list.push_back( HA_d_2_TO_w_sub[iv].distr_list[tcut]);
	    }
	  }
	  HA_tot_tcut= HA_u_tcut - HA_d_tcut;
	  HA_tot_1_TO_tcut= HA_u_1_TO_tcut- HA_d_1_TO_tcut;
	  HA_tot_2_TO_tcut= HA_u_2_TO_tcut - HA_d_2_TO_tcut;
	  if(tcut<= Nts[iens]/2) HA_tot_2_TO_tcut_w_sub= HA_u_2_TO_tcut_w_sub - HA_d_2_TO_tcut_w_sub;
	  //1+2 time orderings
	  Print_To_File({}, { virt_list, HA_u_tcut.ave(), HA_u_tcut.err(), HA_d_tcut.ave(), HA_d_tcut.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"A_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off   Au  Ad");
	  Print_To_File({}, { virt_list, HA_tot_tcut.ave(), HA_tot_tcut.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"A_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off A");
	  //1 time ordering
	  Print_To_File({}, { virt_list, HA_u_1_TO_tcut.ave(), HA_u_1_TO_tcut.err(), HA_d_1_TO_tcut.ave(), HA_d_1_TO_tcut.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_1_A_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off   Au  Ad");
	  Print_To_File({}, { virt_list,  HA_tot_1_TO_tcut.ave(), HA_tot_1_TO_tcut.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_1_A_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off A");
	  //2 time ordering
	  Print_To_File({}, { virt_list, HA_u_2_TO_tcut.ave(), HA_u_2_TO_tcut.err(), HA_d_2_TO_tcut.ave(), HA_d_2_TO_tcut.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_2_A_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off   Au  Ad");
	  Print_To_File({}, { virt_list,  HA_tot_2_TO_tcut.ave(), HA_tot_2_TO_tcut.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"TO_2_A_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off A");
	  //2 time ordering w sub
	  if(tcut<=Nts[iens]/2) {
	   Print_To_File({}, { virt_list, HA_u_2_TO_tcut_w_sub.ave(), HA_u_2_TO_tcut_w_sub.err(), HA_d_2_TO_tcut_w_sub.ave(), HA_d_2_TO_tcut_w_sub.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"sub_A_quark_contr_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_ixg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off   Au  Ad");
	   Print_To_File({}, { virt_list,  HA_tot_2_TO_tcut_w_sub.ave(), HA_tot_2_TO_tcut_w_sub.err()}, "../data/ph_emission_3d/"+ph_type_mes+"/H/"+Ens_tags[iens]+"/"+TAG_CURR+"sub_A_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_xg_"+to_string(ixg)+"_tcut_"+to_string(tcut), "", "#off A");
	  }

	}

		 
	
	//define ax_u and ax_d in 2nd time ordering
	distr_t_list ax_u_TO_2(UseJack);
	distr_t_list ax_d_TO_2(UseJack);
	//bootstrap
	distr_t_list ax_u_TO_2_boot(0);
	distr_t_list ax_d_TO_2_boot(0);

	for(int t=t_weak;t<=Nts[iens]/2;t++) {
	  ax_u_TO_2.distr_list.push_back( ax_u.distr_list[t]);
	  ax_d_TO_2.distr_list.push_back( ax_d.distr_list[t]);
	  ax_u_TO_2_boot.distr_list.push_back( ax_u_boot.distr_list[t]);
	  ax_d_TO_2_boot.distr_list.push_back( ax_d_boot.distr_list[t]);
	}

	int tmax_reco_u= ax_u_TO_2.size();
	int tmax_reco_d= ax_d_TO_2.size();
	 
	//generate covariance matrix
	Vfloat cov_ax_u, cov_ax_d, corr_ax_u, corr_ax_d;
	Vfloat TT, RR;
	for(int tt=0;tt<tmax_reco_u;tt++)
	  for(int rr=0;rr<tmax_reco_u;rr++) {
	    TT.push_back(tt);
	    RR.push_back(rr);
	    cov_ax_u.push_back( ax_u_TO_2_boot.distr_list[tt]%ax_u_TO_2_boot.distr_list[rr]);
	    cov_ax_d.push_back( ax_d_TO_2_boot.distr_list[tt]%ax_d_TO_2_boot.distr_list[rr]);
	    corr_ax_u.push_back( (ax_u_TO_2_boot.distr_list[tt]%ax_u_TO_2_boot.distr_list[rr])/(ax_u_TO_2_boot.err(tt)*ax_u_TO_2_boot.err(rr)));
	    corr_ax_d.push_back( (ax_d_TO_2_boot.distr_list[tt]%ax_d_TO_2_boot.distr_list[rr])/(ax_d_TO_2_boot.err(tt)*ax_d_TO_2_boot.err(rr)));
	  }

	//print covariance matrix
	Print_To_File({},{TT,RR, cov_ax_u, corr_ax_u}, "../data/ph_emission_3d/"+ph_type_mes+"/covariance/"+Ens_tags[iens]+"/"+TAG_CURR+"cov_Au_mu_"+to_string(mu)+"_nu_"+to_string(nu)+".cov", "" , "");
	Print_To_File({},{TT,RR, cov_ax_d, corr_ax_d}, "../data/ph_emission_3d/"+ph_type_mes+"/covariance/"+Ens_tags[iens]+"/"+TAG_CURR+"cov_Ad_mu_"+to_string(mu)+"_nu_"+to_string(nu)+".cov", "" , "");

	if(!Skip_spectral_reconstruction && Reconstruct_axial_part) {
		 
	  //spectral reconstruction for second time ordering
	  for(int isg=0;isg<(signed)sigmas.size();isg++) {
	    cout<<"Calling spectral reconstruction 2nd-TO with sigma= "<<sigmas[isg]<<" GeV, axial channel, (mu,nu) : ("<<mu<<", "<<nu<<")"<<endl<<flush;
	    Vfloat syst_re_u(virt_list.size()), syst_re_d(virt_list.size()), syst_im_u(virt_list.size()), syst_im_d(virt_list.size());

#pragma omp parallel for schedule(dynamic)
	    for(int ie=0;ie<(signed)virt_list.size();ie++) {


	      
	      double mult_re_u=1e-2;
	      double mult_re_d=1e-2;
	      double mult_im_u=1e-2;
	      double mult_im_d=1e-2;
	      double s= sigmas[isg]*a_distr.ave();
	      double Eg_virt= sqrt( Eg*Eg + pow(MP.ave()*virt_list[ie],2));
	      double E0_u_RE= E0_fact_u*sqrt( pow(Mjpsi*a_distr.ave(),2) + pow(kz,2));
	      double E0_d_RE= E0_fact_d*sqrt( pow(Mphi*a_distr.ave(),2) + pow(kz,2));
	      double E0_u_IM= min( Eg + 3*s, E0_u_RE);
	      double E0_d_IM= min( Eg + 3*s, E0_d_RE);
	      cout<<"Eg: "<<Eg_virt<<" sigma: "<<s<<endl;
	      double l_re_u, l_re_d;
	      double l_im_u, l_im_d;

	      cout<<"Computing ixg: "<<ixg<<" xk: "<<virt_list[ie]<<" Eg: "<<Eg_virt<<" E0_fact_u: "<<E0_fact_u<<", E0_fact_d: "<<E0_fact_d<<" sigma: "<<sigmas[isg]<<" SM_TYPE: "<<SM_TYPE<<" CONS CURRENT: "<<CONS_EM_CURR<<endl<<flush;
              cout<<"MV_u*a "<<Mjpsi*a_distr.ave()<<" MV_d*a "<<Mphi*a_distr.ave()<<endl<<flush;
	     

	      //Real part
	      syst_re_u[ie] = 0;
	      RE_HA_sm_u[iens][mu][nu][ixg][isg].distr_list[ie] = 0.0*Get_id_jack_distr(Njacks);
	      //RE_HA_sm_u[iens][mu][nu][ixg][isg].distr_list[ie] = Get_Laplace_transfo(  Eg_virt,  s, E0_u_RE,  Nts[iens], tmax_reco_u-1, prec, SM_TYPE+"_RE",K_RE, ax_u_TO_2, syst_re_u[ie], mult_re_u, l_re_u, MODE_FF, "Ef_"+to_string_with_precision(E0_fact_u,1), TAG_CURR+"Au_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_"+Ens_tags[iens], -1,0, renorm_A, "virtual_FF", cov_ax_u, fake_func,0, fake_func_d ,  Integrate_Up_To_Emax, Emax, beta, 1,0, F_NORM_u, Atr_GEN_u);
	      cout<<"Re HA u, computed"<<endl;
	      RE_HA_sm_d[iens][mu][nu][ixg][isg].distr_list[ie] = Get_Laplace_transfo(  Eg_virt,  s, E0_d_RE,  Nts[iens], tmax_reco_d-1, prec, SM_TYPE+"_RE",K_RE, ax_d_TO_2, syst_re_d[ie], mult_re_d, l_re_d, MODE_FF, "Ef_"+to_string_with_precision(E0_fact_d,1), TAG_CURR+"Ad_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_"+Ens_tags[iens], -1,0, renorm_A, "virtual_FF", cov_ax_d, fake_func,0, fake_func_d ,  Integrate_Up_To_Emax, Emax, beta, 1,0, F_NORM_d, Atr_GEN_d);
	      cout<<"Re HA d, computed"<<endl;
	      
	      //Imag part
	      syst_im_u[ie] = 0;
	      IM_HA_sm_u[iens][mu][nu][ixg][isg].distr_list[ie]= 0.0*Get_id_jack_distr(Njacks);
	      //IM_HA_sm_u[iens][mu][nu][ixg][isg].distr_list[ie] = Get_Laplace_transfo(  Eg_virt,  s, E0_u,  Nts[iens], tmax_reco_u-1, prec, SM_TYPE+"_IM",K_IM, ax_u_TO_2, syst_im_u[ie], mult_im_u, l_im_u, MODE_FF, "Ef_"+to_string_with_precision(E0_fact_u,1), TAG_CURR+"Au_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_"+Ens_tags[iens], -1,0, renorm_A, "virtual_FF", cov_ax_u, fake_func,0, fake_func_d ,  Integrate_Up_To_Emax, Emax, beta, 1,0, F_NORM_u, Atr_GEN_u);
	      cout<<"Im HA u, computed"<<endl;
	      IM_HA_sm_d[iens][mu][nu][ixg][isg].distr_list[ie] = Get_Laplace_transfo(  Eg_virt,  s, E0_d_IM,  Nts[iens], tmax_reco_d-1, prec, SM_TYPE+"_IM",K_IM, ax_d_TO_2, syst_im_d[ie], mult_im_d, l_im_d, MODE_FF, "Ef_"+to_string_with_precision(E0_fact_d,1), TAG_CURR+"Ad_mu_"+to_string(mu)+"_nu_"+to_string(nu)+"_"+Ens_tags[iens], -1,0, renorm_A, "virtual_FF", cov_ax_d, fake_func,0, fake_func_d ,  Integrate_Up_To_Emax, Emax, beta, 1,0, F_NORM_d, Atr_GEN_d);
	      cout<<"Im HA d, computed"<<endl;


	      cout<<"Computed ixg: "<<ixg<<" xk: "<<virt_list[ie]<<" Eg: "<<Eg_virt<<" E0_fact_u: "<<E0_fact_u<<", E0_fact_d: "<<E0_fact_d<<" sigma: "<<sigmas[isg]<<" SM_TYPE: "<<SM_TYPE<<" CONS CURRENT: "<<CONS_EM_CURR<<endl<<flush;

	     	     
	    }

	    //print to file
	    //Real part
	    Print_To_File({}, {virt_list, RE_HA_sm_u[iens][mu][nu][ixg][isg].ave(), RE_HA_sm_u[iens][mu][nu][ixg][isg].err(), syst_re_u,  RE_HA_sm_d[iens][mu][nu][ixg][isg].ave(), RE_HA_sm_d[iens][mu][nu][ixg][isg].err(), syst_re_d, (RE_HA_sm_u[iens][mu][nu][ixg][isg]-RE_HA_sm_d[iens][mu][nu][ixg][isg]).ave(), (RE_HA_sm_u[iens][mu][nu][ixg][isg]-RE_HA_sm_d[iens][mu][nu][ixg][isg]).err() }, "../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+"RE_A_alpha_"+to_string_with_precision(beta,2)+"_E0_"+to_string_with_precision(E0_fact,2)+"_SM_"+SM_TYPE+"_ixg_"+to_string(ixg)+"_sigma_"+to_string_with_precision(sigmas[isg],3)+"_mu_"+to_string(mu)+"_nu_"+to_string(nu)+".dat", "", "#xk  u d  u+d");
	    //Imag part
	    Print_To_File({}, {virt_list, IM_HA_sm_u[iens][mu][nu][ixg][isg].ave(), IM_HA_sm_u[iens][mu][nu][ixg][isg].err(), syst_im_u, IM_HA_sm_d[iens][mu][nu][ixg][isg].ave(), IM_HA_sm_d[iens][mu][nu][ixg][isg].err(), syst_im_d, (IM_HA_sm_u[iens][mu][nu][ixg][isg] - IM_HA_sm_d[iens][mu][nu][ixg][isg]).ave(),  (IM_HA_sm_u[iens][mu][nu][ixg][isg] - IM_HA_sm_d[iens][mu][nu][ixg][isg]).err()   }, "../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+"IM_A_alpha_"+to_string_with_precision(beta,2)+"_E0_"+to_string_with_precision(E0_fact,2)+"_SM_"+SM_TYPE+"_ixg_"+to_string(ixg)+"_sigma_"+to_string_with_precision(sigmas[isg],3)+"_mu_"+to_string(mu)+"_nu_"+to_string(nu)+".dat", "", "#xk  u   d   u+d"); 

	    cout<<"sigma: "<<sigmas[isg]<<" computed!"<<endl;
	   
	  }
	  cout<<"done!"<<endl<<flush;
	}
      }
    }
  }


 
  if( !Skip_spectral_reconstruction) {
    cout<<"Storing jackknives distributions for spectral quantities..."<<endl;
    for(int iens=0;iens<Nens;iens++) {
      string TAG_CURR_NEW= ((CONS_EM_CURR==0)?"LOC":"CONS");
      boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR_NEW);
      boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR_NEW+"/jackknives");
      boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR_NEW+"/jackknives/alpha_"+to_string_with_precision(beta,2)+"_E0_"+to_string_with_precision(E0_fact,2)+"_SM_TYPE_"+SM_TYPE);
      for(int ixg=1;ixg<n_xg;ixg++) {
	boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR_NEW+"/jackknives/alpha_"+to_string_with_precision(beta,2)+"_E0_"+to_string_with_precision(E0_fact,2)+"_SM_TYPE_"+SM_TYPE+"/ixg_"+to_string(ixg));
	for(int isg=0; isg < (signed)sigmas.size();isg++) {
	  boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR_NEW+"/jackknives/alpha_"+to_string_with_precision(beta,2)+"_E0_"+to_string_with_precision(E0_fact,2)+"_SM_TYPE_"+SM_TYPE+"/ixg_"+to_string(ixg)+"/sigma_"+to_string_with_precision(sigmas[isg],3));
	  for(int ixk=0;ixk < (signed)virt_list.size();ixk++) {
	    boost::filesystem::create_directory("../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR_NEW+"/jackknives/alpha_"+to_string_with_precision(beta,2)+"_E0_"+to_string_with_precision(E0_fact,2)+"_SM_TYPE_"+SM_TYPE+"/ixg_"+to_string(ixg)+"/sigma_"+to_string_with_precision(sigmas[isg],3)+"/ixk_"+to_string(ixk));
	    if(Reconstruct_vector_part) {
	      for( auto &pair_V:red_mu_nu_pair_V) {
		int mu= pair_V.first;
		int nu= pair_V.second;
		Print_To_File({}, {RE_HV_sm_d[iens][mu][nu][ixg][isg].distr_list[ixk].distr, IM_HV_sm_d[iens][mu][nu][ixg][isg].distr_list[ixk].distr}, "../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+"/jackknives/alpha_"+to_string_with_precision(beta,2)+"_E0_"+to_string_with_precision(E0_fact,2)+"_SM_TYPE_"+SM_TYPE+"/ixg_"+to_string(ixg)+"/sigma_"+to_string_with_precision(sigmas[isg],3)+"/ixk_"+to_string(ixk)+"/Vd_mu_"+to_string(mu)+"_nu_"+to_string(nu)+".jack", "", "");
	       
	      }
	    }
	    if(Reconstruct_axial_part) {
	       for( auto &pair_A:red_mu_nu_pair_A) {
		int mu=	pair_A.first;
		int nu=	pair_A.second;
		Print_To_File({}, {RE_HA_sm_d[iens][mu][nu][ixg][isg].distr_list[ixk].distr, IM_HA_sm_d[iens][mu][nu][ixg][isg].distr_list[ixk].distr},	"../data/ph_emission_3d/"+ph_type_mes+"/FF/"+Ens_tags[iens]+"/"+TAG_CURR+"/jackknives/alpha_"+to_string_with_precision(beta,2)+"_E0_"+to_string_with_precision(E0_fact,2)+"_SM_TYPE_"+SM_TYPE+"/ixg_"+to_string(ixg)+"/sigma_"+to_string_with_precision(sigmas[isg],3)+"/ixk_"+to_string(ixk)+"/Ad_mu_"+to_string(mu)+"_nu_"+to_string(nu)+".jack",	"", "");

	       }
	    }
	  }
	}
      }
    }
    cout<<"Jackknives distribution printed! "<<endl;
  }
   




  


  cout<<"Bye!"<<endl;
  return;
}
