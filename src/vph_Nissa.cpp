#include "../include/vph_Nissa.h"

using namespace std;

const double M2PiPhys=pow(0.135,2);
const double alpha = 1/137.04;
const double e2 = alpha*4.0*M_PI;
const int Nboots= 100;
const bool UseJack=1;
const int nboots=150;
const int Njacks=40;
const double qu = 2.0/3.0; //electric charge of u-type quark
const double qd = -1.0/3.0; //electric charge of d-type quark
const string Meson="Ds";
double Lambda_QCD= 0.3; //300 MeV
bool Is_reph=true;
bool Perform_continuum_extrapolation=true;
int num_xg=11;
Vfloat xg_t_list;
Vfloat xg_to_spline;
Vfloat xg_to_spline_VMD;
Vfloat a_to_print;


void Get_xg_t_list() {

  for(int ixg=1;ixg<num_xg;ixg++) xg_t_list.push_back( ixg*0.1);
  return;
}



void Get_lattice_spacings_to_print() {

  int Nlat=300;
  double sx= 0.08*1.5/(Nlat-1.0); //fm
  for(int a=0;a<Nlat;a++) a_to_print.push_back( sx*a);


}

void Get_xg_to_spline() {

  int Nxgs=300;
  double sx = 0.9/(Nxgs-1.0);
  for(int ixg=0;ixg<Nxgs;ixg++) { xg_to_spline.push_back(0.1+sx*ixg);}
    
  return;
}

void Get_xg_to_spline_VMD() {

  int Nxgs=300;
  double sx = 1.0/(Nxgs-1.0);
  for(int ixg=0;ixg<Nxgs;ixg++) { xg_to_spline_VMD.push_back(0.0+sx*ixg);}
    
  return;
}


void Get_Tmin_Tmax(string W, int &Tmin, int &Tmax, int ixg, string Ens) {


   if(Ens == "cA211a.12.48") {

    if(ixg==0) {
      if(W=="A") {Tmin= 16; Tmax=22;}
      else {Tmin=20;Tmax=32;}

    }
    else if(ixg==1) {
      if(W=="A") {Tmin= 16; Tmax=22;}
      else {Tmin=20;Tmax=32;}

    }
    else if(ixg==2) {
      if(W=="A") {Tmin= 16; Tmax=23;}
      else {Tmin=12;Tmax=25;}

    }
    else if(ixg==3) {
      if(W=="A") {Tmin= 16; Tmax=23;}
      else {Tmin=18;Tmax=28;}
    
    }
    else if(ixg==4) {
      if(W=="A") {Tmin= 17; Tmax=23;}
      else {Tmin=18;Tmax=28;}
    
    }
    else if(ixg==5) {
      if(W=="A") {Tmin= 17; Tmax=23;}
      else {Tmin=18;Tmax=26;}
    
    }
    else if(ixg==6) {
      if(W=="A") {Tmin= 17; Tmax=26;}
      else {Tmin=17;Tmax=26;}
    
    }
    else if(ixg==7) {
      if(W=="A") {Tmin= 17; Tmax=28;}
      else {Tmin=17;Tmax=26;}
    
    }
    else if(ixg==8) {
      if(W=="A") {Tmin= 17; Tmax=29;}
      else {Tmin=15;Tmax=26;}
    
    }
    else if(ixg==9) {
      if(W=="A") {Tmin= 14; Tmax=29;}
      else {Tmin=15;Tmax=26;}
    
    }
    else if(ixg==10) {
      if(W=="A") {Tmin= 14; Tmax=25;}
      else {Tmin=15;Tmax=27;}

    }
    else crash("ixg: "+to_string(ixg)+" does not have an established fit range");

  }


  

   else if(Ens == "cB211b.072.64") {

    if(ixg==0) {
      if(W=="A") {Tmin= 20; Tmax=30;}
      else {Tmin=20;Tmax=30;}

    }
    else if(ixg==1) {
      if(W=="A") {Tmin= 17; Tmax=24;}
      else {Tmin=21;Tmax=45;}

    }
    else if(ixg==2) {
      if(W=="A") {Tmin= 17; Tmax=29;}
      else {Tmin=23;Tmax=45;}

    }
    else if(ixg==3) {
      if(W=="A") {Tmin= 18; Tmax=28;}
      else {Tmin=23;Tmax=45;}
    
    }
    else if(ixg==4) {
      if(W=="A") {Tmin= 18; Tmax=29;}
      else {Tmin=18;Tmax=38;}
    
    }
    else if(ixg==5) {
      if(W=="A") {Tmin= 17; Tmax=30;}
      else {Tmin=18;Tmax=38;}
    
    }
    else if(ixg==6) {
      if(W=="A") {Tmin= 17; Tmax=29;}
      else {Tmin=18;Tmax=38;}
    
    }
    else if(ixg==7) {
      if(W=="A") {Tmin= 14; Tmax=29;}
      else {Tmin=18;Tmax=38;}
    
    }
    else if(ixg==8) {
      if(W=="A") {Tmin= 10; Tmax=29;}
      else {Tmin=12;Tmax=38;}
    
    }
    else if(ixg==9) {
      if(W=="A") {Tmin= 10; Tmax=30;}
      else {Tmin=26;Tmax=39;}
    
    }
    else if(ixg==10) {
      if(W=="A") {Tmin= 10; Tmax=34;}
      else {Tmin=26;Tmax=39;}

    }
    else if(ixg==11) {
      if(W=="A") {Tmin= 20; Tmax=30;}
      else {Tmin=20;Tmax=30;}

    }
    else if(ixg==12) {
      if(W=="A") {Tmin= 20; Tmax=30;}
      else {Tmin=20;Tmax=30;}

    }
    else crash("ixg: "+to_string(ixg)+" does not have an established fit range");

  }


  //############ ENSEMBLE C80  ####################

  else  if(Ens == "cC211a.06.80") {

    if(ixg==0) {
      if(W=="A") {Tmin= 19; Tmax=29;}
      else {Tmin=28;Tmax=55;}

    }
    else if(ixg==1) {
      if(W=="A") {Tmin= 19; Tmax=29;}
      else {Tmin=28;Tmax=55;}

    }
    else if(ixg==2) {
      if(W=="A") {Tmin= 21; Tmax=29;}
      else {Tmin=28;Tmax=55;}

    }
    else if(ixg==3) {
      if(W=="A") {Tmin= 21; Tmax=29;}
      else {Tmin=28;Tmax=55;}
    
    }
    else if(ixg==4) {
      if(W=="A") {Tmin= 21; Tmax=30;}
      else {Tmin=28;Tmax=55;}
    
    }
    else if(ixg==5) {
      if(W=="A") {Tmin= 21; Tmax=30;}
      else {Tmin=25;Tmax=55;}
    
    }
    else if(ixg==6) {
      if(W=="A") {Tmin= 21; Tmax=30;}
      else {Tmin=28;Tmax=47;}
    
    }
    else if(ixg==7) {
      if(W=="A") {Tmin= 21; Tmax=31;}
      else {Tmin=16;Tmax=40;}
    
    }
    else if(ixg==8) {
      if(W=="A") {Tmin= 21; Tmax=37;}
      else {Tmin=15;Tmax=47;}
    
    }
    else if(ixg==9) {
      if(W=="A") {Tmin= 35; Tmax=47;}
      else {Tmin=24;Tmax=46;}
    
    }
    else if(ixg==10) {
      if(W=="A") {Tmin= 35; Tmax=47;}
      else {Tmin=25;Tmax=46;}

    }

    else crash("ixg: "+to_string(ixg)+" does not have an established fit range");

  }

  else if(Ens=="cD211a.054.96") {

   


    if(ixg==0) {
      if(W=="A") {Tmin= 21; Tmax=28;}
      else {Tmin=28;Tmax=55;}

    }
    else if(ixg==1) {
      if(W=="A") {Tmin= 19; Tmax=35;}
      else {Tmin=29;Tmax=45;}

    }
    else if(ixg==2) {
      if(W=="A") {Tmin= 19; Tmax=35;}
      else {Tmin=29;Tmax=45;}

    }
    else if(ixg==3) {
      if(W=="A") {Tmin= 20; Tmax=35;}
      else {Tmin=29;Tmax=45;}
    
    }
    else if(ixg==4) {
      if(W=="A") {Tmin= 20; Tmax=35;}
      else {Tmin=29;Tmax=49;}
    
    }
    else if(ixg==5) {
      if(W=="A") {Tmin= 20; Tmax=36;}
      else {Tmin=30;Tmax=49;}
    
    }
    else if(ixg==6) {
      if(W=="A") {Tmin= 20; Tmax=35;}
      else {Tmin=30;Tmax=50;}
    
    }
    else if(ixg==7) {
      if(W=="A") {Tmin= 21; Tmax=37;}
      else {Tmin=29;Tmax=49;}
    
    }
    else if(ixg==8) {
      if(W=="A") {Tmin= 19; Tmax=45;}
      else {Tmin=14;Tmax=42;}
    
    }
    else if(ixg==9) {
      if(W=="A") {Tmin= 20; Tmax=44;}
      else {Tmin=33;Tmax=67;}
    
    }
    else if(ixg==10) {
      if(W=="A") {Tmin= 29; Tmax=54;}
      else {Tmin=45;Tmax=67;}

    }

    else crash("ixg: "+to_string(ixg)+" does not have an established fit range");


  }

  else crash("Ensemble: "+Ens+" does not have definite time intervals for FV and FA");


  

  return; 
}




void Compute_form_factors_Nissa() {

  

  Get_xg_t_list();
  Get_xg_to_spline();
  Get_xg_to_spline_VMD();
  Get_lattice_spacings_to_print();
  int size_mu_nu= Is_reph?2:4;
  string ph_type= Is_reph?"rph":"vph";

  //create directories
  boost::filesystem::create_directory("../data/ph_emission");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type);
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson);
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/C");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/H");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/mass");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/decay_const");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF/per_kin");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF_u");
  boost::filesystem::create_directory("../data/ph_emission/"+ph_type+"/"+Meson+"/FF_d");
   
  string ph_type_mes=ph_type+"/"+Meson;
  
  
  //axial
  vector<vector<vector<data_t>>> C_A_Fu_data(size_mu_nu), C_A_Fd_data(size_mu_nu), C_A_Bu_data(size_mu_nu), C_A_Bd_data(size_mu_nu);
  //vector
  vector<vector<vector<data_t>>> C_V_Fu_data(size_mu_nu), C_V_Fd_data(size_mu_nu), C_V_Bu_data(size_mu_nu), C_V_Bd_data(size_mu_nu);

 

  


  //2pts
  data_t data_2pts, data_2pts_SM;
  data_t data_2pts_V1, data_2pts_V2, data_2pts_V3;

  for(int mu=0;mu<size_mu_nu;mu++) {

    C_A_Fu_data[mu].resize(size_mu_nu);
    C_A_Fd_data[mu].resize(size_mu_nu);
    C_A_Bu_data[mu].resize(size_mu_nu);
    C_A_Bd_data[mu].resize(size_mu_nu);
    C_V_Fu_data[mu].resize(size_mu_nu);
    C_V_Fd_data[mu].resize(size_mu_nu);
    C_V_Bu_data[mu].resize(size_mu_nu);
    C_V_Bd_data[mu].resize(size_mu_nu);

    for(int nu=0;nu<size_mu_nu;nu++) {

      C_A_Fu_data[mu][nu].resize(num_xg);
      C_A_Fd_data[mu][nu].resize(num_xg);
      C_A_Bu_data[mu][nu].resize(num_xg);
      C_A_Bd_data[mu][nu].resize(num_xg);
      C_V_Fu_data[mu][nu].resize(num_xg);
      C_V_Fd_data[mu][nu].resize(num_xg);
      C_V_Bu_data[mu][nu].resize(num_xg);
      C_V_Bd_data[mu][nu].resize(num_xg);


    }

  }


  
  int off_i = (Is_reph?1:0);
  

 
  //Read data

  //2pts function

  data_2pts.Read("../new_vph_gpu_data_w_123", "mes_contr_2pts_3", "P5P5");
  data_2pts_SM.Read("../new_vph_gpu_data_w_123", "mes_contr_2pts_SM_3", "P5P5");
  data_2pts_V1.Read("../new_vph_gpu_data_w_123", "mes_contr_2pts_3", "V1V1");
  data_2pts_V2.Read("../new_vph_gpu_data_w_123", "mes_contr_2pts_3", "V2V2");
  data_2pts_V3.Read("../new_vph_gpu_data_w_123", "mes_contr_2pts_3", "V3V3");
  

  //read data
  for(int ixg=0;ixg<num_xg;ixg++) {


    for(int mu=0;mu<size_mu_nu;mu++) {

      for(int nu=0;nu<size_mu_nu;nu++) {

	//vector
	//Fu
	C_V_Fu_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123", "C_mu_"+to_string(mu+off_i)+"_FF_u_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5");
	//Fd
	C_V_Fd_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123", "C_mu_"+to_string(mu+off_i)+"_FF_d_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5");
	//Bu
	C_V_Bu_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123", "C_mu_"+to_string(mu+off_i)+"_BB_u_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5");
	//Bd
	C_V_Bd_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123", "C_mu_"+to_string(mu+off_i)+"_BB_d_ixg_"+to_string(ixg), "V"+to_string(nu+off_i)+"P5");

	//axial
	C_A_Fu_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123", "C_mu_"+to_string(mu+off_i)+"_FF_u_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5");
	//Fd
	C_A_Fd_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123", "C_mu_"+to_string(mu+off_i)+"_FF_d_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5");
	//Bu
	C_A_Bu_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123", "C_mu_"+to_string(mu+off_i)+"_BB_u_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5");
	//Bd
	C_A_Bd_data[mu][nu][ixg].Read("../new_vph_gpu_data_w_123", "C_mu_"+to_string(mu+off_i)+"_BB_d_ixg_"+to_string(ixg), "A"+to_string(nu+off_i)+"P5");

      }
    }
  }


  int Nens = data_2pts.size;
  GaussianMersenne GM(543543);

  //vectors to store FV and FA for all ensembles and values of gamma
  vector<distr_t_list> FA_per_ens(Nens), FV_per_ens(Nens);
  vector<distr_t_list> FA_per_kin(num_xg-1), FV_per_kin(num_xg-1); //num_xg -1 is to exclude xg=0 which is undefined
  distr_t_list a_distr_list(UseJack);
  //M_Ds, F_Ds
  distr_t_list MP_list(UseJack);
  distr_t_list FP_list(UseJack);
  distr_t_list MP_ov_FP_list(UseJack);
  //resize
  
  
  //define lambda function to combine FF and BB

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

 

  for(int iens=0;iens<Nens;iens++) {

    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/C/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/H/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results");
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/mass/"+data_2pts.Tag[iens]);
    boost::filesystem::create_directory("../data/ph_emission/"+ph_type_mes+"/decay_const/"+data_2pts.Tag[iens]);
    
    cout<<"Analyzing ensemble: "<<data_2pts.Tag[iens]<<endl;

    //Lattice info
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(data_2pts.Tag[iens]);
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    Corr.Nt = data_2pts.nrows[iens];
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;

    
    //read theta values and loop over them
    Vfloat thetas, masses_u, masses_d, virts;

    thetas= Read_From_File("../new_vph_gpu_data_w_123/"+data_2pts.Tag[iens]+"/pars_list.dat", 1 , 5);
    virts=  Read_From_File("../new_vph_gpu_data_w_123/"+data_2pts.Tag[iens]+"/pars_list.dat", 2 , 5);
    masses_u= Read_From_File("../new_vph_gpu_data_w_123/"+data_2pts.Tag[iens]+"/pars_list.dat", 3 , 5);
    masses_d= Read_From_File("../new_vph_gpu_data_w_123/"+data_2pts.Tag[iens]+"/pars_list.dat", 4 , 5);

    cout<<"pars_list.dat: Read!"<<endl;

    if((signed)thetas.size() != num_xg) crash("Number of rows in pars_list.dat does not match num_xg"); 


    //RCs
    distr_t Za, Zv, a_distr;
    if(data_2pts.Tag[iens].substr(1,1)=="A") { Za= ZA_A; Zv=ZV_A;a_distr=a_A;}
    else if(data_2pts.Tag[iens].substr(1,1)=="B") { Za= ZA_B; Zv=ZV_B;a_distr=a_B;}
    else if(data_2pts.Tag[iens].substr(1,1)=="C") { Za= ZA_C; Zv=ZV_C; a_distr=a_C;}
    else if(data_2pts.Tag[iens].substr(1,1)=="D") { Za= ZA_D; Zv=ZV_D;a_distr=a_D;}
    else crash("Ensemble: "+data_2pts.Tag[iens]+" not recognised");


    a_distr_list.distr_list.push_back(a_distr);


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

    distr_t_list pt2_distr= Corr.corr_t(data_2pts.col(0)[iens], "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt.dat");
    distr_t_list eff_mass = Corr.effective_mass_t(pt2_distr, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass.dat");
    distr_t_list fp_distr= Corr.decay_constant_t( pow( mu+md,2)*pt2_distr, "../data/ph_emission/"+ph_type_mes+"/"+"decay_const/"+data_2pts.Tag[iens]+"/decay_const.dat");
    distr_t_list pt2_vector_distr= (1.0/3.0)*Corr.corr_t(summ_master(data_2pts_V1.col(0)[iens], data_2pts_V2.col(0)[iens], data_2pts_V3.col(0)[iens]), "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_V.dat");
    distr_t_list eff_mass_V = Corr.effective_mass_t(pt2_vector_distr, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_V.dat");
    distr_t M_P=Corr.Fit_distr(eff_mass);
    distr_t F_P=Corr.Fit_distr(fp_distr);
    MP_list.distr_list.push_back(M_P/a_distr);
    FP_list.distr_list.push_back(F_P/a_distr);
    MP_ov_FP_list.distr_list.push_back( M_P/F_P);
    
    //set time interval for eff_mass_fit SM
    if(data_2pts.Tag[iens].substr(1,1) =="A") {Corr.Tmin=24; Corr.Tmax=35;}
    else if(data_2pts.Tag[iens].substr(1,1) =="B") {Corr.Tmin=20; Corr.Tmax=36;}
    else if(data_2pts.Tag[iens].substr(1,1) == "C")  {Corr.Tmin=33;Corr.Tmax=51;}
    else if(data_2pts.Tag[iens].substr(1,1) == "D")  {Corr.Tmin=41;Corr.Tmax=54;}
    else crash("In fixing [Tmin, Tmax] for smeared MP, Ensemble: "+data_2pts.Tag[iens]+" not recognized");
    distr_t_list pt2_distr_SM= Corr.corr_t(data_2pts_SM.col(0)[iens], "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/corr_2pt_SM.dat");
    distr_t_list eff_mass_SM = Corr.effective_mass_t(pt2_distr_SM, "../data/ph_emission/"+ph_type_mes+"/"+"mass/"+data_2pts.Tag[iens]+"/eff_mass_SM.dat");
    distr_t M_P_SM = Corr.Fit_distr(eff_mass_SM);


    cout<<"aM_P: "<<M_P.ave()<<" +- "<<M_P.err()<<" -> "<< (M_P/(a_distr)).ave()<<" +- "<<(M_P/a_distr).err()<<" GeV"<<endl;
    cout<<"aF_P: "<<F_P.ave()<<" +- "<<F_P.err()<<" -> "<< (F_P/(a_distr)).ave()<<" +- "<<(F_P/a_distr).err()<<" GeV"<<endl;
    cout<<"MP/FP: "<<(M_P/F_P).ave()<<" +- "<<(M_P/F_P).err()<<endl;
    cout<<"aM_P(SM): "<<M_P_SM.ave()<<" +- "<<M_P_SM.err()<<" -> "<< (M_P_SM/(a_distr)).ave()<<" +- "<<(M_P_SM/a_distr).err()<<" GeV"<<endl;
    cout<<"MP/FP: "<<(M_P/F_P).ave()<<" +- "<<(M_P/F_P).err()<<endl;


    //define meson mass exponential to be removed
    auto EXP_MES_FUNC = [&] (double a, double b, double c) { return (b<c/2)?1.0/(exp(-a*b)):1.0/(exp(-a*(c-b)));};

    distr_t_list EXP_MES= distr_t_list::f_of_distr( EXP_MES_FUNC, M_P, Corr.Nt);

   

    vector<vector<vector<distr_t_list>>> Ax_glb, Vec_glb, Ax_u_glb, Ax_d_glb, Vec_u_glb, Vec_d_glb;
    distr_t_list FV(UseJack), FA(UseJack), FV_u(UseJack), FV_d(UseJack), FA_u(UseJack), FA_d(UseJack);


    distr_t_list xg_list(UseJack);

    Vfloat ax_Tmin, ax_Tmax;
    Vfloat vec_Tmin, vec_Tmax;
    

    for(int ixg=0;ixg<num_xg;ixg++) {

      vector<vector<distr_t_list>> Ax_tens(size_mu_nu);
      vector<vector<distr_t_list>> Ax_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Ax_tens_d(size_mu_nu);
      vector<vector<distr_t_list>> Vec_tens(size_mu_nu);
      vector<vector<distr_t_list>> Vec_tens_u(size_mu_nu);
      vector<vector<distr_t_list>> Vec_tens_d(size_mu_nu);
      


      //get xg, Eg, kz from thetas

      double theta=thetas[ixg];

      pt3_momenta pt3_mom(0.0, 0.0, thetas[ixg]/2.0, masses_u[ixg], masses_d[ixg], virts[ixg], L_info.L, L_info.T);

      double Eg= pt3_mom.Egamma();
      distr_t xg= pt3_mom.x_gamma(M_P);
      if(xg.ave() > 1e-10 ) xg_list.distr_list.push_back(xg);
      double kz = pt3_mom.k()[2];
    
   

      //define photon exponential to be removed
      Vfloat EXP_PH(Corr.Nt,0.0);
      for(int t=0; t < Corr.Nt;t++) EXP_PH[t] = exp( Eg*abs( Corr.Nt/2 - t));
      
     


      cout<<"##### Considering kinematic with......"<<endl;
      cout<<"Eg: "<<Eg<<endl;
      cout<<"xg: "<<xg.ave()<<" +- "<<xg.err()<<endl;
      cout<<"thz: "<<theta<<endl;
      cout<<"kz: "<<kz<<endl;

    
    

      //define vectors to combine FF and BB
      Vfloat th_FF(Corr.Nt);
      Vfloat th_BB(Corr.Nt);
      for(int t=0;t<Corr.Nt;t++) {
	if(t<Corr.Nt/2) {th_FF[t]=1.0; th_BB[t] = 0.0;}
	else {th_BB[t] = 1.0; th_FF[t] = 0.0;}
      }

      //loop over mu and nu
      for(int mu=0;mu<size_mu_nu;mu++) {
	for(int nu=0;nu<size_mu_nu;nu++) {

	  cout<<"Analyzing mu: "<<mu+off_i<<" nu: "<<nu+off_i<<endl;
	  int Im_Re;
	  double parity;

	  //vector
	
	  Corr.Reflection_sign = -1;
	  Im_Re=1;
	  parity=1.0;

	  Corr.Perform_Nt_t_average = 0;
	  distr_t_list vec_F_u = qu*Corr.corr_t(C_V_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_B_u = qu*Corr.corr_t(C_V_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_F_d = qd*Corr.corr_t(C_V_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list vec_B_d = qd*Corr.corr_t(C_V_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	
	  distr_t_list vec_u = th_FF*vec_F_u + th_BB*vec_B_u;
	  distr_t_list vec_d = th_FF*vec_F_d + th_BB*vec_B_d;
	  distr_t_list vec = vec_u + vec_d;
	  distr_t_list vec_symm= vec;
	  distr_t_list vec_u_symm= vec_u;
	  distr_t_list vec_d_symm= vec_d;
	  //symmetrize vec
	  for(int t=0; t<Corr.Nt;t++) {
	    vec_symm.distr_list[t] = 0.5*(vec.distr_list[t] + Corr.Reflection_sign*vec.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    vec_u_symm.distr_list[t] = 0.5*(vec_u.distr_list[t] + Corr.Reflection_sign*vec_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    vec_d_symm.distr_list[t] = 0.5*(vec_d.distr_list[t] + Corr.Reflection_sign*vec_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
	  }
	  Corr.Perform_Nt_t_average=1;
	
	

	  //axial
	  if( (mu==0 || nu==0) && (mu != 0 || nu != 0)) {Im_Re=1; Corr.Reflection_sign=-1; parity=1.0;}
	  else { Im_Re=0; Corr.Reflection_sign=1; parity=1.0;}

	  Corr.Perform_Nt_t_average=0;
	  distr_t_list ax_F_u = qu*Corr.corr_t(C_A_Fu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_B_u = qu*Corr.corr_t(C_A_Bu_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_F_d = qd*Corr.corr_t(C_A_Fd_data[mu][nu][ixg].col(Im_Re)[iens],"");
	  distr_t_list ax_B_d = qd*Corr.corr_t(C_A_Bd_data[mu][nu][ixg].col(Im_Re)[iens],"");

	
	  distr_t_list ax_u = parity*(ax_F_u*th_FF  + ax_B_u*th_BB);
	  distr_t_list ax_d = parity*(ax_F_d*th_FF  + ax_B_d*th_BB);
	  distr_t_list ax = ax_u -ax_d;
	  distr_t_list ax_symm=ax;
	  distr_t_list ax_u_symm= ax_u;
	  distr_t_list ax_d_symm= ax_d;
	  //symmetrize ax
	  for(int t=0; t<Corr.Nt;t++) {
	    ax_symm.distr_list[t] = 0.5*(ax.distr_list[t] + Corr.Reflection_sign*ax.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_u_symm.distr_list[t] = 0.5*(ax_u.distr_list[t] + Corr.Reflection_sign*ax_u.distr_list[( Corr.Nt -t)%Corr.Nt]);
	    ax_d_symm.distr_list[t] = 0.5*(ax_d.distr_list[t] + Corr.Reflection_sign*ax_d.distr_list[( Corr.Nt -t)%Corr.Nt]);
	  }
	  Corr.Perform_Nt_t_average=1;

	  //restore standard reflection sign
	  Corr.Reflection_sign=1;


	  //get H tensor
	  distr_t_list HA_u= ax_u*EXP_MES*EXP_PH;
	  distr_t_list HA_d= ax_d*EXP_MES*EXP_PH;
	  distr_t_list HV_u= vec_u*EXP_MES*EXP_PH;
	  distr_t_list HV_d= vec_d*EXP_MES*EXP_PH;
	  distr_t_list HA_u_symm= ax_u_symm*EXP_MES*EXP_PH;
	  distr_t_list HA_d_symm= ax_d_symm*EXP_MES*EXP_PH;
	  distr_t_list HV_u_symm= vec_u_symm*EXP_MES*EXP_PH;
	  distr_t_list HV_d_symm= vec_d_symm*EXP_MES*EXP_PH;
	  distr_t_list HA_symm= ax_symm*EXP_MES*EXP_PH;
	  distr_t_list HV_symm= vec_symm*EXP_MES*EXP_PH;



	  //########################   C TENSOR ###########################


	  //print forward and backward contribution to C
	  Print_To_File({}, {ax_F_u.ave(), ax_F_u.err(), ax_F_d.ave(), ax_F_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_FF_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_F_u    ax_F_d");
	  Print_To_File({}, {ax_B_u.ave(), ax_B_u.err(), ax_B_d.ave(), ax_B_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_BB_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_B_u    ax_B_d");
	  Print_To_File({}, {vec_F_u.ave(), vec_F_u.err(), vec_F_d.ave(), vec_F_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_FF_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_F_u    vec_F_d");
	  Print_To_File({}, {vec_B_u.ave(), vec_B_u.err(), vec_B_d.ave(), vec_B_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_BB_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_B_u    vec_B_d");
	

	  //single contributions to C
	  Print_To_File({}, {ax_u.ave() ,ax_u.err(), ax_d.ave(), ax_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_u    ax_d");
	  Print_To_File({}, {vec_u.ave() ,vec_u.err(), vec_d.ave(), vec_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_u    vec_d");

	  //single contributions to C with symmetrization
	  Print_To_File({}, {ax_u_symm.ave() ,ax_u_symm.err(), ax_d_symm.ave(), ax_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   ax_u    ax_d");
	  Print_To_File({}, {vec_u_symm.ave() ,vec_u_symm.err(), vec_d_symm.ave(), vec_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/quark_contr_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   vec_u    vec_d");

	  //total contribution without symmetrization to C

	  Print_To_File({}, {ax.ave(), ax.err(), vec.ave(), vec.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/wo_T_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    A     V");

	  //total contribution to C
	  Print_To_File({}, {ax_symm.ave(), ax_symm.err(), vec_symm.ave(), vec_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"C/"+data_2pts.Tag[iens]+"/mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    A     V");


	  //###############################################################


	  //########################   H TENSOR ###########################


	  //single contributions to H without symmetrization
	   Print_To_File({}, {HA_u.ave() ,HA_u.err(), HA_d.ave(), HA_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/quark_contr_A_wo_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HA_u    HA_d");
	  Print_To_File({}, {HV_u.ave() ,HV_u.err(), HV_d.ave(), HV_d.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/quark_contr_V_wo_symm_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HV_u    HV_d");
	  

	  
	  //single contributions to H with symmetrization
	  Print_To_File({}, {HA_u_symm.ave() ,HA_u_symm.err(), HA_d_symm.ave(), HA_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/quark_contr_A_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HA_u    HA_d");
	  Print_To_File({}, {HV_u_symm.ave() ,HV_u_symm.err(), HV_d_symm.ave(), HV_d_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/quark_contr_V_mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "#   HV_u    HV_d");

	  //total contribution to H
	  Print_To_File({}, {HA_symm.ave(), HA_symm.err(), HV_symm.ave(), HV_symm.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"H/"+data_2pts.Tag[iens]+"/mu_"+to_string(mu+off_i)+"_nu_"+to_string(nu+off_i)+"_xg_"+to_string_with_precision(xg.ave(),2)+".dat" , "", "#  t    HA     HV");

	  //################################################################

	  //push back
	  Ax_tens[mu].push_back(ax_symm);
	  Ax_tens_u[mu].push_back(ax_u_symm);
	  Ax_tens_d[mu].push_back(ax_d_symm);
	  Vec_tens[mu].push_back(vec_symm);
	  Vec_tens_u[mu].push_back(vec_u_symm);
	  Vec_tens_d[mu].push_back(vec_d_symm);
	  
	

	}
      }

      //push_back Ax_tens and Vec_tens
      Ax_glb.push_back(Ax_tens);
      Ax_u_glb.push_back(Ax_tens_u);
      Ax_d_glb.push_back(Ax_tens_d);
      
      Vec_glb.push_back(Vec_tens);
      Vec_u_glb.push_back(Vec_tens_u);
      Vec_d_glb.push_back(Vec_tens_d);

      //Compute FV and FA
      distr_t_list FA0_distr= 0.5*(Ax_glb[0][1-off_i][1-off_i] + Ax_glb[0][2-off_i][2-off_i]);
      distr_t_list FA_distr = (0.5*(Ax_tens[1-off_i][1-off_i] + Ax_tens[2-off_i][2-off_i])*EXP_PH - FA0_distr)*(1.0/Eg)*F_P/FA0_distr;
      distr_t_list FV_distr = -0.5*(Za/Zv)*(F_P/(-1.0*FA0_distr))*( Vec_tens[1-off_i][2-off_i] - Vec_tens[2-off_i][1-off_i])*EXP_PH/kz;
      distr_t_list FV_sub_distr= -0.5*(Za/Zv)*(F_P/(-1.0*FA0_distr))*( (Vec_tens[1-off_i][2-off_i] - Vec_tens[2-off_i][1-off_i])*EXP_PH - (Vec_glb[0][1-off_i][2-off_i] + Vec_glb[0][2-off_i][1-off_i])*1.0)/kz;
      //compute FV and FA (up-component)
      distr_t_list FA0_u_distr= 0.5*(Ax_u_glb[0][1-off_i][1-off_i] + Ax_u_glb[0][2-off_i][2-off_i]);
      distr_t_list FA_u_distr = (0.5*(Ax_tens_u[1-off_i][1-off_i] + Ax_tens_u[2-off_i][2-off_i])*EXP_PH - FA0_u_distr)*(1.0/Eg)*F_P/FA0_distr;
      distr_t_list FV_u_distr = -0.5*(Za/Zv)*(F_P/(-1.0*FA0_distr))*( Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH/kz;
      distr_t_list FV_u_sub_distr= -0.5*(Za/Zv)*(F_P/(-1.0*FA0_distr))*( (Vec_tens_u[1-off_i][2-off_i] - Vec_tens_u[2-off_i][1-off_i])*EXP_PH - (Vec_u_glb[0][1-off_i][2-off_i] + Vec_u_glb[0][2-off_i][1-off_i])*1.0)/kz;

      //compute FV and FA (d-component)
      distr_t_list FA0_d_distr= 0.5*(Ax_d_glb[0][1-off_i][1-off_i] + Ax_d_glb[0][2-off_i][2-off_i]);
      distr_t_list FA_d_distr = (0.5*(Ax_tens_d[1-off_i][1-off_i] + Ax_tens_d[2-off_i][2-off_i])*EXP_PH - FA0_d_distr)*(1.0/Eg)*F_P/FA0_distr;
      distr_t_list FV_d_distr = -0.5*(Za/Zv)*(F_P/(-1.0*FA0_distr))*( Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH/kz;
      distr_t_list FV_d_sub_distr= -0.5*(Za/Zv)*(F_P/(-1.0*FA0_distr))*( (Vec_tens_d[1-off_i][2-off_i] - Vec_tens_d[2-off_i][1-off_i])*EXP_PH - (Vec_d_glb[0][1-off_i][2-off_i] + Vec_d_glb[0][2-off_i][1-off_i])*1.0)/kz;


      //Print FV and FAt
      Print_To_File({}, {FA_distr.ave(), FA_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FA_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_distr.ave(), FV_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FV_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_sub_distr.ave(), FV_sub_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF/"+data_2pts.Tag[iens]+"/FV_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");

      //Print FV and FA (up-component)
      Print_To_File({}, {FA_u_distr.ave(), FA_u_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FA_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_u_distr.ave(), FV_u_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FV_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {FV_u_sub_distr.ave(), FV_u_sub_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_u/"+data_2pts.Tag[iens]+"/FV_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      
      //Print FV and FA (d-component)
      Print_To_File({}, {FA_d_distr.ave(), FA_d_distr.err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FA_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {(-1.0*FV_d_distr).ave(), (-1.0*FV_d_distr).err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FV_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");
      Print_To_File({}, {(-1.0*FV_d_sub_distr).ave(), (-1.0*FV_d_sub_distr).err()}, "../data/ph_emission/"+ph_type_mes+"/"+"FF_d/"+data_2pts.Tag[iens]+"/FV_sub_xg_"+to_string_with_precision(xg.ave(),2)+".dat", "", "");

      //fit the form factors

      //set interval depending on xg
      int Tmin_V, Tmax_V, Tmin_A, Tmax_A;
      Get_Tmin_Tmax("A", Tmin_A, Tmax_A, ixg, data_2pts.Tag[iens]);
      Get_Tmin_Tmax("V", Tmin_V, Tmax_V, ixg, data_2pts.Tag[iens]);

    


      //axial
      Corr.Tmin= Tmin_A; Corr.Tmax= Tmax_A;
      if(xg.ave() > 1e-10) { //exclude xg=0
      FA.distr_list.push_back( Corr.Fit_distr(FA_distr));
      FA_u.distr_list.push_back(Corr.Fit_distr(FA_u_distr));
      FA_d.distr_list.push_back(Corr.Fit_distr(FA_d_distr));
      FA_per_kin[ixg-1].distr_list.push_back(Corr.Fit_distr(FA_distr));
      ax_Tmin.push_back(Tmin_A);
      ax_Tmax.push_back(Tmax_A);
      }
      //vector
      if(xg.ave() > 1e-10) { //exclude xg=0
      Corr.Tmin= Tmin_V; Corr.Tmax= Tmax_V;
      FV.distr_list.push_back( Corr.Fit_distr(FV_distr));
      FV_u.distr_list.push_back( Corr.Fit_distr(FV_u_distr));
      FV_d.distr_list.push_back( Corr.Fit_distr(FV_d_distr));
      FV_per_kin[ixg-1].distr_list.push_back(Corr.Fit_distr(FV_distr));
      vec_Tmin.push_back(Tmin_V);
      vec_Tmax.push_back(Tmax_V);
      }
      
      
    
      
    }

    FA_per_ens[iens] = FA;
    FV_per_ens[iens] = FV;
    


    //Print fitted form factors
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FV.ave(), FV.err(), vec_Tmin, vec_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FV.dat", "", "#xg Dxg  FV DFV Tmin Tmax");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA.ave(), FA.err(), ax_Tmin, ax_Tmax}, "../data/ph_emission/"+ph_type_mes+"/FF/"+data_2pts.Tag[iens]+"/fit_results/FA.dat", "", "#xg Dxg  FA DFA Tmin Tmax");

    Print_To_File({}, {xg_list.ave(), xg_list.err(), FV_u.ave(), FV_u.err()}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FV.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_u.ave(), FA_u.err()}, "../data/ph_emission/"+ph_type_mes+"/FF_u/"+data_2pts.Tag[iens]+"/fit_results/FA.dat", "", "#xg Dxg  FA DFA");

    Print_To_File({}, {xg_list.ave(), xg_list.err(), (-1.0*FV_d).ave(), (-1.0*FV_d).err()}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FV.dat", "", "#xg Dxg  FV DFV");
    Print_To_File({}, {xg_list.ave(), xg_list.err(), FA_d.ave(), FA_d.err()}, "../data/ph_emission/"+ph_type_mes+"/FF_d/"+data_2pts.Tag[iens]+"/fit_results/FA.dat", "", "#xg Dxg  FA DFA");


  }


  


  if(num_xg==1) crash("Exiting...No interpolation nor continuum extrapolation to perform, only xg available is zero");


  //Print all ensembles for given xg
  for(int ixg=1;ixg<num_xg;ixg++) {
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FA_per_kin[ixg-1].ave(), FA_per_kin[ixg-1].err()}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FA_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err");
    Print_To_File({}, {(a_distr_list/fmTGeV).ave(), FV_per_kin[ixg-1].ave(), FV_per_kin[ixg-1].err()}, "../data/ph_emission/"+ph_type_mes+"/FF/per_kin/FV_ixg_"+to_string(ixg)+".dat", "", "#a[fm]  FA FA_err");
  }

  //################################################################################################

  //interpolate form factors for each ensemble

  vector<vector<boost::math::interpolators::cardinal_cubic_b_spline<double>>> FA_interpol_jacks(Nens);
  vector<vector<boost::math::interpolators::cardinal_cubic_b_spline<double>>> FV_interpol_jacks(Nens);

  for(int iens=0;iens<Nens;iens++) {

    for(int ijack=0;ijack<Njacks;ijack++) {

      Vfloat FA_jacks, FV_jacks;
      for(int ixg=1;ixg<num_xg;ixg++) { FA_jacks.push_back( FA_per_ens[iens].distr_list[ixg-1].distr[ijack]); FV_jacks.push_back( FV_per_ens[iens].distr_list[ixg-1].distr[ijack]);}
      FA_interpol_jacks[iens].emplace_back( FA_jacks.begin(), FA_jacks.end(), 0.1, 0.1);
      FV_interpol_jacks[iens].emplace_back( FV_jacks.begin(), FV_jacks.end(), 0.1, 0.1);
    }
  }

  auto FA_interpol_distr = [&FA_interpol_jacks](double xg, int iens) -> distr_t {

			     distr_t return_distr(UseJack);

			     for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FA_interpol_jacks[iens][ijack](xg));}

			     return return_distr;

			   };

  auto FV_interpol_distr = [&FV_interpol_jacks](double xg, int iens) -> distr_t {

			     distr_t return_distr(UseJack);

			     for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FV_interpol_jacks[iens][ijack](xg));}

			     return return_distr;

			     } ;



  for(int iens=0;iens<Nens;iens++) {

    distr_t_list FA_interpol_to_print_distr(UseJack);
    distr_t_list FV_interpol_to_print_distr(UseJack);
    for(auto &X: xg_to_spline) {
      FA_interpol_to_print_distr.distr_list.push_back( FA_interpol_distr(X, iens));
      FV_interpol_to_print_distr.distr_list.push_back( FV_interpol_distr(X, iens));
    }

    Print_To_File({}, {xg_to_spline, FA_interpol_to_print_distr.ave(), FA_interpol_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/"+data_2pts.Tag[iens]+"/FA_interpol.dat", "", "#xg FA FA_err");
    Print_To_File({}, {xg_to_spline, FV_interpol_to_print_distr.ave(), FV_interpol_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/"+data_2pts.Tag[iens]+"/FV_interpol.dat", "", "#xg FV FV_err");


  }


  //Print MP, FP, MP_ov_FP
  Print_To_File({}, {a_distr_list.ave(), MP_list.ave(), MP_list.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/mass/masses.list", "", "#a MP MP_err");
  Print_To_File({}, {a_distr_list.ave(), FP_list.ave(), FP_list.err(), MP_ov_FP_list.ave(), MP_ov_FP_list.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/decay_const/fP.list", "", "#a FP FP_err MP/FP MP/FP_err  ");


  //#################################################################################################

  //continuum extrapolation

  
 

  if(Perform_continuum_extrapolation) {
  
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FA_cont_interpol_jacks;
    vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> FV_cont_interpol_jacks;


  class ipar_FF_Nissa {
  public:
    ipar_FF_Nissa() : FF(0.0), FF_err(0.0) {}
    double FF, FF_err, a;
  };
  
  class fpar_FF_Nissa {
  public:
    fpar_FF_Nissa() {}
    fpar_FF_Nissa(const Vfloat &par) {
      if((signed)par.size() != 2) crash("In class fpar_FF_Nissa  class constructor Vfloat par has size != 2");
      F0=par[0];
      D1=par[1];
    }
    double F0, D1;
  };

  
  //init bootstrap fit
  bootstrap_fit<fpar_FF_Nissa,ipar_FF_Nissa> bf_FF(Njacks);
  bf_FF.set_warmup_lev(0); //sets warmup
  bf_FF.Set_number_of_measurements(Nens);
  bf_FF.Set_verbosity(1);
  bf_FF.Add_par("F0", 0.07, 0.001);
  bf_FF.Add_par("D1", 1.0, 0.1);
  //fit on mean values to get ch2
  bootstrap_fit<fpar_FF_Nissa,ipar_FF_Nissa> bf_FF_ch2(1);
  bf_FF_ch2.set_warmup_lev(0); //sets warmup
  bf_FF_ch2.Set_number_of_measurements(Nens);
  bf_FF_ch2.Set_verbosity(1);
  bf_FF_ch2.Add_par("F0", 0.07, 0.001);
  bf_FF_ch2.Add_par("D1", 1.0, 0.1);
  

  //ansatz
  bf_FF.ansatz=  [ ](const fpar_FF_Nissa &p, const ipar_FF_Nissa &ip) {

		  
		   return p.F0 + p.D1*pow(ip.a*Lambda_QCD,2);
		 };
  bf_FF.measurement=  [ ](const fpar_FF_Nissa &p, const ipar_FF_Nissa &ip) {

		 return ip.FF;
		 };
  bf_FF.error=  [ ](const fpar_FF_Nissa &p, const ipar_FF_Nissa &ip) {

		 return ip.FF_err;
		 };

  bf_FF_ch2.ansatz= bf_FF.ansatz;
  bf_FF_ch2.measurement = bf_FF.measurement;
  bf_FF_ch2.error = bf_FF.error;
  
  //FIT FA for all xg
  Vfloat ch2_FA;
  distr_t_list F0_A_list(UseJack);
  distr_t_list D1_A_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {
  
  //fill the data
  vector<vector<ipar_FF_Nissa>> data(Njacks);
  vector<vector<ipar_FF_Nissa>> data_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
  for(auto &data_iboot: data) data_iboot.resize(Nens);
  for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data[ijack][iens].FF= FA_per_ens[iens].distr_list[ixg-1].distr[ijack];
      data[ijack][iens].FF_err= FA_per_ens[iens].err(ixg-1);
      if(data_2pts.Tag[iens] == "cA211a.12.48") data[ijack][iens].a = a_A.distr[ijack];
      else if(data_2pts.Tag[iens] == "cB211b.072.64") data[ijack][iens].a = a_B.distr[ijack];
      else if(data_2pts.Tag[iens] == "cC211a.06.80") data[ijack][iens].a = a_C.distr[ijack];
      else if(data_2pts.Tag[iens] == "cD211a.054.96") data[ijack][iens].a = a_D.distr[ijack];
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      //mean values
      if(ijack==0) {
	data_ch2[ijack][iens].FF= FA_per_ens[iens].ave(ixg-1);
	data_ch2[ijack][iens].FF_err= FA_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") data_ch2[ijack][iens].a = a_A.distr[ijack];
	else if(data_2pts.Tag[iens] == "cB211b.072.64") data_ch2[ijack][iens].a = a_B.distr[ijack];
	else if(data_2pts.Tag[iens] == "cC211a.06.80") data_ch2[ijack][iens].a = a_C.distr[ijack];
	else if(data_2pts.Tag[iens] == "cD211a.054.96") data_ch2[ijack][iens].a = a_D.distr[ijack];
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FA, xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();

  
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1);}
  //push_back retrieved parameters
  F0_A_list.distr_list.push_back(F0);
  D1_A_list.distr_list.push_back(D1);
  //push_back ch2
  ch2_FA.push_back( Bt_fit_ch2.get_ch2_ave());

 
  //print fit func
  distr_t_list FA_xg_to_print(UseJack);
  for(auto &a: a_to_print) FA_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2));
  Print_To_File({}, {a_to_print, FA_xg_to_print.ave(), FA_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FA FA_err");
 
  }


  
  //FIT FV for all xg
  Vfloat ch2_FV;
  bf_FF.Set_par_val("F0", -0.1, 0.001);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", -0.1, 0.001);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  distr_t_list F0_V_list(UseJack);
  distr_t_list D1_V_list(UseJack);
  for(int ixg=1;ixg<num_xg;ixg++) {

    //fill the data
    vector<vector<ipar_FF_Nissa>> data(Njacks);
    vector<vector<ipar_FF_Nissa>> data_ch2(1);
    //allocate space for output result
    boot_fit_data<fpar_FF_Nissa> Bt_fit;
    boot_fit_data<fpar_FF_Nissa> Bt_fit_ch2;
    for(auto &data_iboot: data) data_iboot.resize(Nens);
    for(auto &data_iboot: data_ch2) data_iboot.resize(Nens);
    for(int ijack=0;ijack<Njacks;ijack++) {
      for(int iens=0;iens<Nens;iens++) {
	data[ijack][iens].FF= FV_per_ens[iens].distr_list[ixg-1].distr[ijack];
	data[ijack][iens].FF_err= FV_per_ens[iens].err(ixg-1);
	if(data_2pts.Tag[iens] == "cA211a.12.48") data[ijack][iens].a = a_A.distr[ijack];
	else if(data_2pts.Tag[iens] == "cB211b.072.64") data[ijack][iens].a = a_B.distr[ijack];
	else if(data_2pts.Tag[iens] == "cC211a.06.80") data[ijack][iens].a = a_C.distr[ijack];
	else if(data_2pts.Tag[iens] == "cD211a.054.96") data[ijack][iens].a = a_D.distr[ijack];
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

	if(ijack==0) {
	  	data_ch2[ijack][iens].FF= FV_per_ens[iens].ave(ixg-1);
		data_ch2[ijack][iens].FF_err= FV_per_ens[iens].err(ixg-1);
		if(data_2pts.Tag[iens] == "cA211a.12.48") data_ch2[ijack][iens].a = a_A.distr[ijack];
		else if(data_2pts.Tag[iens] == "cB211b.072.64") data_ch2[ijack][iens].a = a_B.distr[ijack];
		else if(data_2pts.Tag[iens] == "cC211a.06.80") data_ch2[ijack][iens].a = a_C.distr[ijack];
		else if(data_2pts.Tag[iens] == "cD211a.054.96") data_ch2[ijack][iens].a = a_D.distr[ijack];
		else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	}
      }
    }
    
  //append
  bf_FF.Append_to_input_par(data);
  bf_FF_ch2.Append_to_input_par(data_ch2);
  //fit
  cout<<"Fitting FV, xg: "<<ixg<<endl;
  Bt_fit= bf_FF.Perform_bootstrap_fit();
  Bt_fit_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  distr_t F0(UseJack), D1(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { F0.distr.push_back( Bt_fit.par[ijack].F0); D1.distr.push_back( Bt_fit.par[ijack].D1);}
  //push_back retrieved parameters
  F0_V_list.distr_list.push_back(F0);
  D1_V_list.distr_list.push_back(D1);
  //push_back ch2
  ch2_FV.push_back( Bt_fit_ch2.get_ch2_ave());

  //print fit func
  distr_t_list FV_xg_to_print(UseJack);
  for(auto &a: a_to_print) FV_xg_to_print.distr_list.push_back( F0 + D1*pow(a*fmTGeV*Lambda_QCD,2));
  Print_To_File({}, {a_to_print, FV_xg_to_print.ave(), FV_xg_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_xg_"+to_string_with_precision(0.10*ixg,2)+".fit_func", "", "#a[fm] FV FV_err");
  }


          
  

 

  //Print continuum extrapolated form factors
  Print_To_File({}, {xg_t_list, F0_A_list.ave(), F0_A_list.err(), D1_A_list.ave(), D1_A_list.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_cont.dat", "", "#xg  F0   D1");
  Print_To_File({}, {xg_t_list, F0_V_list.ave(), F0_V_list.err(), D1_V_list.ave(), D1_V_list.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_cont.dat", "", "#xg  F0   D1");



  //VMD-like fit for FV and FA
  //###############################################################################################################

  class ipar_VMD {
  public:
    ipar_VMD() : FF(0.0), FF_err(0.0) {}
    double FF, FF_err, xg;
  };
  
  class fpar_VMD {
  public:
    fpar_VMD() {}
    fpar_VMD(const Vfloat &par) {
      if((signed)par.size() != 2) crash("In class fpar_VMD  class constructor Vfloat par has size != 2");
      A=par[0];
      M=par[1];
    }
    double A, M;
  };

  
  //init bootstrap fit
  bootstrap_fit<fpar_VMD,ipar_VMD> bf_VMD(Njacks);
  bf_VMD.set_warmup_lev(0); //sets warmup
  bf_VMD.Set_number_of_measurements(xg_t_list.size());
  bf_VMD.Set_verbosity(1);
  bf_VMD.Add_par("A", 0.07, 0.001);
  bf_VMD.Add_par("M", 1.3, 0.1);
  //fit on mean values to get ch2
  bootstrap_fit<fpar_VMD,ipar_VMD> bf_VMD_ch2(1);
  bf_VMD_ch2.set_warmup_lev(0); //sets warmup
  bf_VMD_ch2.Set_number_of_measurements(xg_t_list.size());
  bf_VMD_ch2.Set_verbosity(1);
  bf_VMD_ch2.Add_par("A", 0.07, 0.001);
  bf_VMD_ch2.Add_par("M", 1.3, 0.1);
  //bf_VMD.Fix_par("M", 1.25);
  //bf_VMD_ch2.Fix_par("M", 1.25);


  //ansatz
  bf_VMD.ansatz=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

		    double E_res= sqrt( p.M*p.M + ip.xg*ip.xg/4.0);
		    return p.A/( E_res*( E_res - (1.0 - ip.xg/2.0)) );
		 };
  bf_VMD.measurement=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

		 return ip.FF;
		 };
  bf_VMD.error=  [ ](const fpar_VMD &p, const ipar_VMD &ip) {

		 return ip.FF_err;
		 };

  bf_VMD_ch2.ansatz= bf_VMD.ansatz;
  bf_VMD_ch2.measurement = bf_VMD.measurement;
  bf_VMD_ch2.error = bf_VMD.error;


  //start fitting FA
  //fill the data
  vector<vector<ipar_VMD>> data_VMD(Njacks);
  vector<vector<ipar_VMD>> data_VMD_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_VMD> Bt_fit_FA_VMD;
  boot_fit_data<fpar_VMD> Bt_fit_FA_VMD_ch2;
  for(auto &data_iboot: data_VMD) data_iboot.resize(xg_t_list.size());
  for(auto &data_iboot: data_VMD_ch2) data_iboot.resize(xg_t_list.size());
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int ix=0;ix<(signed)xg_t_list.size();ix++) {
      data_VMD[ijack][ix].FF = F0_A_list.distr_list[ix].distr[ijack];
      data_VMD[ijack][ix].FF_err= F0_A_list.err(ix);
      data_VMD[ijack][ix].xg= xg_t_list[ix];
      if(ijack==0) {
	data_VMD_ch2[ijack][ix].FF = F0_A_list.ave(ix);
	data_VMD_ch2[ijack][ix].FF_err= F0_A_list.err(ix);
	data_VMD_ch2[ijack][ix].xg= xg_t_list[ix];

      }
    }
  }

  //append
  bf_VMD.Append_to_input_par(data_VMD);
  bf_VMD_ch2.Append_to_input_par(data_VMD_ch2);
  //fit
  cout<<"Fitting FA using VMD ansatz"<<endl;
  Bt_fit_FA_VMD= bf_VMD.Perform_bootstrap_fit();
  Bt_fit_FA_VMD_ch2= bf_VMD_ch2.Perform_bootstrap_fit();
  double ch2_red_FA_VMD= Bt_fit_FA_VMD_ch2.get_ch2_ave()/( xg_t_list.size() -2.0);

  //retrieve params
  distr_t Ampl_FA(UseJack), pole_FA(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { Ampl_FA.distr.push_back( Bt_fit_FA_VMD.par[ijack].A); pole_FA.distr.push_back( Bt_fit_FA_VMD.par[ijack].M);}


  //start fitting FV

  bf_VMD.Set_par_val("A", -0.07, 0.001);
  bf_VMD_ch2.Set_par_val("A", -0.07, 0.001);
  //bf_VMD.Fix_par("M", 1.073);
  //bf_VMD_ch2.Fix_par("M", 1.073);

  //allocate space for output result
  boot_fit_data<fpar_VMD> Bt_fit_FV_VMD;
  boot_fit_data<fpar_VMD> Bt_fit_FV_VMD_ch2;

  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int ix=0;ix<(signed)xg_t_list.size();ix++) {
      data_VMD[ijack][ix].FF = F0_V_list.distr_list[ix].distr[ijack];
      data_VMD[ijack][ix].FF_err= F0_V_list.err(ix);
      data_VMD[ijack][ix].xg= xg_t_list[ix];
      if(ijack==0) {
	data_VMD_ch2[ijack][ix].FF = F0_V_list.ave(ix);
	data_VMD_ch2[ijack][ix].FF_err= F0_V_list.err(ix);
	data_VMD_ch2[ijack][ix].xg= xg_t_list[ix];

      }
    }
  }

  //append
  bf_VMD.Append_to_input_par(data_VMD);
  bf_VMD_ch2.Append_to_input_par(data_VMD_ch2);
  //fit
  cout<<"Fitting FV using VMD ansatz"<<endl;
  Bt_fit_FV_VMD= bf_VMD.Perform_bootstrap_fit();
  Bt_fit_FV_VMD_ch2= bf_VMD_ch2.Perform_bootstrap_fit();

  double ch2_red_FV_VMD= Bt_fit_FV_VMD_ch2.get_ch2_ave()/( xg_t_list.size() -2.0);
 

  //retrieve params
  distr_t Ampl_FV(UseJack), pole_FV(UseJack);
  for(int ijack=0;ijack<Njacks;ijack++) { Ampl_FV.distr.push_back( Bt_fit_FV_VMD.par[ijack].A); pole_FV.distr.push_back( Bt_fit_FV_VMD.par[ijack].M);}

  
  distr_t_list FA_VMD_fit(UseJack), FV_VMD_fit(UseJack);
  //plot fit function
   for(auto &X: xg_to_spline_VMD) {
     ipar_VMD pp_VMD;
     pp_VMD.xg=X;
     distr_t FA_VMD_xg(UseJack), FV_VMD_xg(UseJack);
     for(int ijack=0;ijack<Njacks;ijack++) {
       FA_VMD_xg.distr.push_back( bf_VMD.ansatz( Bt_fit_FA_VMD.par[ijack], pp_VMD));
       FV_VMD_xg.distr.push_back( bf_VMD.ansatz( Bt_fit_FV_VMD.par[ijack], pp_VMD));
     }

     FA_VMD_fit.distr_list.push_back( FA_VMD_xg);
     FV_VMD_fit.distr_list.push_back( FV_VMD_xg);
  
   }

   //print

   string header_FA= "Ampl: "+to_string_with_precision(Ampl_FA.ave(),5)+" +- "+to_string_with_precision(Ampl_FA.err(),5)+" M^res/Mp: "+to_string_with_precision(pole_FA.ave(), 5)+" +- "+to_string_with_precision(pole_FA.err(), 5)+" ch2/dof: "+to_string_with_precision(ch2_red_FA_VMD  ,5);
   string header_FV= "Ampl: "+to_string_with_precision(Ampl_FV.ave(),5)+" +- "+to_string_with_precision(Ampl_FV.err(),5)+" M^res/Mp: "+to_string_with_precision(pole_FV.ave(), 5)+" +- "+to_string_with_precision(pole_FV.err(), 5)+" ch2/dof: "+to_string_with_precision(ch2_red_FV_VMD  ,5);
   Print_To_File({},{ xg_to_spline_VMD, FA_VMD_fit.ave(), FA_VMD_fit.err()} , "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_VMD.fit" , "", header_FA);
   Print_To_File({},{ xg_to_spline_VMD, FV_VMD_fit.ave(), FV_VMD_fit.err()} , "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_VMD.fit" , "", header_FV);







  //#################################################################################################################

  
  
  
  
  
  //interpolate continuum extrapolated form factors
  for(int ijack=0;ijack<Njacks;ijack++) {

    Vfloat FV_cont_ij, FA_cont_ij;
    for(int ixg=1;ixg<num_xg;ixg++) {
      FA_cont_ij.push_back( F0_A_list.distr_list[ixg-1].distr[ijack]);
      FV_cont_ij.push_back( F0_V_list.distr_list[ixg-1].distr[ijack]);
    }

    FA_cont_interpol_jacks.emplace_back( FA_cont_ij.begin(), FA_cont_ij.end(), 0.1, 0.1);
    FV_cont_interpol_jacks.emplace_back( FV_cont_ij.begin(), FV_cont_ij.end(), 0.1, 0.1);

  }


  auto FA_cont_interpol_distr= [&FA_cont_interpol_jacks](double xg) -> distr_t {


				 distr_t return_distr(UseJack);

				 for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FA_cont_interpol_jacks[ijack](xg));}

				 return return_distr;


			       };
  auto FV_cont_interpol_distr= [&FV_cont_interpol_jacks](double xg) -> distr_t {

				 
				 distr_t return_distr(UseJack);

				 for(int ijack=0; ijack<Njacks;ijack++) { return_distr.distr.push_back( FV_cont_interpol_jacks[ijack](xg));}

				 return return_distr;

			       };

 
  distr_t_list FA_interpol_cont_to_print_distr(UseJack);
  distr_t_list FV_interpol_cont_to_print_distr(UseJack);
  for(auto &X: xg_to_spline) {
    FA_interpol_cont_to_print_distr.distr_list.push_back( FA_cont_interpol_distr(X));
    FV_interpol_cont_to_print_distr.distr_list.push_back( FV_cont_interpol_distr(X));
  }

  Print_To_File({}, {xg_to_spline, FA_interpol_cont_to_print_distr.ave(), FA_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FA_interpol.dat", "", "#xg FA FA_err");
  Print_To_File({}, {xg_to_spline, FV_interpol_cont_to_print_distr.ave(), FV_interpol_cont_to_print_distr.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/FF/continuum/FV_interpol.dat", "", "#xg FV FV_err");


  //#################################################################################






  //##########################################################################################
  //##########################################################################################

  //perform continuum extrapolation of MP, fP, and MP/fP
  //MP
  bf_FF.Set_par_val("F0", MP_list.ave(0), MP_list.ave(0)/10);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", MP_list.ave(0), MP_list.ave(0)/10);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  double ch2_MP;
  distr_t F0_MP(UseJack);
  distr_t D1_MP(UseJack);
  vector<vector<ipar_FF_Nissa>> data_MP(Njacks);
  vector<vector<ipar_FF_Nissa>> data_MP_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit_MP;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_MP_ch2;
  for(auto &data_iboot: data_MP) data_iboot.resize(Nens);
  for(auto &data_iboot: data_MP_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data_MP[ijack][iens].FF= MP_list.distr_list[iens].distr[ijack];
      data_MP[ijack][iens].FF_err= MP_list.err(iens);
      if(data_2pts.Tag[iens] == "cA211a.12.48") data_MP[ijack][iens].a = a_A.distr[ijack];
      else if(data_2pts.Tag[iens] == "cB211b.072.64") data_MP[ijack][iens].a = a_B.distr[ijack];
      else if(data_2pts.Tag[iens] == "cC211a.06.80") data_MP[ijack][iens].a = a_C.distr[ijack];
      else if(data_2pts.Tag[iens] == "cD211a.054.96") data_MP[ijack][iens].a = a_D.distr[ijack];
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      if(ijack==0) {
	data_MP_ch2[ijack][iens].FF= MP_list.ave(iens);
	data_MP_ch2[ijack][iens].FF_err= MP_list.err(iens);
	if(data_2pts.Tag[iens] == "cA211a.12.48") data_MP_ch2[ijack][iens].a = a_A.distr[ijack];
	else if(data_2pts.Tag[iens] == "cB211b.072.64") data_MP_ch2[ijack][iens].a = a_B.distr[ijack];
	else if(data_2pts.Tag[iens] == "cC211a.06.80") data_MP_ch2[ijack][iens].a = a_C.distr[ijack];
	else if(data_2pts.Tag[iens] == "cD211a.054.96") data_MP_ch2[ijack][iens].a = a_D.distr[ijack];
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	
      }
    }
  }
  //append
  bf_FF.Append_to_input_par(data_MP);
  bf_FF_ch2.Append_to_input_par(data_MP_ch2);
  //fit
  cout<<"Fitting MP"<<endl;
  Bt_fit_MP= bf_FF.Perform_bootstrap_fit();
  Bt_fit_MP_ch2 = bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) { F0_MP.distr.push_back( Bt_fit_MP.par[ijack].F0); D1_MP.distr.push_back( Bt_fit_MP.par[ijack].D1);}
  //get ch2
  ch2_MP= Bt_fit_MP_ch2.get_ch2_ave();


  //FP
  double ch2_FP;
  bf_FF.Set_par_val("F0", FP_list.ave(0), FP_list.ave(0)/10);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", FP_list.ave(0), FP_list.ave(0)/10);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  distr_t F0_FP(UseJack);
  distr_t D1_FP(UseJack);
  vector<vector<ipar_FF_Nissa>> data_FP(Njacks);
  vector<vector<ipar_FF_Nissa>> data_FP_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit_FP;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_FP_ch2;
  for(auto &data_iboot: data_FP) data_iboot.resize(Nens);
  for(auto &data_iboot: data_FP_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data_FP[ijack][iens].FF= FP_list.distr_list[iens].distr[ijack];
      data_FP[ijack][iens].FF_err= FP_list.err(iens);
      if(data_2pts.Tag[iens] == "cA211a.12.48") data_FP[ijack][iens].a = a_A.distr[ijack];
      else if(data_2pts.Tag[iens] == "cB211b.072.64") data_FP[ijack][iens].a = a_B.distr[ijack];
      else if(data_2pts.Tag[iens] == "cC211a.06.80") data_FP[ijack][iens].a = a_C.distr[ijack];
      else if(data_2pts.Tag[iens] == "cD211a.054.96") data_FP[ijack][iens].a = a_D.distr[ijack];
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      if(ijack==0) {
	data_FP_ch2[ijack][iens].FF= FP_list.ave(iens);
	data_FP_ch2[ijack][iens].FF_err= FP_list.err(iens);
	if(data_2pts.Tag[iens] == "cA211a.12.48") data_FP_ch2[ijack][iens].a = a_A.distr[ijack];
	else if(data_2pts.Tag[iens] == "cB211b.072.64") data_FP_ch2[ijack][iens].a = a_B.distr[ijack];
	else if(data_2pts.Tag[iens] == "cC211a.06.80") data_FP_ch2[ijack][iens].a = a_C.distr[ijack];
	else if(data_2pts.Tag[iens] == "cD211a.054.96") data_FP_ch2[ijack][iens].a = a_D.distr[ijack];
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
	
      }
    }
  }
  //append
  bf_FF.Append_to_input_par(data_FP);
  bf_FF_ch2.Append_to_input_par(data_FP_ch2);
  //fit
  cout<<"Fitting FP"<<endl;
  Bt_fit_FP= bf_FF.Perform_bootstrap_fit();
  Bt_fit_FP_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) { F0_FP.distr.push_back( Bt_fit_FP.par[ijack].F0); D1_FP.distr.push_back( Bt_fit_FP.par[ijack].D1);}
  //get ch2
  ch2_FP = Bt_fit_FP_ch2.get_ch2_ave();



  //MP/FP
  double ch2_MP_ov_FP;
  bf_FF.Set_par_val("F0", MP_ov_FP_list.ave(0), MP_ov_FP_list.ave(0)/10);
  bf_FF.Set_par_val("D1", 1.0, 0.1);
  bf_FF_ch2.Set_par_val("F0", MP_ov_FP_list.ave(0), MP_ov_FP_list.ave(0)/10);
  bf_FF_ch2.Set_par_val("D1", 1.0, 0.1);
  distr_t F0_MP_ov_FP(UseJack);
  distr_t D1_MP_ov_FP(UseJack);
  vector<vector<ipar_FF_Nissa>> data_MP_ov_FP(Njacks);
  vector<vector<ipar_FF_Nissa>> data_MP_ov_FP_ch2(1);
  //allocate space for output result
  boot_fit_data<fpar_FF_Nissa> Bt_fit_MP_ov_FP;
  boot_fit_data<fpar_FF_Nissa> Bt_fit_MP_ov_FP_ch2;
  for(auto &data_iboot: data_MP_ov_FP) data_iboot.resize(Nens);
   for(auto &data_iboot: data_MP_ov_FP_ch2) data_iboot.resize(Nens);
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int iens=0;iens<Nens;iens++) {
      data_MP_ov_FP[ijack][iens].FF= MP_ov_FP_list.distr_list[iens].distr[ijack];
      data_MP_ov_FP[ijack][iens].FF_err= MP_ov_FP_list.err(iens);
      if(data_2pts.Tag[iens] == "cA211a.12.48") data_MP_ov_FP[ijack][iens].a = a_A.distr[ijack];
      else if(data_2pts.Tag[iens] == "cB211b.072.64") data_MP_ov_FP[ijack][iens].a = a_B.distr[ijack];
      else if(data_2pts.Tag[iens] == "cC211a.06.80") data_MP_ov_FP[ijack][iens].a = a_C.distr[ijack];
      else if(data_2pts.Tag[iens] == "cD211a.054.96") data_MP_ov_FP[ijack][iens].a = a_D.distr[ijack];
      else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");

      if(ijack==0) {
	data_MP_ov_FP_ch2[ijack][iens].FF= MP_ov_FP_list.ave(iens);
	data_MP_ov_FP_ch2[ijack][iens].FF_err= MP_ov_FP_list.err(iens);
	if(data_2pts.Tag[iens] == "cA211a.12.48") data_MP_ov_FP_ch2[ijack][iens].a = a_A.distr[ijack];
	else if(data_2pts.Tag[iens] == "cB211b.072.64") data_MP_ov_FP_ch2[ijack][iens].a = a_B.distr[ijack];
	else if(data_2pts.Tag[iens] == "cC211a.06.80") data_MP_ov_FP_ch2[ijack][iens].a = a_C.distr[ijack];
	else if(data_2pts.Tag[iens] == "cD211a.054.96") data_MP_ov_FP_ch2[ijack][iens].a = a_D.distr[ijack];
	else crash("Ens_tag: "+data_2pts.Tag[iens]+" not recognized");
      }
    }
  }
  //append
  bf_FF.Append_to_input_par(data_MP_ov_FP);
  bf_FF_ch2.Append_to_input_par(data_MP_ov_FP_ch2);
  //fit
  cout<<"Fitting MP/FP"<<endl;
  Bt_fit_MP_ov_FP= bf_FF.Perform_bootstrap_fit();
  Bt_fit_MP_ov_FP_ch2= bf_FF_ch2.Perform_bootstrap_fit();
  //retrieve parameters
  for(int ijack=0;ijack<Njacks;ijack++) { F0_MP_ov_FP.distr.push_back( Bt_fit_MP_ov_FP.par[ijack].F0); D1_MP_ov_FP.distr.push_back( Bt_fit_MP_ov_FP.par[ijack].D1);}
  //get ch2
  ch2_MP_ov_FP = Bt_fit_MP_ov_FP_ch2.get_ch2_ave();


  //print fitting functions for MP, FP, MP_ov_FP
  distr_t_list MP_to_print(UseJack);
  distr_t_list FP_to_print(UseJack);
  distr_t_list MP_ov_FP_to_print(UseJack);
  for(auto &a: a_to_print) {
    MP_to_print.distr_list.push_back( F0_MP + D1_MP*pow(a*fmTGeV*Lambda_QCD,2));
    FP_to_print.distr_list.push_back( F0_FP + D1_FP*pow(a*fmTGeV*Lambda_QCD,2));
    MP_ov_FP_to_print.distr_list.push_back( F0_MP_ov_FP + D1_MP_ov_FP*pow(a*fmTGeV*Lambda_QCD,2));
  }
  Print_To_File({}, {a_to_print, MP_to_print.ave(), MP_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/mass/masses.fit_func", "", "#a[fm]  MP MP_err");
  Print_To_File({}, {a_to_print, FP_to_print.ave(), FP_to_print.err(), MP_ov_FP_to_print.ave(), MP_ov_FP_to_print.err()}, "../data/ph_emission/"+ph_type+"/"+Meson+"/decay_const/fP.fit_func", "", "#a[fm]  FP FP_err   MP/FP  MP/FP_err");


  //print ch2 summary
  cout<<"#######   SUMMARY OF Ch2 ##########"<<endl;
  cout<<"### FA ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { cout<<"ch2(xg: "<<xg_t_list[ixg-1]<<"): "<<ch2_FA[ixg-1]<<endl; }
  cout<<"### FV ###"<<endl;
  for(int ixg=1;ixg<num_xg;ixg++) { cout<<"ch2(xg: "<<xg_t_list[ixg-1]<<"): "<<ch2_FV[ixg-1]<<endl; }
  cout<<"##########"<<endl;
  cout<<"ch2(MP): "<< ch2_MP<<endl;
  cout<<"ch2(FP): "<< ch2_FP<<endl;
  cout<<"ch2(MP/FP): "<<ch2_MP_ov_FP<<endl;
  cout<<"###################################"<<endl;


 
  
  //##########################################################################################
  //##########################################################################################
  
		
  


 

  }

  
 


  return;




}
