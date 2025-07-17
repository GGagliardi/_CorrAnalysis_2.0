#include "../include/K3lnu.h"
#include "RC_WI_analysis.h"
#include "numerics.h"
#include "scale_setting_main.h"
#include "stat.h"

using namespace std;


const double alpha = 1.0/137.035999;
const double MK_FLAG = 0.494600000;
Vfloat xg_list({0.1, 0.3, 0.5, 0.7, 0.9});
const bool UseJack = 1;
const int Njacks = 100;
const int Nboots=100;
const double qu = 2.0/3.0;
const double qd = -qu / 2.0;


An_cont_ret Analytic_continuation(const distr_t_list& C_in, const distr_t& E, const exc_state_info& exc_info, int tw, int TO,  string out) {

  An_cont_ret RET;
  if(TO != 1 && TO != 2) crash("Only time orderings implemented are 1, 2. Asked: "+to_string(TO));
  //cout<<"TO: "<<TO<<endl;
  int T= C_in.size();
  double r= (TO==1)?1.0:-1.0;
  int tmax=tw;
  int ts_opt=tw;
  distr_t_list C(UseJack); 
  //build correlation function
  for(int t=0;t<T;t++) {
    if(TO==1 && t<=tw)  C.distr_list.push_back( C_in.distr_list[tw-t] ); //first TO
    if(TO==2 && t<(T-tw)) { tmax=T-tw-1; ts_opt=T/2 ; (t==0)? (C.distr_list.push_back(0.0*Get_id_jack_distr(Njacks))):(C.distr_list.push_back( C_in.distr_list[tw+t])); } //second TO
  }
  //integrate
  distr_t_list F0(UseJack), F_meff(UseJack), F_m0(UseJack), F_mRES(UseJack);
  distr_t meff=1e3*Get_id_jack_distr(Njacks) ; //initial guess for effective mass (equivalent to meff=infty)
  distr_t x= EXP_D( -r*E);
  distr_t gap_1= EXP_D(-1.0*(exc_info.mGS + r*E));
  distr_t gap_2= EXP_D(-1.0*(exc_info.mRES + r*E));
  for(int t=0;t<=tmax;t++) {
    distr_t gap_3= EXP_D(-1.0*(meff + r*E));
    F0.distr_list.push_back(   (t==0?(0.0*Get_id_jack_distr(Njacks)):F0[t-1]) + C[t]*POW_D(x,t) );
    F_m0.distr_list.push_back(   F0[t] + C[t]*POW_D(x,t)*gap_1/(1-gap_1));
    F_mRES.distr_list.push_back( F0[t] + C[t]*POW_D(x,t)*gap_2/(1-gap_2));
    F_meff.distr_list.push_back( F0[t] + C[t]*POW_D(x,t)*gap_3/(1-gap_3));
    //update meff if good enough
    if(t != tmax)
      if( ((C[t]/C[t+1]).ave() > 1.0 ) && (LOG_D( C[t]/C[t+1]).err() < 0.1*LOG_D( C[t]/C[t+1]).ave() ) ) meff= LOG_D( C[t]/C[t+1]);

    //cout<<"meff: "<<meff.ave()<<" "<<meff.err()<<endl;
  }
  //populate return class
  RET.unsub_FF = F0[ts_opt];
  RET.mGS_FF = F_m0[ts_opt-4];
  RET.mRES_FF = F_mRES[ts_opt-4];
  RET.meff_FF = F_meff[ts_opt-4];
  RET.unsub_FF_list = F0;
  RET.mGS_FF_list= F_m0;
  RET.mRES_FF_list = F_mRES;
  RET.meff_FF_list = F_meff;

  //print
  if(out != "") RET.Print(out);
  
  return RET;
}

An_cont_ret operator+(const An_cont_ret &A, const An_cont_ret &B) {

  An_cont_ret ret;

  An_cont_ret A1=A;
  An_cont_ret A2=B;

  if(A.mGS_FF.size() != B.mGS_FF.size()) crash("operator A+B in An_cont_ret class called with A e B with different Njacks");

  int dim_A1 = A1.unsub_FF_list.size();
  int dim_A2 = A2.unsub_FF_list.size();
  
  //check if they have same dimension
  if(dim_A1 != dim_A2) {
    int r= dim_A2-dim_A1;
    for(int i= 0 ; i < abs(r); i++) {
      if(r > 0) { A1.unsub_FF_list.distr_list.push_back( A1.unsub_FF_list[dim_A1-1] );  A1.meff_FF_list.distr_list.push_back( A1.meff_FF_list[dim_A1-1] );  A1.mGS_FF_list.distr_list.push_back( A1.mGS_FF_list[dim_A1-1] );  A1.mRES_FF_list.distr_list.push_back( A1.mRES_FF_list[dim_A1-1] ); }
      else  { A2.unsub_FF_list.distr_list.push_back( A2.unsub_FF_list[dim_A2-1] );  A2.meff_FF_list.distr_list.push_back( A2.meff_FF_list[dim_A2-1] );  A2.mGS_FF_list.distr_list.push_back( A2.mGS_FF_list[dim_A2-1] );  A2.mRES_FF_list.distr_list.push_back( A2.mRES_FF_list[dim_A2-1] ); }

    }
  }

  ret.mGS_FF = A1.mGS_FF + A2.mGS_FF;
  ret.meff_FF = A1.meff_FF + A2.meff_FF;
  ret.mRES_FF = A1.mRES_FF + A2.mRES_FF;
  ret.unsub_FF = A1.unsub_FF + A2.unsub_FF;

  ret.mGS_FF_list = A1.mGS_FF_list + A2.mGS_FF_list;
  ret.meff_FF_list = A1.meff_FF_list + A2.meff_FF_list;
  ret.mRES_FF_list = A1.mRES_FF_list + A2.mRES_FF_list;
  ret.unsub_FF_list = A1.unsub_FF_list + A2.unsub_FF_list;

  return ret;
}

An_cont_ret operator-(const An_cont_ret &A, const An_cont_ret &B) {


  An_cont_ret ret;
  An_cont_ret A1=A;
  An_cont_ret A2=B;
    
  if(A.mGS_FF.size() != B.mGS_FF.size()) crash("operator A-B in An_cont_ret class called with A e B with different Njacks");
    
  int dim_A1 = A1.unsub_FF_list.size();
  int dim_A2 = A2.unsub_FF_list.size();
    
  //check if they have same dimension
  if(dim_A1 != dim_A2) {
    int r= dim_A2-dim_A1;
    for(int i= 0 ; i < abs(r); i++) {
      if(r > 0) { A1.unsub_FF_list.distr_list.push_back( A1.unsub_FF_list[dim_A1-1] );  A1.meff_FF_list.distr_list.push_back( A1.meff_FF_list[dim_A1-1] );  A1.mGS_FF_list.distr_list.push_back( A1.mGS_FF_list[dim_A1-1] );  A1.mRES_FF_list.distr_list.push_back( A1.mRES_FF_list[dim_A1-1] ); }
      else  { A2.unsub_FF_list.distr_list.push_back( A2.unsub_FF_list[dim_A2-1] );  A2.meff_FF_list.distr_list.push_back( A2.meff_FF_list[dim_A2-1] );  A2.mGS_FF_list.distr_list.push_back( A2.mGS_FF_list[dim_A2-1] );  A2.mRES_FF_list.distr_list.push_back( A2.mRES_FF_list[dim_A2-1] ); }

    }
  }
    
  ret.mGS_FF = A1.mGS_FF - A2.mGS_FF;
  ret.meff_FF = A1.meff_FF - A2.meff_FF;
  ret.mRES_FF = A1.mRES_FF - A2.mRES_FF;
  ret.unsub_FF = A1.unsub_FF - A2.unsub_FF;
  ret.mGS_FF_list = A1.mGS_FF_list - A2.mGS_FF_list;
  ret.meff_FF_list = A1.meff_FF_list - A2.meff_FF_list;
  ret.mRES_FF_list = A1.mRES_FF_list - A2.mRES_FF_list;
  ret.unsub_FF_list = A1.unsub_FF_list - A2.unsub_FF_list;
    
  return ret;
}

An_cont_ret operator+(const An_cont_ret &A, const distr_t &B) {


  An_cont_ret ret;
        
  if(A.mGS_FF.size() != B.size()) crash("operator A-B in An_cont_ret class called with A e B with different Njacks");
      
  ret.mGS_FF = A.mGS_FF +B;
  ret.meff_FF = A.meff_FF+B;
  ret.mRES_FF = A.mRES_FF +B;
  ret.unsub_FF = A.unsub_FF +B;
  ret.mGS_FF_list = A.mGS_FF_list +B;
  ret.meff_FF_list = A.meff_FF_list+B;
  ret.mRES_FF_list = A.mRES_FF_list+B;
  ret.unsub_FF_list = A.unsub_FF_list+B;
    
  return ret;
    
}

An_cont_ret operator-(const An_cont_ret &A, const distr_t &B) {


  An_cont_ret ret;
        
  if(A.mGS_FF.size() != B.size()) crash("operator A-B in An_cont_ret class called with A e B with different Njacks");
      
  ret.mGS_FF = A.mGS_FF - B;
  ret.meff_FF = A.meff_FF- B;
  ret.mRES_FF = A.mRES_FF - B;
  ret.unsub_FF = A.unsub_FF - B;
  ret.mGS_FF_list = A.mGS_FF_list - B;
  ret.meff_FF_list = A.meff_FF_list- B;
  ret.mRES_FF_list = A.mRES_FF_list- B;
  ret.unsub_FF_list = A.unsub_FF_list- B;
    
  return ret;
}


An_cont_ret operator*(const double &a, const An_cont_ret &A) {


  An_cont_ret ret;
        
   
      
  ret.mGS_FF = a*A.mGS_FF;
  ret.meff_FF = a*A.meff_FF;
  ret.mRES_FF = a*A.mRES_FF;
  ret.unsub_FF = a*A.unsub_FF;
  ret.mGS_FF_list = a*A.mGS_FF_list;
  ret.meff_FF_list = a*A.meff_FF_list;
  ret.mRES_FF_list = a*A.mRES_FF_list;
  ret.unsub_FF_list = a*A.unsub_FF_list;
    
  return ret;
    
}

An_cont_ret operator*(const An_cont_ret &A, const double &a) {

  return a*A;
}


An_cont_ret operator-(const distr_t &A, const An_cont_ret &B) {

  return  (-1.0*B) + A;
}

An_cont_ret operator+(const distr_t &A, const An_cont_ret &B) {

  return B+A;
}




void An_cont_ret::Print(string out) {
 
  Print_To_File({}, {this->unsub_FF_list.ave(), this->unsub_FF_list.err(), this->meff_FF_list.ave(), this->meff_FF_list.err(), this->mGS_FF_list.ave(), this->mGS_FF_list.err(), this->mRES_FF_list.ave(), this->mRES_FF_list.err()}, out, "", "");
  return;
}



void K3lnu() {

  Do_HLT_virtual();

  Analyze_Aprime();

  Get_electrounquenching();

  exit(-1);

  //get lattice spacing, quark masses and RCs
  scale_setting_info SCALE_SETTING_INFO= Get_scale_setting_info();
  RCs_info RCs_INFO= Get_RCs("ss");

  auto Sort_light_confs = [](string A, string B) {
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


  for(int ixg=0;ixg<(signed)xg_list.size();ixg++) {

    
    //read data

    data_t V12_l_tw1, A11_l_tw1, V12_s_tw1, A11_s_tw1;
    data_t V12_l_tw2, A11_l_tw2, V12_s_tw2, A11_s_tw2;
    data_t V12_0l_tw1, A11_0l_tw1, V12_0s_tw1, A11_0s_tw1;
    data_t V12_0l_tw2, A11_0l_tw2, V12_0s_tw2, A11_0s_tw2;
    data_t pt2_K, pt2_PI, pt2_K_loc;
    data_t A33_0l_tw1, A33_0s_tw1, A33_0l_tw2, A33_0s_tw2;


    string xg_string=to_string_with_precision(xg_list[ixg],1);

    V12_l_tw1.Read("../K3lnu_BB/tw1_l", "mes_contr_PHOTON_EMISSION_tw1_xg_"+xg_string+"_V12", "S0P5", Sort_light_confs);
    A11_l_tw1.Read("../K3lnu_BB/tw1_l", "mes_contr_PHOTON_EMISSION_tw1_xg_"+xg_string+"_A11", "S0P5", Sort_light_confs);
    V12_s_tw1.Read("../K3lnu_BB/tw1_s", "mes_contr_PHOTON_EMISSION_tw1_xg_"+xg_string+"_V12", "S0P5", Sort_light_confs);
    A11_s_tw1.Read("../K3lnu_BB/tw1_s", "mes_contr_PHOTON_EMISSION_tw1_xg_"+xg_string+"_A11", "S0P5", Sort_light_confs);

    V12_l_tw2.Read("../K3lnu_BB/tw2_l", "mes_contr_PHOTON_EMISSION_tw2_xg_"+xg_string+"_V12", "S0P5", Sort_light_confs);
    A11_l_tw2.Read("../K3lnu_BB/tw2_l", "mes_contr_PHOTON_EMISSION_tw2_xg_"+xg_string+"_A11", "S0P5", Sort_light_confs);
    V12_s_tw2.Read("../K3lnu_BB/tw2_s", "mes_contr_PHOTON_EMISSION_tw2_xg_"+xg_string+"_V12", "S0P5", Sort_light_confs);
    A11_s_tw2.Read("../K3lnu_BB/tw2_s", "mes_contr_PHOTON_EMISSION_tw2_xg_"+xg_string+"_A11", "S0P5", Sort_light_confs);

    V12_0l_tw1.Read("../K3lnu_BB/tw1_l", "mes_contr_PHOTON_EMISSION_tw1_xg_0_V12", "S0P5", Sort_light_confs);
    A11_0l_tw1.Read("../K3lnu_BB/tw1_l", "mes_contr_PHOTON_EMISSION_tw1_xg_0_A11", "S0P5", Sort_light_confs);
    A33_0l_tw1.Read("../K3lnu_BB/tw1_l", "mes_contr_PHOTON_EMISSION_tw1_xg_0_A33", "S0P5", Sort_light_confs);
    V12_0s_tw1.Read("../K3lnu_BB/tw1_s", "mes_contr_PHOTON_EMISSION_tw1_xg_0_V12", "S0P5", Sort_light_confs);
    A11_0s_tw1.Read("../K3lnu_BB/tw1_s", "mes_contr_PHOTON_EMISSION_tw1_xg_0_A11", "S0P5", Sort_light_confs);
    A33_0s_tw1.Read("../K3lnu_BB/tw1_s", "mes_contr_PHOTON_EMISSION_tw1_xg_0_A33", "S0P5", Sort_light_confs);
    
    V12_0l_tw2.Read("../K3lnu_BB/tw2_l", "mes_contr_PHOTON_EMISSION_tw2_xg_0_V12", "S0P5", Sort_light_confs);
    A11_0l_tw2.Read("../K3lnu_BB/tw2_l", "mes_contr_PHOTON_EMISSION_tw2_xg_0_A11", "S0P5", Sort_light_confs);
    A33_0l_tw2.Read("../K3lnu_BB/tw2_l", "mes_contr_PHOTON_EMISSION_tw2_xg_0_A33", "S0P5", Sort_light_confs);
    V12_0s_tw2.Read("../K3lnu_BB/tw2_s", "mes_contr_PHOTON_EMISSION_tw2_xg_0_V12", "S0P5", Sort_light_confs);
    A11_0s_tw2.Read("../K3lnu_BB/tw2_s", "mes_contr_PHOTON_EMISSION_tw2_xg_0_A11", "S0P5", Sort_light_confs);
    A33_0s_tw2.Read("../K3lnu_BB/tw2_s", "mes_contr_PHOTON_EMISSION_tw2_xg_0_A33", "S0P5", Sort_light_confs);
    
    pt2_K.Read("../K3lnu_BB/tw1_l", "mes_contr_MES_2PT_K", "P5P5", Sort_light_confs);
    pt2_K_loc.Read("../K3lnu_BB/tw1_s", "mes_contr_MES_2PT_K_loc", "P5P5", Sort_light_confs);
    pt2_PI.Read("../K3lnu_BB/tw1_l", "mes_contr_MES_2PT_pi", "P5P5", Sort_light_confs);


    boost::filesystem::create_directory("../data/K3lnu");

    int Nens=pt2_K.Tag.size();

    for(int iens=0;iens<Nens;iens++) {

      string Ens=pt2_K.Tag[iens];
      cout<<"Analyzing ensemble: "<<Ens<<" xg: "<<xg_string<<endl;

      boost::filesystem::create_directory("../data/K3lnu/"+Ens);
      boost::filesystem::create_directory("../data/K3lnu/"+Ens+"/tw1");
      boost::filesystem::create_directory("../data/K3lnu/"+Ens+"/tw2");
      

      distr_t a_distr= SCALE_SETTING_INFO.a_B;
      int RC_ens=-1;
      for(int b=0;b<(signed)RCs_INFO.Ens.size();b++) { if(RCs_INFO.Ens[b]== "cB211b.072.64") RC_ens=b;}
      if(RC_ens==-1) crash("RC ens not found");
      distr_t Zv= RCs_INFO.Zv[RC_ens];
      distr_t Za= RCs_INFO.Za[RC_ens];
      int tw1=28; int tw2=-1;
      int tw2_ens=-1;
      if(Ens=="cB211b.072.48" || Ens=="cB211b.072.64") {
	for(int b=0;b<(signed)V12_l_tw2.Tag.size();b++) { if (V12_l_tw2.Tag[b] == Ens) tw2_ens=b;}
	if(tw2_ens==-1) crash("Ens not found");
	tw1=25; tw2=32;
      }

      CorrAnalysis Corr(UseJack,Njacks,100);
      Corr.Nt = pt2_K.nrows[iens];
      int T= Corr.Nt;
      int L= Corr.Nt/2;
      double ml=0.00072;
      double ms=0.01825;
      double xg=xg_list[ixg];
      double k1= 2*M_PI/L;

      cout<<"Analyzing correlators!"<<endl;


      bool Analyze_tw2= (Ens != "cB211b.072.96");
      //3pt
      distr_t_list V12_l_tw1_distr(UseJack), A11_l_tw1_distr(UseJack), V12_s_tw1_distr(UseJack), A11_s_tw1_distr(UseJack);
      distr_t_list V12_l_tw2_distr(UseJack), A11_l_tw2_distr(UseJack), V12_s_tw2_distr(UseJack), A11_s_tw2_distr(UseJack);
      //3pt 0-mom
      distr_t_list V12_0l_tw1_distr(UseJack), A11_0l_tw1_distr(UseJack), V12_0s_tw1_distr(UseJack), A11_0s_tw1_distr(UseJack), A33_0l_tw1_distr(UseJack), A33_0s_tw1_distr(UseJack);
      distr_t_list V12_0l_tw2_distr(UseJack), A11_0l_tw2_distr(UseJack), V12_0s_tw2_distr(UseJack), A11_0s_tw2_distr(UseJack), A33_0l_tw2_distr(UseJack), A33_0s_tw2_distr(UseJack);
      //2pt
      distr_t_list PT2_K_distr(UseJack), PT2_PI_distr(UseJack), PT2_K_loc_distr(UseJack);
      Corr.Reflection_sign=1;
      PT2_K_distr= Corr.corr_t(pt2_K.col(0)[iens],"");
      PT2_K_loc_distr = Corr.corr_t(pt2_K_loc.col(0)[iens], "");
      PT2_PI_distr= Corr.corr_t(pt2_PI.col(0)[iens],"");
      //Analyze 2PT
      distr_t_list MK_eff= Corr.effective_mass_t(PT2_K_distr, "../data/K3lnu/"+Ens+"/eff_K");
      distr_t_list MK_eff_loc= Corr.effective_mass_t(PT2_K_loc_distr, "../data/K3lnu/"+Ens+"/eff_K_loc");
      distr_t_list MPI_eff=Corr.effective_mass_t(PT2_PI_distr, "../data/K3lnu/"+Ens+"/eff_PI");
      Corr.Tmin=22;
      Corr.Tmax=35;
      if(Ens != "cB211b.072.48") Corr.Tmax=50;
      if(Ens=="cB211b.072.96") {Corr.Tmin=25; Corr.Tmax=60;}
      distr_t_list MEM_K= Corr.matrix_element_t(PT2_K_distr, "../data/K3lnu/"+Ens+"/MEM_K");
      distr_t MK= Corr.Fit_distr(MK_eff);
      distr_t MK_loc= Corr.Fit_distr(MK_eff_loc);
      distr_t_list FK_eff= (ms+ml)*Corr.residue_t(PT2_K_loc_distr, "")/(MK_loc*SINH_D(MK_loc)*Corr.matrix_element_t(PT2_K_distr, ""));
      distr_t FK= Corr.Fit_distr(FK_eff);
      distr_t Mpi= Corr.Fit_distr(MPI_eff);
      distr_t Eg= 0.5*xg*MK;
      distr_t Eg0 = 0.0*Get_id_jack_distr(Njacks);


      cout<<"MK: "<<(MK/a_distr).ave()<<" "<<(MK/a_distr).err()<<endl;
      cout<<"FK: "<<(FK/a_distr).ave()<<" "<<(FK/a_distr).err()<<endl;
      cout<<"Mpi: "<<(Mpi/a_distr).ave()<<" "<<(Mpi/a_distr).err()<<endl;
      cout<<"Zv: "<<Zv.ave()<<" "<<Zv.err()<<" Za: "<<Za.ave()<<" "<<Za.err()<<endl;
      cout<<"Eg: "<<Eg.ave()<<" "<<Eg.err()<<endl;
      
      distr_t FK_v_Eg=FK/Eg;
    
      
      distr_t_list Z_tw1= 1.0/( MEM_K*EXP_D(-1.0*MK*tw1)/(2.0*MK));
      distr_t_list Z_tw2(UseJack);
      if(Analyze_tw2) Z_tw2= 1.0/( MEM_K*EXP_D(-1.0*MK*tw2)/(2.0*MK));

      D(1);
      
      
            
      //analyze PT3
      //non-zero mom
      Corr.Perform_Nt_t_average=0;
      V12_l_tw1_distr= Za*Z_tw1*Corr.corr_t( V12_l_tw1.col(1)[iens],  "../data/K3lnu/"+Ens+"/tw1/V12_3PT_l_xg_"+xg_string)/Eg;
      A11_l_tw1_distr= Zv*Z_tw1*Corr.corr_t( A11_l_tw1.col(0)[iens], "")/Eg;
      V12_s_tw1_distr= Za*Z_tw1*Corr.corr_t( V12_s_tw1.col(1)[iens], "")/Eg;
      A11_s_tw1_distr= Zv*Z_tw1*Corr.corr_t( A11_s_tw1.col(0)[iens], "")/Eg;
      //zero mom
      V12_0l_tw1_distr= Za*Z_tw1*Corr.corr_t( V12_0l_tw1.col(1)[iens], "")/Eg;
      A11_0l_tw1_distr= Zv*Z_tw1*Corr.corr_t( A11_0l_tw1.col(0)[iens], "")/Eg;
      A33_0l_tw1_distr= Zv*Z_tw1*Corr.corr_t( A33_0l_tw1.col(0)[iens], "")/Eg;
      V12_0s_tw1_distr= Za*Z_tw1*Corr.corr_t( V12_0s_tw1.col(1)[iens], "")/Eg;
      A11_0s_tw1_distr= Zv*Z_tw1*Corr.corr_t( A11_0s_tw1.col(0)[iens], "")/Eg;
      A33_0s_tw1_distr= Zv*Z_tw1*Corr.corr_t( A33_0s_tw1.col(0)[iens], "")/Eg;
      
      if(Analyze_tw2) {
	//non-zero mom
	V12_l_tw2_distr= Za*Z_tw2*Corr.corr_t( V12_l_tw2.col(1)[tw2_ens], "")/Eg;
	A11_l_tw2_distr= Zv*Z_tw2*Corr.corr_t( A11_l_tw2.col(0)[tw2_ens], "")/Eg;
	V12_s_tw2_distr= Za*Z_tw2*Corr.corr_t( V12_s_tw2.col(1)[tw2_ens], "")/Eg;
	A11_s_tw2_distr= Zv*Z_tw2*Corr.corr_t( A11_s_tw2.col(0)[tw2_ens], "")/Eg;
	//zero mom
	V12_0l_tw2_distr= Za*Z_tw2*Corr.corr_t( V12_0l_tw2.col(1)[tw2_ens], "")/Eg;
	A11_0l_tw2_distr= Zv*Z_tw2*Corr.corr_t( A11_0l_tw2.col(0)[tw2_ens], "")/Eg;
	A33_0l_tw2_distr= Zv*Z_tw2*Corr.corr_t( A33_0l_tw2.col(0)[tw2_ens], "")/Eg;
	V12_0s_tw2_distr= Za*Z_tw2*Corr.corr_t( V12_0s_tw2.col(1)[tw2_ens], "")/Eg;
	A11_0s_tw2_distr= Zv*Z_tw2*Corr.corr_t( A11_0s_tw2.col(0)[tw2_ens], "")/Eg;
	A33_0s_tw2_distr= Zv*Z_tw2*Corr.corr_t( A33_0s_tw2.col(0)[tw2_ens], "")/Eg;
      }

      Print_To_File({}, { ((V12_l_tw1_distr- V12_0l_tw1_distr)*Eg).ave(), ((V12_l_tw1_distr-V12_0l_tw1_distr)*Eg).err() }, "../data/K3lnu/"+Ens+"/tw1/V12_l_xg_"+xg_string, "", "");
      if(Analyze_tw2 )    Print_To_File({}, { ((V12_l_tw2_distr- V12_0l_tw2_distr)*Eg).ave(), ((V12_l_tw2_distr-V12_0l_tw2_distr)*Eg).err() }, "../data/K3lnu/"+Ens+"/tw2/V12_l_xg_"+xg_string, "", "");
      
      
      //get FK from 3pt
      distr_t FK_3pt_l(UseJack);
      distr_t FK_3pt_s(UseJack);
      distr_t FK_3pt_l_tw2(UseJack);
      distr_t FK_3pt_s_tw2(UseJack);

      distr_t FK_3pt_MOT_l(UseJack);
      distr_t FK_3pt_MOT_s(UseJack);
      distr_t FK_3pt_MOT_l_tw2(UseJack);
      distr_t FK_3pt_MOT_s_tw2(UseJack);
      
      for(int t=0;t<T;t++) {
	FK_3pt_l = (t==0)?(0.0*Get_id_jack_distr(Njacks)):(FK_3pt_l + A11_0l_tw1_distr[t]);
	FK_3pt_s = (t==0)?(0.0*Get_id_jack_distr(Njacks)):(FK_3pt_s + A11_0s_tw1_distr[t]);
	FK_3pt_MOT_l = (t==0)?(0.0*Get_id_jack_distr(Njacks)):(FK_3pt_MOT_l + A11_l_tw1_distr[t]);
	FK_3pt_MOT_s = (t==0)?(0.0*Get_id_jack_distr(Njacks)):(FK_3pt_MOT_s + A11_s_tw1_distr[t]);
	if(Analyze_tw2) {
	  FK_3pt_l_tw2 = (t==0)?(0.0*Get_id_jack_distr(Njacks)):(FK_3pt_l_tw2 + A11_0l_tw2_distr[t]);
	  FK_3pt_s_tw2 = (t==0)?(0.0*Get_id_jack_distr(Njacks)):(FK_3pt_s_tw2 + A11_0s_tw2_distr[t]);
	  FK_3pt_MOT_l_tw2 = (t==0)?(0.0*Get_id_jack_distr(Njacks)):(FK_3pt_MOT_l_tw2 + A11_l_tw2_distr[t]);
	  FK_3pt_MOT_s_tw2 = (t==0)?(0.0*Get_id_jack_distr(Njacks)):(FK_3pt_MOT_s_tw2 + A11_s_tw2_distr[t]);
	}
      }


      distr_t Zr_l_tw1 = FK_v_Eg/FK_3pt_l;
      distr_t Zr_s_tw1 = FK_v_Eg/FK_3pt_s;
      distr_t Zr_l_tw2(UseJack), Zr_s_tw2(UseJack);
      if(Analyze_tw2) { Zr_l_tw2= FK_v_Eg/FK_3pt_l_tw2; Zr_s_tw2= FK_v_Eg/FK_3pt_s_tw2; }
      cout<<"rew_l: "<<Zr_l_tw1.ave()<<" "<<Zr_l_tw1.err()<<endl;
      cout<<"rew_s: "<<Zr_s_tw1.ave()<<" "<<Zr_s_tw1.err()<<endl;
      if(Analyze_tw2) {
	cout<<"rew_l(tw2): "<<Zr_l_tw2.ave()<<" "<<Zr_l_tw2.err()<<endl;
	cout<<"rew_s(tw2): "<<Zr_s_tw2.ave()<<" "<<Zr_s_tw2.err()<<endl;
      }
      
      
      //noise subtraction in vector component
      distr_t_list CFV_tw1_l = V12_l_tw1_distr -V12_0l_tw1_distr;
      distr_t_list CFV_tw1_s = V12_s_tw1_distr -V12_0s_tw1_distr;
      distr_t_list CFV_tw2_l(UseJack), CFV_tw2_s(UseJack);
      if(Analyze_tw2) {
	CFV_tw2_l = V12_l_tw2_distr -V12_0l_tw2_distr;
	CFV_tw2_s = V12_s_tw2_distr -V12_0s_tw2_distr;
      }

      
      //######### give info on excited state ############//
      distr_t PI_PI = 2.0*SQRT_D( Mpi*Mpi + k1*k1 );
      distr_t PI_PI_MOT = SQRT_D( POW_D(PI_PI,2) + Eg*Eg);
      distr_t K_PI  = SQRT_D(MK*MK + k1*k1) + SQRT_D(Mpi*Mpi + k1*k1);
      distr_t K_PI_MOT = SQRT_D( POW_D(K_PI,2) + Eg*Eg);
      distr_t K_PI_PI = SQRT_D(MK*MK + 2*k1*k1) + 2*SQRT_D(Mpi*Mpi + k1*k1);
      distr_t K_PI_PI_MOT = SQRT_D( POW_D(K_PI_PI,2) + Eg*Eg);
      distr_t Kstar = (0.89167*a_distr);
      distr_t Kstar_MOT= SQRT_D( POW_D(Kstar,2) + Eg*Eg);
      distr_t K1 = (1.253*a_distr);
      distr_t K1_MOT= SQRT_D( POW_D(K1,2) + Eg*Eg);
      distr_t rho= (0.775*a_distr);
      distr_t rho_MOT= SQRT_D( POW_D(rho,2) + Eg*Eg);
      distr_t phi= (1.019*a_distr);
      distr_t phi_MOT= SQRT_D( POW_D(phi,2) + Eg*Eg);
      exc_state_info TO1_V, TO1_A, TO1_A_0mom,  TO2_l, TO2_s, TO2_l_0mom, TO2_s_0mom  , null;
      TO1_V.mGS = K_PI_MOT-MK; TO1_V.mRES = Kstar_MOT-MK;
      TO1_A_0mom.mGS = K_PI_PI-MK; TO1_A_0mom.mRES = K1-MK;
      TO1_A.mGS = K_PI_PI_MOT-MK; TO1_A.mRES = K1_MOT-MK;
      TO2_l_0mom.mGS = PI_PI; TO2_l_0mom.mRES= rho;
      TO2_l.mGS = PI_PI_MOT; TO2_l.mRES = rho_MOT;
      TO2_s_0mom.mGS = phi; TO2_s_0mom.mRES= phi;
      TO2_s.mGS= phi_MOT; TO2_s.mRES=phi_MOT;
      null.mGS = 1e4*Get_id_jack_distr(Njacks); null.mRES=null.mGS;
      //################################################//

    
    
      cout<<"######### Performing analytic continuation tw1 ############"<<endl;
      //FV
      An_cont_ret FV_I_l_tw1= Analytic_continuation( CFV_tw1_l, Eg, TO1_V, tw1, 1, "../data/K3lnu/"+Ens+"/tw1/FV_I_l_xg_"+xg_string);
      An_cont_ret FV_I_s_tw1= Analytic_continuation( CFV_tw1_s, Eg, TO1_V, tw1, 1, "../data/K3lnu/"+Ens+"/tw1/FV_I_s_xg_"+xg_string);
      An_cont_ret FV_I_tot_tw1= Analytic_continuation( qu*CFV_tw1_l + qd*CFV_tw1_s, Eg, TO1_V, tw1, 1, "../data/K3lnu/"+Ens+"/tw1/FV_I_TOT_xg_"+xg_string);
      An_cont_ret FV_II_l_tw1= Analytic_continuation( CFV_tw1_l, Eg, TO2_l, tw1, 2, "../data/K3lnu/"+Ens+"/tw1/FV_II_l_xg_"+xg_string);
      An_cont_ret FV_II_s_tw1= Analytic_continuation( CFV_tw1_s, Eg, TO2_s, tw1, 2, "../data/K3lnu/"+Ens+"/tw1/FV_II_s_xg_"+xg_string);
      //FA
      An_cont_ret FA_I_l_tw1= Analytic_continuation( -1.0*A11_l_tw1_distr, Eg, TO1_A, tw1, 1, "");
      An_cont_ret FA_I_s_tw1= Analytic_continuation( -1.0*A11_s_tw1_distr, Eg, TO1_A, tw1, 1, "");
      An_cont_ret FA_I_tot_tw1= Analytic_continuation( -qu*A11_l_tw1_distr - (-1)*qd*A11_s_tw1_distr, Eg, TO1_A, tw1, 1, "");
      An_cont_ret FA_II_l_tw1= Analytic_continuation( -1.0*A11_l_tw1_distr, Eg, TO2_l, tw1, 2, "");
      An_cont_ret FA_II_s_tw1= Analytic_continuation( -1.0*A11_s_tw1_distr, Eg, TO2_s, tw1, 2, "");
      //FA zero momentum
      An_cont_ret FA_I_0l_tw1= Analytic_continuation( -1.0*A11_l_tw1_distr, Eg0, TO1_A, tw1, 1, "");
      An_cont_ret FA_I_0s_tw1= Analytic_continuation( -1.0*A11_s_tw1_distr, Eg0, TO1_A, tw1, 1, "");
      An_cont_ret FA_0I_tot_tw1= Analytic_continuation( -qu*A11_l_tw1_distr - (-1)*qd*A11_s_tw1_distr, Eg0, TO1_A, tw1, 1, "");
      An_cont_ret FA_II_0l_tw1= Analytic_continuation( -1.0*A11_l_tw1_distr, Eg0, TO2_l, tw1, 2, "");
      An_cont_ret FA_II_0s_tw1= Analytic_continuation( -1.0*A11_s_tw1_distr, Eg0, TO2_s, tw1, 2, "");
      //print
      (FA_I_l_tw1 -FA_I_0l_tw1).Print( "../data/K3lnu/"+Ens+"/tw1/FA_I_l_xg_"+xg_string);
      (FA_I_s_tw1-FA_I_0s_tw1).Print( "../data/K3lnu/"+Ens+"/tw1/FA_I_s_xg_"+xg_string);
      (FA_II_l_tw1-FA_II_0l_tw1).Print( "../data/K3lnu/"+Ens+"/tw1/FA_II_l_xg_"+xg_string);
      (FA_II_s_tw1-FA_II_0s_tw1).Print( "../data/K3lnu/"+Ens+"/tw1/FA_II_s_xg_"+xg_string);
      //sum first and second TO contributions
      (FA_I_l_tw1+ FA_II_l_tw1-FA_I_0l_tw1-FA_II_0l_tw1  ).Print("../data/K3lnu/"+Ens+"/tw1/FA_l_xg_"+xg_string);
      (FA_I_s_tw1+ FA_II_s_tw1-FA_I_0s_tw1-FA_II_0s_tw1  ).Print("../data/K3lnu/"+Ens+"/tw1/FA_s_xg_"+xg_string);
      //zero energy
      An_cont_ret zero_E_I_l= Analytic_continuation( 0.5*(A11_0l_tw1_distr + A11_0l_tw1_distr) -A11_l_tw1_distr, Eg0, null, tw1, 1, "../data/K3lnu/"+Ens+"/tw1/FA_null_I_l_xg_"+xg_string); 
      An_cont_ret zero_E_II_l= Analytic_continuation( 0.5*(A11_0l_tw1_distr + A11_0l_tw1_distr)-A11_l_tw1_distr, Eg0, null, tw1, 2, "../data/K3lnu/"+Ens+"/tw1/FA_null_II_l_xg_"+xg_string);
      An_cont_ret zero_E_I_s= Analytic_continuation( 0.5*(A11_0s_tw1_distr+A11_0s_tw1_distr)-A11_s_tw1_distr, Eg0, null, tw1, 1, "../data/K3lnu/"+Ens+"/tw1/FA_null_I_s_xg_"+xg_string); 
      An_cont_ret zero_E_II_s= Analytic_continuation( 0.5*(A11_0s_tw1_distr+A11_0s_tw1_distr)-A11_s_tw1_distr, Eg0, null, tw1, 2, "../data/K3lnu/"+Ens+"/tw1/FA_null_II_s_xg_"+xg_string); 
      (zero_E_I_l+zero_E_II_l).Print("../data/K3lnu/"+Ens+"/tw1/FA_null_l_xg_"+xg_string);
      (zero_E_I_s+zero_E_II_s).Print("../data/K3lnu/"+Ens+"/tw1/FA_null_s_xg_"+xg_string);

    
       
      if(Analyze_tw2) {
	cout<<"######### Performing analytic continuation tw2 ############"<<endl;
	//FV
	An_cont_ret FV_I_l_tw2= Analytic_continuation( CFV_tw2_l, Eg, TO1_V, tw2, 1, "../data/K3lnu/"+Ens+"/tw2/FV_I_l_xg_"+xg_string);
	An_cont_ret FV_I_s_tw2= Analytic_continuation( CFV_tw2_s, Eg, TO1_V, tw2, 1, "../data/K3lnu/"+Ens+"/tw2/FV_I_s_xg_"+xg_string);
	An_cont_ret FV_I_tot_tw2= Analytic_continuation( qu*CFV_tw2_l + qd*CFV_tw2_s, Eg, TO1_V, tw2, 1, "../data/K3lnu/"+Ens+"/tw2/FV_I_TOT_xg_"+xg_string);
	An_cont_ret FV_II_l_tw2= Analytic_continuation( CFV_tw2_l, Eg, TO2_l, tw2, 2, "../data/K3lnu/"+Ens+"/tw2/FV_II_l_xg_"+xg_string);
	An_cont_ret FV_II_s_tw2= Analytic_continuation( CFV_tw2_s, Eg, TO2_s, tw2, 2, "../data/K3lnu/"+Ens+"/tw2/FV_II_s_xg_"+xg_string);
	//FA
	An_cont_ret FA_I_l_tw2= Analytic_continuation( -1.0*A11_l_tw2_distr, Eg, TO1_A, tw2, 1, "");
	An_cont_ret FA_I_s_tw2= Analytic_continuation( -1.0*A11_s_tw2_distr, Eg, TO1_A, tw2, 1, "");
	An_cont_ret FA_I_tot_tw2= Analytic_continuation( -qu*A11_l_tw2_distr - (-1)*qd*A11_s_tw2_distr, Eg, TO1_A, tw2, 1, "");
	An_cont_ret FA_II_l_tw2= Analytic_continuation( -1.0*A11_l_tw2_distr, Eg, TO2_l, tw2, 2, "");
	An_cont_ret FA_II_s_tw2= Analytic_continuation( -1.0*A11_s_tw2_distr, Eg, TO2_s, tw2, 2, "");
	//FA zero momentum
	An_cont_ret FA_I_0l_tw2= Analytic_continuation( -1.0*A11_l_tw2_distr, Eg0, TO1_A, tw2, 1, "");
	An_cont_ret FA_I_0s_tw2= Analytic_continuation( -1.0*A11_s_tw2_distr, Eg0, TO1_A, tw2, 1, "");
	An_cont_ret FA_0I_tot_tw2= Analytic_continuation( -qu*A11_l_tw2_distr - (-1)*qd*A11_s_tw2_distr, Eg0, TO1_A, tw2, 1, "");
	An_cont_ret FA_II_0l_tw2= Analytic_continuation( -1.0*A11_l_tw2_distr, Eg0, TO2_l, tw2, 2, "");
	An_cont_ret FA_II_0s_tw2= Analytic_continuation( -1.0*A11_s_tw2_distr, Eg0, TO2_s, tw2, 2, "");
	//print
	(FA_I_l_tw2 -FA_I_0l_tw2).Print( "../data/K3lnu/"+Ens+"/tw2/FA_I_l_xg_"+xg_string);
	(FA_I_s_tw2-FA_I_0s_tw2).Print( "../data/K3lnu/"+Ens+"/tw2/FA_I_s_xg_"+xg_string);
	(FA_II_l_tw2-FA_II_0l_tw2).Print( "../data/K3lnu/"+Ens+"/tw2/FA_II_l_xg_"+xg_string);
	(FA_II_s_tw2-FA_II_0s_tw2).Print( "../data/K3lnu/"+Ens+"/tw2/FA_II_s_xg_"+xg_string);
	//sum first and second TO contributions
	(FA_I_l_tw2+ FA_II_l_tw2-FA_I_0l_tw2-FA_II_0l_tw2 ).Print("../data/K3lnu/"+Ens+"/tw2/FA_l_xg_"+xg_string);
	(FA_I_s_tw2+ FA_II_s_tw2-FA_I_0s_tw2-FA_II_0s_tw2 ).Print("../data/K3lnu/"+Ens+"/tw2/FA_s_xg_"+xg_string);
	//zero energy
	An_cont_ret zero_E_I_l= Analytic_continuation( 0.5*(A11_0l_tw2_distr + A11_0l_tw2_distr) -A11_l_tw2_distr, Eg0, null, tw2, 1, "../data/K3lnu/"+Ens+"/tw2/FA_null_I_l_xg_"+xg_string); 
	An_cont_ret zero_E_II_l= Analytic_continuation( 0.5*(A11_0l_tw2_distr + A11_0l_tw2_distr)-A11_l_tw2_distr, Eg0, null, tw2, 2, "../data/K3lnu/"+Ens+"/tw2/FA_null_II_l_xg_"+xg_string);
	An_cont_ret zero_E_I_s= Analytic_continuation( A11_0s_tw2_distr-A11_s_tw2_distr, Eg0, null, tw2, 1, "../data/K3lnu/"+Ens+"/tw2/FA_null_I_s_xg_"+xg_string); 
	An_cont_ret zero_E_II_s= Analytic_continuation( A11_0s_tw2_distr-A11_s_tw2_distr, Eg0, null, tw2, 2, "../data/K3lnu/"+Ens+"/tw2/FA_null_II_s_xg_"+xg_string); 
	(zero_E_I_l+zero_E_II_l).Print("../data/K3lnu/"+Ens+"/tw2/FA_null_l_xg_"+xg_string);
	(zero_E_I_s+zero_E_II_s).Print("../data/K3lnu/"+Ens+"/tw2/FA_null_s_xg_"+xg_string);
	
	

      }
      cout<<"DONE!"<<endl;
    }
  }


  cout<<"Exiting!"<<endl;

  return;

}


void Get_electrounquenching() {

  int Nhits=10;

  //only FV for the moment

  auto Sort_light = [](string A, string B) {
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


  vector<data_t> PT2_P5(Nhits), PT2_V1(Nhits), PT2_V2(Nhits);
  vector<data_t> PT2_P5_TH(Nhits), PT2_V1_TH(Nhits), PT2_V2_TH(Nhits);
  vector<data_t> PT2_P5_MTH(Nhits), PT2_V1_MTH(Nhits), PT2_V2_MTH(Nhits);
  vector<data_t> BUB_V1_R0(Nhits), BUB_V1_R1(Nhits), BUB_V2_R0(Nhits), BUB_V2_R1(Nhits);
  vector<data_t> BUB_V1_R0_TH(Nhits), BUB_V1_R1_TH(Nhits), BUB_V2_R0_TH(Nhits), BUB_V2_R1_TH(Nhits);


  for(int ihit=0;ihit<Nhits;ihit++) {

    PT2_P5[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_2PTH_0_hit"+to_string(ihit), "P5P5", Sort_light);
    PT2_P5_TH[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_2PT_TH_hit"+to_string(ihit), "P5P5", Sort_light);
    PT2_P5_MTH[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_2PT_MTH_hit"+to_string(ihit), "P5P5", Sort_light);

     
    PT2_V1[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_2PTH_0_hit"+to_string(ihit), "V1P5", Sort_light);
    PT2_V1_TH[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_2PT_TH_hit"+to_string(ihit), "V1P5", Sort_light);
    PT2_V1_MTH[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_2PT_MTH_hit"+to_string(ihit), "V1P5", Sort_light);

     
     
    PT2_V2[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_2PTH_0_hit"+to_string(ihit), "V2P5", Sort_light);
    PT2_V2_TH[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_2PT_TH_hit"+to_string(ihit), "V2P5", Sort_light);
    PT2_V2_MTH[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_2PT_MTH_hit"+to_string(ihit), "V2P5", Sort_light);

     
     
    BUB_V1_R0[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_OET_ls_R0_hit"+to_string(ihit), "A1P5", Sort_light);
    BUB_V1_R0_TH[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_OET_ls_R0_MPH_hit"+to_string(ihit), "A1P5", Sort_light);
    BUB_V2_R0[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_OET_ls_R0_hit"+to_string(ihit), "A2P5", Sort_light);
    BUB_V2_R0_TH[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_OET_ls_R0_MPH_hit"+to_string(ihit), "A2P5", Sort_light);

    BUB_V1_R1[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_OET_ls_R1_hit"+to_string(ihit), "A1P5", Sort_light);
    BUB_V1_R1_TH[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_OET_ls_R1_PH_hit"+to_string(ihit), "A1P5", Sort_light);
    BUB_V2_R1[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_OET_ls_R1_hit"+to_string(ihit), "A2P5", Sort_light);
    BUB_V2_R1_TH[ihit].Read("../electrounquenching_Klnugamma/ens", "mes_contr_OET_ls_R1_PH_hit"+to_string(ihit), "A2P5", Sort_light);
    
     

  }


  //analyze data

  int iens=0;

  const double ms = 0.01827820;
  const double ml = 0.00072;
  const double fm_to_inv_Gev= 1.0/0.197327;

  string Ens=BUB_V1_R0[0].Tag[iens];
  cout<<"Analyzing ensemble: "<<Ens<<endl;

  boost::filesystem::create_directory("../data/K3lnu/ELECTRO_UNQUENCHING");
  boost::filesystem::create_directory("../data/K3lnu/ELECTRO_UNQUENCHING/"+Ens);

  //read confs list
  vector<string> confs_list;
  ifstream READ_CONFS("../electrounquenching_Klnugamma/ens/cB211b.072.96/confs_list") ;
  while(!READ_CONFS.eof()) {
    string a;
    READ_CONFS >> a;
    if(!READ_CONFS.eof()) confs_list.push_back(a);
  }
  READ_CONFS.close();
  //read source position list
  vector<vector<int>> tt(confs_list.size());
  vector<vector<int>> tz(confs_list.size());

  for(unsigned int i=0;i<confs_list.size();i++) {
    ifstream READ_POS("../electrounquenching_Klnugamma/ens/cB211b.072.96/"+confs_list[i]+"/source_pos");
    while(!READ_POS.eof()) {

      int t, z;
      READ_POS >> t >> z;
      if(!READ_POS.eof()) {
	tt[i].push_back(t); tz[i].push_back(z);
      }
    }
    READ_POS.close();
  }
   
 

  distr_t a_distr(UseJack);

  double Zv= 0.70637654;
  double Za= 0.74278317;
  int Nconfs = BUB_V1_R0[0].Nconfs[0];
   
  GaussianMersenne GM(43112);

  double aB_ave= 0.07948*fm_to_inv_Gev;
  double aB_err =0.00011*fm_to_inv_Gev;

  CorrAnalysis Corr(UseJack,Njacks,100);
  Corr.Nt = BUB_V1_R0[0].nrows[iens];
  int T= Corr.Nt;
  int L= Corr.Nt/2;
  double k= 2*M_PI/L;


  VVVfloat BUB_OET_V1(Nhits), BUB_OET_V2(Nhits);
  VVVfloat BUB_OET_V1_MOM_RE(Nhits), BUB_OET_V2_MOM_RE(Nhits);
  VVVfloat BUB_OET_V1_MOM_IM(Nhits), BUB_OET_V2_MOM_IM(Nhits);

  for(int ihit=0;ihit<Nhits;ihit++) {

    BUB_OET_V1[ihit] = Multiply_Vvector_by_scalar( summ_master( BUB_V1_R0[ihit].col(1)[0], BUB_V1_R1[ihit].col(1)[0]), ms+ml);
    BUB_OET_V1_MOM_RE[ihit] = Multiply_Vvector_by_scalar( summ_master( Multiply_Vvector_by_scalar(BUB_V1_R0_TH[ihit].col(0)[0],-1.0), BUB_V1_R1_TH[ihit].col(0)[0]), ms+ml);
    BUB_OET_V1_MOM_IM[ihit] = Multiply_Vvector_by_scalar( summ_master( BUB_V1_R0_TH[ihit].col(1)[0], BUB_V1_R1_TH[ihit].col(1)[0]), ms+ml);
     
  }

  vector<vector<vector<complex<double>>>> BUB_OET_V1_MOM_Complex(Nhits), BUB_OET_V2_MOM_Complex(Nhits);
  vector<vector<vector<complex<double>>>> PT2_V1_MOM_Complex(Nhits), PT2_V2_MOM_Complex(Nhits);

  for(int ihit=0;ihit<Nhits;ihit++) {
    BUB_OET_V1_MOM_Complex[ihit].resize(T); BUB_OET_V2_MOM_Complex[ihit].resize(T);
    PT2_V1_MOM_Complex[ihit].resize(T); PT2_V2_MOM_Complex[ihit].resize(T);
    for(int t=0;t<T;t++) {
      for(int iconf=0;iconf<Nconfs;iconf++) {
	BUB_OET_V1_MOM_Complex[ihit][t].push_back( BUB_OET_V1_MOM_RE[ihit][t][iconf] + 1i*BUB_OET_V1_MOM_IM[ihit][t][iconf]);
	BUB_OET_V2_MOM_Complex[ihit][t].push_back( BUB_OET_V2_MOM_RE[ihit][t][iconf] + 1i*BUB_OET_V2_MOM_IM[ihit][t][iconf]);

	PT2_V1_MOM_Complex[ihit][t].push_back( PT2_V1_TH[ihit].col(0)[0][t][iconf] + 1i*PT2_V1_TH[ihit].col(0)[0][t][iconf]);
	PT2_V2_MOM_Complex[ihit][t].push_back( PT2_V2_TH[ihit].col(0)[0][t][iconf] + 1i*PT2_V2_TH[ihit].col(0)[0][t][iconf]);
	 
      }
    }
  }

     

  //convolute

  double FAC= pow(L,3)*T*Zv*Za*qd;
   
   
   

  VVVfloat FV(Nhits);

  for(int ihit=0;ihit<Nhits;ihit++) {


    for(int tw=0; tw < Corr.Nt/2; tw++) {


      vector<complex<double>> F(Nconfs,0.0);
      for(int jhit=0;jhit<Nhits;jhit++) {

	 
	 
	for(int t=0;t<Corr.Nt/2;tw++) {
	  for(int iconf=0;iconf<Nconfs;iconf++) {

	    int tj= tt[iconf][jhit];
	    int zj= tz[iconf][jhit];
	    int ti= tt[iconf][ihit];
	    int zi= tz[iconf][ihit];

	    int off_t=(tj-ti+T);
	    
	    F[iconf] += (FAC/k)*(0.5/Nhits)*(BUB_OET_V1_MOM_Complex[jhit][(t+off_t)%T][iconf]*PT2_V2_MOM_Complex[ihit][tw][iconf] -BUB_OET_V2_MOM_Complex[jhit][(t+off_t)%T][iconf]*PT2_V1_MOM_Complex[ihit][tw][iconf] - (BUB_OET_V1[jhit][(t+off_t)%T][iconf]*PT2_V2[ihit].col(0)[0][tw][iconf] -BUB_OET_V2[jhit][(t+off_t)%T][iconf]*PT2_V1[ihit].col(0)[0][tw][iconf]))*exp(k*(t-tw))*exp(1i*(k*zj))*exp(1i*(-1.0*k*zi));
	  }
	}
      }
       

    }

     
  }

  

  


  return;
}



void Analyze_Aprime() {

  
  auto Sort_light_confs = [](string A, string B) {

			   

    int conf_length_A= A.length();
    int conf_length_B= B.length();

    int pos_a_slash=-1;
    int pos_b_slash=-1;
    for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
    for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;

    string A_bis= A.substr(pos_a_slash+1);
    string B_bis= B.substr(pos_b_slash+1);

    //A_bis=A;
    //B_bis=B;

			     
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




  bool Get_ASCII_Mix= true;
  
  if(Get_ASCII_Mix) {
    //read binary files                                                                                                                                                                                                                                                          
    boost::filesystem::create_directory("../Aprime");
    

    vector<string> Ens_T1({"cZ.85.56"});
    vector<string> Ens_TT1({ "cZ211a.085.56"});
      
    
      
    for( int it=0; it<(signed)Ens_T1.size(); it++) {
	
	

      vector<string> channels({"mix_l1_l1", "mix_l2_l2"});
	
	
	
	
      for(auto &channel : channels) {
	boost::filesystem::create_directory("../Aprime/"+channel);
	boost::filesystem::create_directory("../Aprime/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary                                                                                                                                                                                                                                                              
      vector<string> Corr_tags({ "TM_P5P5"});
	
	
      for(int id=0; id<(signed)Corr_tags.size(); id++) {
	  
	  
	  
	for( auto &channel: channels) {
	    
	    
	  vector<string> C_list;
	  ifstream CONFS_LIST;
	    
	  CONFS_LIST.open("../Aprime_bin/"+Ens_T1[it]+"/confsListMix");
	  while(!CONFS_LIST.eof()) {
	    string a;
	    CONFS_LIST>>a;
	    if(!CONFS_LIST.eof()) C_list.push_back(a);
	  }
	    
	    
	  CONFS_LIST.close();
	    
	  FILE *stream = fopen( ("../Aprime_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
	  size_t Nconfs, T, Nhits;
	  bin_read(Nconfs, stream);
	  bin_read(Nhits, stream);
	  bin_read(T, stream);
	  cout<<"Nconfs: "<<Nconfs<<endl;
	  if(Nconfs != C_list.size()) crash("number of configs different than expected: "+to_string(C_list.size()));
	  cout<<"T: "<<T<<" "<<T/2+1<<endl;
	  cout<<"Nhits: "<<Nhits<<endl;
	  for(size_t iconf=0;iconf<Nconfs;iconf++) {
	    vector<double> C(T/2+1);
	    for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	    boost::filesystem::create_directory("../Aprime/"+channel+"/"+Ens_TT1[it]+"/"+C_list[iconf]);
	    ofstream PrintCorr("../Aprime/"+channel+"/"+Ens_TT1[it]+"/"+C_list[iconf]+"/mes_contr_"+channel+"_"+Corr_tags[id]);
	    PrintCorr.precision(16);
	    PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	    for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	    if(Corr_tags[id].substr(3,4) == "V0P5" || Corr_tags[id].substr(3,4) == "TKVK") { for(size_t t=T/2+1; t<T;t++) PrintCorr<<-1*C[T-t]<<endl;   }
	    else  {for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;}
	    PrintCorr.close();
	      
	  }
	    
	  fclose(stream);
	    
	}
	  
      }
    }
  }


  data_t P5P5_l1, P5P5_l2;


  P5P5_l1.Read("../Aprime/mix_l1_l1", "mes_contr_mix_l1_l1_TM_P5P5", "P5P5", Sort_light_confs);
  P5P5_l2.Read("../Aprime/mix_l2_l2", "mes_contr_mix_l2_l2_TM_P5P5", "P5P5", Sort_light_confs);



  int iens=0;


  CorrAnalysis Corr(UseJack, Njacks,Nboots);
  Corr.Nt = P5P5_l1.nrows[iens];
  int L= Corr.Nt/2;
  int T= Corr.Nt;

  Corr.Tmin=30;
  Corr.Tmax=45;


  boost::filesystem::create_directory("../data/Aprime");

  distr_t_list P5P5_l1_distr = Corr.corr_t(P5P5_l1.col(0)[iens], "");
    
  distr_t_list P5P5_l2_distr = Corr.corr_t(P5P5_l2.col(0)[iens], "");


  distr_t_list Mpi_l1 = Corr.effective_mass_t(P5P5_l1_distr, "../data/Aprime/Mpi_l1");
  distr_t_list Mpi_l2 = Corr.effective_mass_t(P5P5_l2_distr, "../data/Aprime/Mpi_l2");

  distr_t_list fpi_l1 = Corr.decay_constant_t( P5P5_l1_distr*pow(2*0.00077,2), "../data/Aprime/fpi_l1");
  distr_t_list fpi_l2 = Corr.decay_constant_t( P5P5_l2_distr*pow(2*0.00085,2), "../data/Aprime/fpi_l2");



  distr_t Mpi_l1_fit= Corr.Fit_distr(Mpi_l1);
  distr_t Mpi_l2_fit= Corr.Fit_distr(Mpi_l2);

  distr_t fpi_l1_fit= Corr.Fit_distr(fpi_l1);
  distr_t fpi_l2_fit= Corr.Fit_distr(fpi_l2);


  cout<<"Mpi(0.00077): "<<Mpi_l1_fit.ave()<<" "<<Mpi_l1_fit.err()<<endl;
  cout<<"Mpi(0.00085): "<<Mpi_l2_fit.ave()<<" "<<Mpi_l2_fit.err()<<endl;
  cout<<"Mpi(expected from ChPT): "<<Mpi_l2_fit.ave()*sqrt(0.00077/0.00085)<<" "<<Mpi_l2_fit.err()*sqrt(0.00077/0.00085)<<endl;

  cout<<"fpi(0.00077): "<<fpi_l1_fit.ave()<<" "<<fpi_l1_fit.err()<<endl;
  cout<<"fpi(0.00085): "<<fpi_l2_fit.ave()<<" "<<fpi_l2_fit.err()<<endl;
    
    



  exit(-1);


  return;
}


void Do_HLT_virtual() {



  

  
  //read data for H1 on B48, B64 and B96 ensembles

  int Nj=50;
  double fm_to_inv_Gev=1.0/0.197327;
  bool Skip_HLT=false;


  //generate lattice spacings and RCs
  //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
  GaussianMersenne GM(36551294);
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack), a_Z(UseJack), a_E(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack), ZV_E(UseJack);
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack), ZA_E(UseJack);
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err, a_Z_ave, a_Z_err, a_E_ave, a_E_err;
  double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err, ZV_E_ave, ZV_E_err;
  double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err, ZA_E_ave, ZA_E_err;
  a_info.LatInfo_new_ens("cA211a.12.48");
  a_A_ave= a_info.a_from_afp_FLAG;
  a_A_err= a_info.a_from_afp_FLAG_err;
  ZA_A_ave = a_info.Za_WI_strange;
  ZA_A_err = a_info.Za_WI_strange_err;
  ZV_A_ave = a_info.Zv_WI_strange;
  ZV_A_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp_FLAG;
  a_B_err= a_info.a_from_afp_FLAG_err;
  ZA_B_ave = a_info.Za_WI_FLAG;
  ZA_B_err = a_info.Za_WI_FLAG_err;
  ZV_B_ave = a_info.Zv_WI_FLAG;
  ZV_B_err = a_info.Zv_WI_FLAG_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp_FLAG;
  a_C_err= a_info.a_from_afp_FLAG_err;
  ZA_C_ave = a_info.Za_WI_FLAG;
  ZA_C_err = a_info.Za_WI_FLAG_err;
  ZV_C_ave = a_info.Zv_WI_FLAG;
  ZV_C_err = a_info.Zv_WI_FLAG_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp_FLAG;
  a_D_err= a_info.a_from_afp_FLAG_err;
  ZA_D_ave = a_info.Za_WI_FLAG;
  ZA_D_err = a_info.Za_WI_FLAG_err;
  ZV_D_ave = a_info.Zv_WI_FLAG;
  ZV_D_err = a_info.Zv_WI_FLAG_err;
  a_info.LatInfo_new_ens("cZ211a.077.64");
  a_Z_ave= a_info.a_from_afp_FLAG;
  a_Z_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cE211a.044.112");
  a_E_ave= a_info.a_from_afp_FLAG;
  a_E_err= a_info.a_from_afp_FLAG_err;
  ZA_E_ave = a_info.Za_WI_FLAG;
  ZA_E_err = a_info.Za_WI_FLAG_err;
  ZV_E_ave = a_info.Zv_WI_FLAG;
  ZV_E_err = a_info.Zv_WI_FLAG_err;


  if(UseJack) {
   
    for(int ijack=0;ijack<Nj;ijack++) {
     
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Nj-1.0))));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Nj-1.0))));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Nj-1.0))));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Nj-1.0))));
      a_Z.distr.push_back( fm_to_inv_Gev*( a_Z_ave + GM()*a_Z_err*(1.0/sqrt(Nj-1.0))));
      a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err*(1.0/sqrt(Nj-1.0))));
      ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err*(1.0/sqrt(Nj -1.0)));
      ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err*(1.0/sqrt(Nj -1.0)));
      ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err*(1.0/sqrt(Nj -1.0)));
      ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err*(1.0/sqrt(Nj -1.0)));
      ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err*(1.0/sqrt(Nj -1.0)));
      ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err*(1.0/sqrt(Nj -1.0)));
      ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err*(1.0/sqrt(Nj -1.0)));
      ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err*(1.0/sqrt(Nj -1.0)));
      ZA_E.distr.push_back(  ZA_E_ave + GM()*ZA_E_err*(1.0/sqrt(Nj -1.0)));
      ZV_E.distr.push_back(  ZV_E_ave + GM()*ZV_E_err*(1.0/sqrt(Nj -1.0)));
       
    }
  }

  else {

    for(int ijack=0;ijack<Nj;ijack++) {
     
      a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err));
      a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err));
      a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err));
      a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err));
      a_Z.distr.push_back( fm_to_inv_Gev*( a_Z_ave + GM()*a_Z_err));
      a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err));
      ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err);
      ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err);
      ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err);
      ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err);
      ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err);
      ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err);
      ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err);
      ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err);
      ZA_E.distr.push_back(  ZA_E_ave + GM()*ZA_E_err);
      ZV_E.distr.push_back(  ZV_E_ave + GM()*ZV_E_err);
	
    }
  }

    




  //###########################################################################################################à




  

  vector<string> Ens_Tag{"B48", "B64", "B96", "C80", "D96"};
  vector<int> L_list{48,64,96,80,96};
  
  vector<int> tw_list{25,25,28,33,39};
  vector<string> xk_list{"60","65","70","75","80","85","90","95"};
  vector<double> xk_list_double{0.6,0.65,0.70,0.75,0.80,0.85,0.9,0.95};
  vector<string> xgamma_list{"0.1","0.3","0.5","0.7"};
  //vector<string> xgamma_list{"0.5"};

  int Nens=Ens_Tag.size();

  //################################# READ DATA #########################################//

  vector<vector<vector<distr_t_list>>> C_H1(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> C_H1_boots(Ens_Tag.size());
  

  //initialize vectors where to store the final results

  vector<vector<vector<distr_t_list>>> H1_RE(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> H1_IM(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> H1_RE_odg(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> H1_IM_odg(Ens_Tag.size());

  vector<vector<vector<distr_t_list>>> H1_RE_naive(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> H1_IM_naive(Ens_Tag.size());

 

  
  int iens=-1;
  for(auto &c: C_H1) {
    iens++;
    string Ens=Ens_Tag[iens];
    int id_xg=-1;
    int Ntimes;
    if(Ens.substr(0,1) == "B") Ntimes=40;
    else if(Ens.substr(0,1) =="C") Ntimes=47;
    else Ntimes=56;
    c.resize(xgamma_list.size());
    for(auto &cc: c) {
      id_xg++;
      for(int i=0;i<(signed)xk_list.size();i++) {
	cc.emplace_back(UseJack);
	for(int t=0;t<Ntimes;t++) {
	  cc[i].distr_list.emplace_back( UseJack,  Read_From_File("../HLT_data_virtual/HLT_data/"+Ens+"_tw_"+to_string(tw_list[iens])+"_dH1_xg_"+xgamma_list[id_xg]+"_xk_"+xk_list[i]+"_tgamma_"+to_string(t)+".dat", 0,1)); 
	}
      }
    }

  }


  iens=-1;
  for(auto &c: C_H1_boots) {
    iens++;
    string Ens=Ens_Tag[iens];
    int id_xg=-1;
    int Ntimes;
    if(Ens.substr(0,1) == "B") Ntimes=40;
    else if(Ens.substr(0,1) =="C") Ntimes=47;
    else Ntimes=56;
    c.resize(xgamma_list.size());
    for(auto &cc: c) {
      id_xg++;
      for(int i=0;i<(signed)xk_list.size();i++) {
	cc.emplace_back(UseJack);
	for(int t=0;t<Ntimes;t++) {
	  cc[i].distr_list.emplace_back( UseJack,  Read_From_File("../HLT_data_virtual/HLT_data_boots/"+Ens+"_tw_"+to_string(tw_list[iens])+"_dH1_xg_"+xgamma_list[id_xg]+"_xk_"+xk_list[i]+"_tgamma_"+to_string(t)+".dat", 0,1)); 
	}
      }
    }

  }

  //#################################

  vector<vector<vector<distr_t_list>>> C_FA(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> C_FA_boots(Ens_Tag.size());
  

  //initialize vectors where to store the final results

  vector<vector<vector<distr_t_list>>> FA_RE(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> FA_IM(Ens_Tag.size());

  vector<vector<vector<distr_t_list>>> FA_RE_odg(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> FA_IM_odg(Ens_Tag.size());

  vector<vector<vector<distr_t_list>>> FA_RE_naive(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> FA_IM_naive(Ens_Tag.size());

  
  iens=-1;
  for(auto &c: C_FA) {
    iens++;
    string Ens=Ens_Tag[iens];
    int id_xg=-1;
    int Ntimes;
    if(Ens.substr(0,1) == "B") Ntimes=40;
    else if(Ens.substr(0,1) =="C") Ntimes=47;
    else Ntimes=56;
    c.resize(xgamma_list.size());
    for(auto &cc: c) {
      id_xg++;
      for(int i=0;i<(signed)xk_list.size();i++) {
	cc.emplace_back(UseJack);
	for(int t=0;t<Ntimes;t++) {
	  cc[i].distr_list.emplace_back( UseJack,  Read_From_File("../HLT_data_virtual/HLT_data/"+Ens+"_tw_"+to_string(tw_list[iens])+"_dFA_xg_"+xgamma_list[id_xg]+"_xk_"+xk_list[i]+"_tgamma_"+to_string(t)+".dat", 0,1)); 
	}
      }
    }

  }


  iens=-1;
  for(auto &c: C_FA_boots) {
    iens++;
    string Ens=Ens_Tag[iens];
    int id_xg=-1;
    int Ntimes;
    if(Ens.substr(0,1) == "B") Ntimes=40;
    else if(Ens.substr(0,1) =="C") Ntimes=47;
    else Ntimes=56;
    c.resize(xgamma_list.size());
    for(auto &cc: c) {
      id_xg++;
      for(int i=0;i<(signed)xk_list.size();i++) {
	cc.emplace_back(UseJack);
	for(int t=0;t<Ntimes;t++) {
	  cc[i].distr_list.emplace_back( UseJack,  Read_From_File("../HLT_data_virtual/HLT_data_boots/"+Ens+"_tw_"+to_string(tw_list[iens])+"_dFA_xg_"+xgamma_list[id_xg]+"_xk_"+xk_list[i]+"_tgamma_"+to_string(t)+".dat", 0,1)); 
	}
      }
    }

  }

  //#############################

  vector<vector<vector<distr_t_list>>> C_H2(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> C_H2_boots(Ens_Tag.size());
  

  //initialize vectors where to store the final results

  vector<vector<vector<distr_t_list>>> H2_RE(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> H2_IM(Ens_Tag.size());

  vector<vector<vector<distr_t_list>>> H2_RE_odg(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> H2_IM_odg(Ens_Tag.size());

  vector<vector<vector<distr_t_list>>> H2_RE_naive(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> H2_IM_naive(Ens_Tag.size());
  
  iens=-1;
  for(auto &c: C_H2) {
    iens++;
    string Ens=Ens_Tag[iens];
    int id_xg=-1;
    int Ntimes;
    if(Ens.substr(0,1) == "B") Ntimes=40;
    else if(Ens.substr(0,1) =="C") Ntimes=47;
    else Ntimes=56;
    c.resize(xgamma_list.size());
    for(auto &cc: c) {
      id_xg++;
      for(int i=0;i<(signed)xk_list.size();i++) {
	cc.emplace_back(UseJack);
	for(int t=0;t<Ntimes;t++) {
	  cc[i].distr_list.emplace_back( UseJack,  Read_From_File("../HLT_data_virtual/HLT_data/"+Ens+"_tw_"+to_string(tw_list[iens])+"_dH2_xg_"+xgamma_list[id_xg]+"_xk_"+xk_list[i]+"_tgamma_"+to_string(t)+".dat", 0,1)); 
	}
      }
    }

  }


  iens=-1;
  for(auto &c: C_H2_boots) {
    iens++;
    string Ens=Ens_Tag[iens];
    int id_xg=-1;
    int Ntimes;
    if(Ens.substr(0,1) == "B") Ntimes=40;
    else if(Ens.substr(0,1) =="C") Ntimes=47;
    else Ntimes=56;
    c.resize(xgamma_list.size());
    for(auto &cc: c) {
      id_xg++;
      for(int i=0;i<(signed)xk_list.size();i++) {
	cc.emplace_back(UseJack);
	for(int t=0;t<Ntimes;t++) {
	  cc[i].distr_list.emplace_back( UseJack,  Read_From_File("../HLT_data_virtual/HLT_data_boots/"+Ens+"_tw_"+to_string(tw_list[iens])+"_dH2_xg_"+xgamma_list[id_xg]+"_xk_"+xk_list[i]+"_tgamma_"+to_string(t)+".dat", 0,1)); 
	}
      }
    }

  }


  //#############################

  vector<vector<vector<distr_t_list>>> C_FV(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> C_FV_boots(Ens_Tag.size());
  

  //initialize vectors where to store the final results

  vector<vector<vector<distr_t_list>>> FV_RE(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> FV_IM(Ens_Tag.size());

  vector<vector<vector<distr_t_list>>> FV_RE_odg(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> FV_IM_odg(Ens_Tag.size());

  vector<vector<vector<distr_t_list>>> FV_RE_naive(Ens_Tag.size());
  vector<vector<vector<distr_t_list>>> FV_IM_naive(Ens_Tag.size());

  
  iens=-1;
  for(auto &c: C_FV) {
    iens++;
    string Ens=Ens_Tag[iens];
    int id_xg=-1;
    int Ntimes;
    if(Ens.substr(0,1) == "B") Ntimes=40;
    else if(Ens.substr(0,1) =="C") Ntimes=47;
    else Ntimes=56;
    c.resize(xgamma_list.size());
    for(auto &cc: c) {
      id_xg++;
      for(int i=0;i<(signed)xk_list.size();i++) {
	cc.emplace_back(UseJack);
	for(int t=0;t<Ntimes;t++) {
	  cc[i].distr_list.emplace_back( UseJack,  Read_From_File("../HLT_data_virtual/HLT_data/"+Ens+"_tw_"+to_string(tw_list[iens])+"_dFV_xg_"+xgamma_list[id_xg]+"_xk_"+xk_list[i]+"_tgamma_"+to_string(t)+".dat", 0,1)); 
	}
      }
    }

  }


  iens=-1;
  for(auto &c: C_FV_boots) {
    iens++;
    string Ens=Ens_Tag[iens];
    int id_xg=-1;
    int Ntimes;
    if(Ens.substr(0,1) == "B") Ntimes=40;
    else if(Ens.substr(0,1) =="C") Ntimes=47;
    else Ntimes=56;
    c.resize(xgamma_list.size());
    for(auto &cc: c) {
      id_xg++;
      for(int i=0;i<(signed)xk_list.size();i++) {
	cc.emplace_back(UseJack);
	for(int t=0;t<Ntimes;t++) {
	  cc[i].distr_list.emplace_back( UseJack,  Read_From_File("../HLT_data_virtual/HLT_data_boots/"+Ens+"_tw_"+to_string(tw_list[iens])+"_dFV_xg_"+xgamma_list[id_xg]+"_xk_"+xk_list[i]+"_tgamma_"+to_string(t)+".dat", 0,1)); 
	}
      }
    }

  }



  

  cout<<"Data read"<<endl;

  

  vector<double> eps_list { 0.1, 0.125, 0.15, 0.175, 0.20, 0.25, 0.3, 0.35}; //GeV


  for(int iens=0;iens<Nens;iens++) {


    double Mpi=0.0;

    string Ens= Ens_Tag[iens];

    distr_t a_distr(UseJack);
    
    if(Ens.substr(0,1) == "B") { a_distr = a_B; Mpi=0.1402;}
    else if(Ens.substr(0,1) == "C") { a_distr= a_C; Mpi = 0.1367; }
    else if(Ens.substr(0,1) == "D") { a_distr= a_D; Mpi=0.141; } 
    else crash("Ensemble: "+Ens+" not found");

    int L= L_list[iens];

    H1_RE[iens].resize(xgamma_list.size());
    H1_IM[iens].resize(xgamma_list.size());
    H1_RE_odg[iens].resize(xgamma_list.size());
    H1_IM_odg[iens].resize(xgamma_list.size());
    H1_RE_naive[iens].resize(xgamma_list.size());
    H1_IM_naive[iens].resize(xgamma_list.size());

    H2_RE[iens].resize(xgamma_list.size());
    H2_IM[iens].resize(xgamma_list.size());
    H2_RE_odg[iens].resize(xgamma_list.size());
    H2_IM_odg[iens].resize(xgamma_list.size());
    H2_RE_naive[iens].resize(xgamma_list.size());
    H2_IM_naive[iens].resize(xgamma_list.size());

    FA_RE[iens].resize(xgamma_list.size());
    FA_IM[iens].resize(xgamma_list.size());
    FA_RE_odg[iens].resize(xgamma_list.size());
    FA_IM_odg[iens].resize(xgamma_list.size());
    FA_RE_naive[iens].resize(xgamma_list.size());
    FA_IM_naive[iens].resize(xgamma_list.size());


    FV_RE[iens].resize(xgamma_list.size());
    FV_IM[iens].resize(xgamma_list.size());
    FV_RE_odg[iens].resize(xgamma_list.size());
    FV_IM_odg[iens].resize(xgamma_list.size());
    FV_RE_naive[iens].resize(xgamma_list.size());
    FV_IM_naive[iens].resize(xgamma_list.size());

    
    
    
    for(int xg=0;xg<(signed)xgamma_list.size() ; xg++) {

      Vfloat xk_list_tailored;
      if(xgamma_list[xg]=="0.5") {xk_list_tailored = {0.6,0.65,0.7,0.75} ; }
      else if(xgamma_list[xg]=="0.7") {xk_list_tailored = {0.6} ; }
      else if(xgamma_list[xg]=="0.3") {xk_list_tailored = {0.6,0.65,0.7,0.75,0.8,0.85};}
      else if(xgamma_list[xg]=="0.1") {xk_list_tailored = {0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95} ; }
      else crash("xgamma value: "+xgamma_list[xg]+" not allowed!");

      
		 

      for(int xk=0;xk<(signed)xk_list_tailored.size() ; xk++) {


	H1_RE[iens][xg].emplace_back( UseJack, eps_list.size());
	H1_IM[iens][xg].emplace_back( UseJack, eps_list.size());
	H1_RE_odg[iens][xg].emplace_back( UseJack, eps_list.size());
	H1_IM_odg[iens][xg].emplace_back( UseJack, eps_list.size());
	H1_RE_naive[iens][xg].emplace_back( UseJack, eps_list.size());
	H1_IM_naive[iens][xg].emplace_back( UseJack, eps_list.size());
	
	H2_RE[iens][xg].emplace_back( UseJack, eps_list.size());
	H2_IM[iens][xg].emplace_back( UseJack, eps_list.size());
	H2_RE_odg[iens][xg].emplace_back( UseJack, eps_list.size());
	H2_IM_odg[iens][xg].emplace_back( UseJack, eps_list.size());
	H2_RE_naive[iens][xg].emplace_back( UseJack, eps_list.size());
	H2_IM_naive[iens][xg].emplace_back( UseJack, eps_list.size());
	
	FA_RE[iens][xg].emplace_back( UseJack, eps_list.size());
	FA_IM[iens][xg].emplace_back( UseJack, eps_list.size());
	FA_RE_odg[iens][xg].emplace_back( UseJack, eps_list.size());
	FA_IM_odg[iens][xg].emplace_back( UseJack, eps_list.size());
	FA_RE_naive[iens][xg].emplace_back( UseJack, eps_list.size());
	FA_IM_naive[iens][xg].emplace_back( UseJack, eps_list.size());
	
	FV_RE[iens][xg].emplace_back( UseJack, eps_list.size());
	FV_IM[iens][xg].emplace_back( UseJack, eps_list.size());
	FV_RE_odg[iens][xg].emplace_back( UseJack, eps_list.size());
	FV_IM_odg[iens][xg].emplace_back( UseJack, eps_list.size());
	FV_RE_naive[iens][xg].emplace_back( UseJack, eps_list.size());
	FV_IM_naive[iens][xg].emplace_back( UseJack, eps_list.size());
	
	distr_t_list C =  C_H1[iens][xg][xk] ;
	distr_t_list C_2 = C_H2[iens][xg][xk];
	distr_t_list C_A = C_FA[iens][xg][xk];
	distr_t_list C_V = C_FV[iens][xg][xk];
	distr_t_list C_boots = C_H1_boots[iens][xg][xk];
	distr_t_list C_2_boots = C_H2_boots[iens][xg][xk];
	distr_t_list C_A_boots = C_FA_boots[iens][xg][xk];
	distr_t_list C_V_boots = C_FV_boots[iens][xg][xk];

	


	boost::filesystem::create_directory("../data/HLT_virtual");
	boost::filesystem::create_directory("../data/HLT_virtual/"+Ens_Tag[iens]);
	boost::filesystem::create_directory("../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]);
	Print_To_File({},{C.ave(), C.err()}, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/C_H1_xk_"+to_string_with_precision(xk_list_tailored[xk],3),"","");
	Print_To_File({},{C_2.ave(), C_2.err()}, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/C_H2_xk_"+to_string_with_precision(xk_list_tailored[xk],3),"","");
	Print_To_File({},{C_A.ave(), C_A.err()}, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/C_FA_xk_"+to_string_with_precision(xk_list_tailored[xk],3),"","");
	Print_To_File({},{C_V.ave(), C_V.err()}, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/C_FV_xk_"+to_string_with_precision(xk_list_tailored[xk],3),"","");

	//compute covariance matrix (to be updated by loading the correlation matrix from file)

	Vfloat Corr;	Vfloat Cov;
	Vfloat Corr_2;	Vfloat Cov_2;
	Vfloat Corr_A;	Vfloat Cov_A;
	Vfloat Corr_V;	Vfloat Cov_V;
	int id=0;

	int TT = C_boots.size();
	
	for(int tt=0;tt<TT;tt++) {
	  for(int rr=0;rr<TT;rr++) {
	    Corr.push_back(  C_boots.distr_list[tt]%C_boots.distr_list[rr]/( C.err(tt)*C.err(rr) ));
	    Cov.push_back( Corr[id]*C.err(tt)*C.err(rr) );

	    Corr_2.push_back(  C_2_boots.distr_list[tt]%C_2_boots.distr_list[rr]/( C_2.err(tt)*C_2.err(rr) ));
	    Cov_2.push_back( Corr_2[id]*C_2.err(tt)*C_2.err(rr) );

	    Corr_A.push_back(  C_A_boots.distr_list[tt]%C_A_boots.distr_list[rr]/( C_A.err(tt)*C_A.err(rr) ));
	    Cov_A.push_back( Corr_A[id]*C_A.err(tt)*C_A.err(rr) );

	    Corr_V.push_back(  C_V_boots.distr_list[tt]%C_V_boots.distr_list[rr]/( C_V.err(tt)*C_V.err(rr) ));
	    Cov_V.push_back( Corr_V[id]*C_V.err(tt)*C_V.err(rr) );
	    
	    id++;
	  }
	}

	//do the HLT
	//define lambda functions
	
	auto K_RE= [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

	  PrecFloat x= (E-m)/2;
	  PrecFloat cosh_ov_sinh_half= (exp(x) + exp(-3*x))/(1-exp(-2*x));
	  PrecFloat cosh_ov_cosh_half= (exp(x) + exp(-3*x))/(1+exp(-2*x)); 
	  
	  return  exp(-x)*cos(s/2)*cos(s/2)/(cosh_ov_sinh_half - cos(s)/sinh(x)) -exp(-x)*sin(s/2)*sin(s/2)/(cosh_ov_cosh_half - cos(s)/cosh(x));

	  
	  
	};
	
	auto K_IM = [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {

	    	  
	  PrecFloat x= (E-m)/2;
	  PrecFloat cosh_ov_sinh_half= (exp(x) + exp(-3*x))/(1-exp(-2*x));
	  PrecFloat cosh_ov_cosh_half= (exp(x) + exp(-3*x))/(1+exp(-2*x)); 

	  return exp(-x)*cos(s/2)*sin(s/2)/(cosh_ov_cosh_half - cos(s)/cosh(x)) +  exp(-x)*sin(s/2)*cos(s/2)/(cosh_ov_sinh_half - cos(s)/sinh(x));
		  
	};


	//parallel loop

	if(!Skip_HLT) {
	  

#pragma omp parallel for
	  for(int ieps=0; ieps < (signed)eps_list.size(); ieps++) {
	    
	    double Emin=0;
	    double k= 0.5*MK_FLAG*stod(xgamma_list[xg]);
	    //Emin= 2*sqrt( pow(Mpi,2) + pow(2*M_PI/(L*a_distr.ave()),2));
	    //Emin= 0.93*sqrt( pow(Emin,2) + pow(0.5*MK_FLAG*stod(xgamma_list[xg]),2));
	    double Ea = sqrt( pow(Mpi,2) + pow( 2*M_PI/(L*a_distr.ave()) -k,2)) +  sqrt( pow(Mpi,2) + pow( 2*M_PI/(L*a_distr.ave()),2));
	    double Eb =  sqrt( pow(Mpi,2) + pow( 2*M_PI/(L*a_distr.ave()) -k/2,2)) +  sqrt( pow(Mpi,2) + pow( 2*M_PI/(L*a_distr.ave()) +k/2,2));
	    Emin = 0.95*min(Ea,Eb);
	    Emin = 2*sqrt( pow(Mpi,2) + pow(k/2,2));

	    
	    
	    cout<<"##### XGAMMA: "<<xgamma_list[xg]<<endl;
	    
	    double aE0 = Emin*a_distr.ave();
	    double erg= MK_FLAG*a_distr.ave()*sqrt( pow(xk_list_tailored[xk],2) + pow(0.5*stod(xgamma_list[xg]),2));
	    double erg_GeV = erg/a_distr.ave();
	    double sigma= eps_list[ieps]*a_distr.ave();
	    double sigma_GeV= sigma/a_distr.ave();

	    double E_p0 = ( Mpi + sqrt( pow(Mpi,2) + pow(k,2)))*a_distr.ave();

	    E_p0 += 2*(E_p0-aE0);
	    double E_p1 = 0.95*min(Ea,Eb)*a_distr.ave();

	    int gamma=3;
	    
	    double r= pow(E_p0/E_p1,2)*pow((1-pow(aE0/E_p0,2)),gamma)/pow((1-pow(aE0/E_p1,2)),gamma);
	    	  
	    cout<<"r: "<<r<<endl;
	    
	    if(E_p0 < aE0) crash("E_p0 < aE0");
	    if(E_p1 < E_p0)  { aE0 = E_p1; r=-1.0;}
	    
	    
	    double syst, l;
	    
	    //######### HLT PARS ##########
	    double mult_IM = 5e-5;
	    double mult_RE = 1e-4;

	    mult_IM *= 10.0;
	    mult_RE *= 10.0;
	    
	    double Ag_target= 1e-3;
	    int tmax;
	    if(Ens.substr(0,1)=="B" ) tmax=38;
	    else tmax = (int)(  38.0*0.0795*fm_to_inv_Gev/a_distr.ave() );
	    double Emax=0.0;
	    int Is_Emax_Finite=(Emax > 1e-10)?1:0;
	    double alpha=0.0;
	    int prec=90;
	    //############################
	    PrecFloat::setDefaultPrecision(prec);
	    
	    
	    distr_t C_0 = 0.0*Get_id_jack_distr(Nj);
	    distr_t C_0_2 = 0.0*Get_id_jack_distr(Nj);
	    distr_t C_0_V = 0.0*Get_id_jack_distr(Nj);
	    distr_t C_0_A = 0.0*Get_id_jack_distr(Nj);

	    
	    for(int t=1;t<TT;t++) {
	      C_0 = C_0 + C.distr_list[t];
	      C_0_2 = C_0_2 + C_2.distr_list[t];
	      C_0_V = C_0_V + C_V.distr_list[t];
	      C_0_A = C_0_A + C_A.distr_list[t];
	    }

	    //########## H1
	    distr_t F_RE =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0, E_p0, E_p1, r,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C , syst, mult_RE, l, "TANT", "H1_"+Ens_Tag[iens], "RE", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1) - C_0;
	    F_RE = F_RE.ave() + (F_RE-F_RE.ave())*sqrt( 1.0 + pow(syst/F_RE.err(),2));
	    distr_t F_IM =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0, E_p0, E_p1, r,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C , syst, mult_IM, l, "TANT", "H1_"+Ens_Tag[iens], "IM", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);
	    F_IM = F_IM.ave() + (F_IM-F_IM.ave())*sqrt( 1.0 + pow(syst/F_IM.err(),2));

	    if(Ens.substr(0,1)=="B" ) tmax=30;
	    else tmax = (int)(  30.0*0.0795*fm_to_inv_Gev/a_distr.ave() );
	    
	    //#########  H2
	    distr_t F_2_RE =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0, E_p0, E_p1, r,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C_2 , syst, mult_RE, l, "TANT", "H2_"+Ens_Tag[iens], "RE", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_2, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1) - C_0_2;
	    F_2_RE = F_2_RE.ave() + (F_2_RE-F_2_RE.ave())*sqrt( 1.0 + pow(syst/F_2_RE.err(),2));
	    distr_t F_2_IM =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0,  E_p0, E_p1, r,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C_2 , syst, mult_IM, l, "TANT", "H2_"+Ens_Tag[iens], "IM", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_2, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);
	    F_2_IM = F_2_IM.ave() + (F_2_IM-F_2_IM.ave())*sqrt( 1.0 + pow(syst/F_2_IM.err(),2));
	    //#########  FA
	    distr_t F_A_RE =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0,  E_p0, E_p1, r,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C_A , syst, mult_RE, l, "TANT", "FA_"+Ens_Tag[iens], "RE", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_A, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1) - C_0_A;
	    F_A_RE = F_A_RE.ave() + (F_A_RE-F_A_RE.ave())*sqrt( 1.0 + pow(syst/F_A_RE.err(),2));
	    distr_t F_A_IM =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0,  E_p0, E_p1, r,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C_A , syst, mult_IM, l, "TANT", "FA_"+Ens_Tag[iens], "IM", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_A, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);
	    F_A_IM = F_A_IM.ave() + (F_A_IM-F_A_IM.ave())*sqrt( 1.0 + pow(syst/F_A_IM.err(),2));
	    //#########  FV	    
	    distr_t F_V_RE =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0,  E_p0, E_p1, r,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C_V , syst, mult_RE, l, "TANT", "FV_"+Ens_Tag[iens], "RE", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_V, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1) - C_0_V;
	    F_V_RE = F_V_RE.ave() + (F_V_RE-F_V_RE.ave())*sqrt( 1.0 + pow(syst/F_V_RE.err(),2));
	    distr_t F_V_IM =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0,  E_p0, E_p1, r,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C_V , syst, mult_IM, l, "TANT", "FV_"+Ens_Tag[iens], "IM", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_V, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);
	    F_V_IM = F_V_IM.ave() + (F_V_IM-F_V_IM.ave())*sqrt( 1.0 + pow(syst/F_V_IM.err(),2));
												 

	    
	    
	    //push back the results
	    H1_RE[iens][xg][xk].distr_list[ieps] = F_RE;
	    H1_IM[iens][xg][xk].distr_list[ieps] = F_IM;
	    //push back the results
	    H2_RE[iens][xg][xk].distr_list[ieps] = F_2_RE;
	    H2_IM[iens][xg][xk].distr_list[ieps] = F_2_IM;
	    //push back the results
	    FA_RE[iens][xg][xk].distr_list[ieps] = F_A_RE;
	    FA_IM[iens][xg][xk].distr_list[ieps] = F_A_IM;
	    //push back the results
	    FV_RE[iens][xg][xk].distr_list[ieps] = F_V_RE;
	    FV_IM[iens][xg][xk].distr_list[ieps] = F_V_IM;
	    	    
	  }
#pragma omp barrier


	  
#pragma omp parallel for
	  for(int ieps=0; ieps < (signed)eps_list.size(); ieps++) {
	    
	    double Emin=0;
	    double k= 0.5*MK_FLAG*stod(xgamma_list[xg]);
	    //Emin= 2*sqrt( pow(Mpi,2) + pow(2*M_PI/(L*a_distr.ave()),2));
	    //Emin= 0.93*sqrt( pow(Emin,2) + pow(0.5*MK_FLAG*stod(xgamma_list[xg]),2));
	    double Ea = sqrt( pow(Mpi,2) + pow( 2*M_PI/(L*a_distr.ave()) -k,2)) +  sqrt( pow(Mpi,2) + pow( 2*M_PI/(L*a_distr.ave()),2));
	    double Eb =  sqrt( pow(Mpi,2) + pow( 2*M_PI/(L*a_distr.ave()) -k/2,2)) +  sqrt( pow(Mpi,2) + pow( 2*M_PI/(L*a_distr.ave()) +k/2,2));
	    Emin = 0.95*min(Ea,Eb);
	    Emin = 2*sqrt( pow(Mpi,2) + pow(k/2,2));

	    
	    cout<<"##### XGAMMA: "<<xgamma_list[xg]<<endl;
	    
	    double aE0 = Emin*a_distr.ave();
	    double erg= MK_FLAG*a_distr.ave()*sqrt( pow(xk_list_tailored[xk],2) + pow(0.5*stod(xgamma_list[xg]),2));
	    double erg_GeV = erg/a_distr.ave();
	    double sigma= eps_list[ieps]*a_distr.ave();
	    double sigma_GeV= sigma/a_distr.ave();
	    double E_p0 = ( Mpi + sqrt( pow(Mpi,2) + pow(k,2)))*a_distr.ave();
	    E_p0 += 2*(E_p0-aE0);
	    
	    double E_p1 = 0.95*min(Ea,Eb)*a_distr.ave();

	    double aE0_odg = E_p1;

	    	    
	    double r2 = -1.0;
	    
	    
	    double syst, l;
	    
	    //######### HLT PARS ##########
	    double mult_IM = 5e-5;
	    double mult_RE = 1e-4;

	    mult_IM *= 10.0;
	    mult_RE *= 10.0;
	    
	    double Ag_target= 1e-3;
	    int tmax;
	    if(Ens.substr(0,1)=="B" ) tmax=38;
	    else tmax = (int)(  38.0*0.0795*fm_to_inv_Gev/a_distr.ave() );
	    double Emax=0.0;
	    int Is_Emax_Finite=(Emax > 1e-10)?1:0;
	    double alpha=0.0;
	    int prec=90;
	    //############################
	    PrecFloat::setDefaultPrecision(prec);
	    
	    
	    distr_t C_0 = 0.0*Get_id_jack_distr(Nj);
	    distr_t C_0_2 = 0.0*Get_id_jack_distr(Nj);
	    distr_t C_0_V = 0.0*Get_id_jack_distr(Nj);
	    distr_t C_0_A = 0.0*Get_id_jack_distr(Nj);
	    

	    
	    for(int t=1;t<TT;t++) {
	      C_0 = C_0 + C.distr_list[t];
	      C_0_2 = C_0_2 + C_2.distr_list[t];
	      C_0_V = C_0_V + C_V.distr_list[t];
	      C_0_A = C_0_A + C_A.distr_list[t];
	    }

	    if(Ens.substr(0,1)=="B" ) tmax=38;
	    else tmax = (int)(  38.0*0.0795*fm_to_inv_Gev/a_distr.ave() );

	    
	    //#######  H1
	    distr_t F_RE_odg =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0_odg,  E_p0, E_p1, r2,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C , syst, mult_RE, l, "TANT", "H1_odg_"+Ens_Tag[iens], "RE", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1) - C_0;
	    F_RE_odg = F_RE_odg.ave() + (F_RE_odg-F_RE_odg.ave())*sqrt( 1.0 + pow(syst/F_RE_odg.err(),2));
	    distr_t F_IM_odg =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0_odg,  E_p0, E_p1, r2,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C , syst, mult_IM, l, "TANT", "H1_odg_"+Ens_Tag[iens], "IM", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);
	    F_IM_odg = F_IM_odg.ave() + (F_IM_odg-F_IM_odg.ave())*sqrt( 1.0 + pow(syst/F_IM_odg.err(),2));

	    if(Ens.substr(0,1)=="B" ) tmax=30;
	    else tmax = (int)(  30.0*0.0795*fm_to_inv_Gev/a_distr.ave() );
	    
	    //#########  H2
	    distr_t F_2_RE_odg =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0_odg,  E_p0, E_p1, r2,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C_2 , syst, mult_RE, l, "TANT", "H2_odg_"+Ens_Tag[iens], "RE", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_2, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1) - C_0_2;
	    F_2_RE_odg = F_2_RE_odg.ave() + (F_2_RE_odg-F_2_RE_odg.ave())*sqrt( 1.0 + pow(syst/F_2_RE_odg.err(),2));
	    distr_t F_2_IM_odg =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0_odg,  E_p0, E_p1, r2,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C_2 , syst, mult_IM, l, "TANT", "H2_odg_"+Ens_Tag[iens], "IM", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_2, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha,1);
	    F_2_IM_odg = F_2_IM_odg.ave() + (F_2_IM_odg-F_2_IM_odg.ave())*sqrt( 1.0 + pow(syst/F_2_IM_odg.err(),2));
	    //#########  FA
	    distr_t F_A_RE_odg =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0_odg,  E_p0, E_p1, r2,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C_A , syst, mult_RE, l, "TANT", "FA_odg_"+Ens_Tag[iens], "RE", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_A, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1) - C_0_A;
	    F_A_RE_odg = F_A_RE_odg.ave() + (F_A_RE_odg-F_A_RE_odg.ave())*sqrt( 1.0 + pow(syst/F_A_RE_odg.err(),2));
	    distr_t F_A_IM_odg =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0_odg,  E_p0, E_p1, r2,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C_A , syst, mult_IM, l, "TANT", "FA_odg_"+Ens_Tag[iens], "IM", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_A, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);
	    F_A_IM_odg = F_A_IM_odg.ave() + (F_A_IM_odg-F_A_IM_odg.ave())*sqrt( 1.0 + pow(syst/F_A_IM_odg.err(),2));
	    //#########  FV
	    distr_t F_V_RE_odg =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0_odg,  E_p0, E_p1, r2,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C_V , syst, mult_RE, l, "TANT", "FV_odg_"+Ens_Tag[iens], "RE", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_V, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1) - C_0_V;
	    F_V_RE_odg = F_V_RE_odg.ave() + (F_V_RE_odg-F_V_RE_odg.ave())*sqrt( 1.0 + pow(syst/F_V_RE_odg.err(),2));
	    distr_t F_V_IM_odg =  Get_Laplace_transfo_piecewise(  erg,  sigma, aE0_odg,  E_p0, E_p1, r2,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C_V , syst, mult_IM, l, "TANT", "FV_odg_"+Ens_Tag[iens], "IM", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_V, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);
	    F_V_IM_odg = F_V_IM_odg.ave() + (F_V_IM_odg-F_V_IM_odg.ave())*sqrt( 1.0 + pow(syst/F_V_IM_odg.err(),2));
	    
	    
	    //push back the results
	    H1_RE_odg[iens][xg][xk].distr_list[ieps] = F_RE_odg;
	    H1_IM_odg[iens][xg][xk].distr_list[ieps] = F_IM_odg;

	    //push back the results
	    H2_RE_odg[iens][xg][xk].distr_list[ieps] = F_2_RE_odg;
	    H2_IM_odg[iens][xg][xk].distr_list[ieps] = F_2_IM_odg;

	    //push back the results
	    FA_RE_odg[iens][xg][xk].distr_list[ieps] = F_A_RE_odg;
	    FA_IM_odg[iens][xg][xk].distr_list[ieps] = F_A_IM_odg;

	    //push back the results
	    FV_RE_odg[iens][xg][xk].distr_list[ieps] = F_V_RE_odg;
	    FV_IM_odg[iens][xg][xk].distr_list[ieps] = F_V_IM_odg;
	    
	  }
#pragma omp barrier



#pragma omp parallel for
	  for(int ieps=0; ieps < (signed)eps_list.size(); ieps++) {
	    
	    double Emin=0;
	    double k= 0.5*MK_FLAG*stod(xgamma_list[xg]);
	    //Emin= 2*sqrt( pow(Mpi,2) + pow(2*M_PI/(L*a_distr.ave()),2));
	    //Emin= 0.93*sqrt( pow(Emin,2) + pow(0.5*MK_FLAG*stod(xgamma_list[xg]),2));
	    double Ea = sqrt( pow(Mpi,2) + pow( 2*M_PI/(L*a_distr.ave()) -k,2)) +  sqrt( pow(Mpi,2) + pow( 2*M_PI/(L*a_distr.ave()),2));
	    double Eb =  sqrt( pow(Mpi,2) + pow( 2*M_PI/(L*a_distr.ave()) -k/2,2)) +  sqrt( pow(Mpi,2) + pow( 2*M_PI/(L*a_distr.ave()) +k/2,2));
	    Emin = 0.95*min(Ea,Eb);
	    Emin = 2*sqrt( pow(Mpi,2) + pow(k/2,2));

	    
	    cout<<"##### XGAMMA: "<<xgamma_list[xg]<<endl;

	    double aE0_naive =  ( Mpi + sqrt( pow(Mpi,2) + pow(k,2)))*a_distr.ave();
	    double erg= MK_FLAG*a_distr.ave()*sqrt( pow(xk_list_tailored[xk],2) + pow(0.5*stod(xgamma_list[xg]),2));
	    double erg_GeV = erg/a_distr.ave();
	    double sigma= eps_list[ieps]*a_distr.ave();
	    double sigma_GeV= sigma/a_distr.ave();
	    	    
	    double syst, l;
	    
	    //######### HLT PARS ##########
	    double mult_IM = 5e-5;
	    double mult_RE = 1e-4;

	    mult_IM *= 10.0;
	    mult_RE *= 10.0;
	    
	    double Ag_target= 1e-3;
	    int tmax;
	    if(Ens.substr(0,1)=="B" ) tmax=38;
	    else tmax = (int)(  38.0*0.0795*fm_to_inv_Gev/a_distr.ave() );
	    double Emax=0.0;
	    int Is_Emax_Finite=(Emax > 1e-10)?1:0;
	    double alpha=0.0;
	    int prec=90;
	    //############################
	    PrecFloat::setDefaultPrecision(prec);
	    
	    
	    distr_t C_0 = 0.0*Get_id_jack_distr(Nj);
	    distr_t C_0_2 = 0.0*Get_id_jack_distr(Nj);
	    distr_t C_0_V = 0.0*Get_id_jack_distr(Nj);
	    distr_t C_0_A = 0.0*Get_id_jack_distr(Nj);
	    
	    
	    for(int t=1;t<TT;t++) {
	      C_0 = C_0 + C.distr_list[t];
	      C_0_2 = C_0_2 + C_2.distr_list[t];
	      C_0_V = C_0_V + C_V.distr_list[t];
	      C_0_A = C_0_A + C_A.distr_list[t];
	    }

	    if(Ens.substr(0,1)=="B" ) tmax=38;
	    else tmax = (int)(  38.0*0.0795*fm_to_inv_Gev/a_distr.ave() );

	    
	    //#######  H1
	    distr_t F_RE_naive =  Get_Laplace_transfo(  erg,  sigma, aE0_naive,   TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C , syst, mult_RE, l, "TANT", "H1_naive_"+Ens_Tag[iens], "RE", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1) - C_0;
	    F_RE_naive = F_RE_naive.ave() + (F_RE_naive-F_RE_naive.ave())*sqrt( 1.0 + pow(syst/F_RE_naive.err(),2));
	    distr_t F_IM_naive =  Get_Laplace_transfo(  erg,  sigma, aE0_naive,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C , syst, mult_IM, l, "TANT", "H1_naive_"+Ens_Tag[iens], "IM", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);
	    F_IM_naive = F_IM_naive.ave() + (F_IM_naive-F_IM_naive.ave())*sqrt( 1.0 + pow(syst/F_IM_naive.err(),2));

	    if(Ens.substr(0,1)=="B" ) tmax=30;
	    else tmax = (int)(  30.0*0.0795*fm_to_inv_Gev/a_distr.ave() );
	    
	    //#########  H2
	    distr_t F_2_RE_naive =  Get_Laplace_transfo(  erg,  sigma, aE0_naive,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C_2 , syst, mult_RE, l, "TANT", "H2_naive_"+Ens_Tag[iens], "RE", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_2, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1) - C_0_2;
	    F_2_RE_naive = F_2_RE_naive.ave() + (F_2_RE_naive-F_2_RE_naive.ave())*sqrt( 1.0 + pow(syst/F_2_RE_naive.err(),2));
	    distr_t F_2_IM_naive =  Get_Laplace_transfo(  erg,  sigma, aE0_naive,   TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C_2 , syst, mult_IM, l, "TANT", "H2_naive_"+Ens_Tag[iens], "IM", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_2, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha,1);
	    F_2_IM_naive = F_2_IM_naive.ave() + (F_2_IM_naive-F_2_IM_naive.ave())*sqrt( 1.0 + pow(syst/F_2_IM_naive.err(),2));
	    //#########  FA
	    distr_t F_A_RE_naive =  Get_Laplace_transfo(  erg,  sigma, aE0_naive,   TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C_A , syst, mult_RE, l, "TANT", "FA_naive_"+Ens_Tag[iens], "RE", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_A, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1) - C_0_A;
	    F_A_RE_naive = F_A_RE_naive.ave() + (F_A_RE_naive-F_A_RE_naive.ave())*sqrt( 1.0 + pow(syst/F_A_RE_naive.err(),2));
	    distr_t F_A_IM_naive =  Get_Laplace_transfo(  erg,  sigma, aE0_naive,   TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C_A , syst, mult_IM, l, "TANT", "FA_naive_"+Ens_Tag[iens], "IM", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_A, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);
	    F_A_IM_naive = F_A_IM_naive.ave() + (F_A_IM_naive-F_A_IM_naive.ave())*sqrt( 1.0 + pow(syst/F_A_IM_naive.err(),2));
	    //#########  FV
	    distr_t F_V_RE_naive =  Get_Laplace_transfo(  erg,  sigma, aE0_naive,  TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C_V , syst, mult_RE, l, "TANT", "FV_naive_"+Ens_Tag[iens], "RE", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_V, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1) - C_0_V;
	    F_V_RE_naive = F_V_RE_naive.ave() + (F_V_RE_naive-F_V_RE_naive.ave())*sqrt( 1.0 + pow(syst/F_V_RE_naive.err(),2));
	    distr_t F_V_IM_naive =  Get_Laplace_transfo(  erg,  sigma, aE0_naive,   TT, tmax , prec, "xk_"+to_string_with_precision(xk_list_tailored[xk],3)+"_xg_"+xgamma_list[xg]+"_Erg_"+to_string_with_precision(erg_GeV,3)+"_Emin_"+to_string_with_precision(Emin,3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C_V , syst, mult_IM, l, "TANT", "FV_naive_"+Ens_Tag[iens], "IM", Ag_target,0, Get_id_distr(Nj,UseJack) , 0.0 , "K_virtual", Cov_V, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);
	    F_V_IM_naive = F_V_IM_naive.ave() + (F_V_IM_naive-F_V_IM_naive.ave())*sqrt( 1.0 + pow(syst/F_V_IM_naive.err(),2));
	    
	    
	    //push back the results
	    H1_RE_naive[iens][xg][xk].distr_list[ieps] = F_RE_naive;
	    H1_IM_naive[iens][xg][xk].distr_list[ieps] = F_IM_naive;

	    //push back the results
	    H2_RE_naive[iens][xg][xk].distr_list[ieps] = F_2_RE_naive;
	    H2_IM_naive[iens][xg][xk].distr_list[ieps] = F_2_IM_naive;

	    //push back the results
	    FA_RE_naive[iens][xg][xk].distr_list[ieps] = F_A_RE_naive;
	    FA_IM_naive[iens][xg][xk].distr_list[ieps] = F_A_IM_naive;

	    //push back the results
	    FV_RE_naive[iens][xg][xk].distr_list[ieps] = F_V_RE_naive;
	    FV_IM_naive[iens][xg][xk].distr_list[ieps] = F_V_IM_naive;
	    
	  }
#pragma omp barrier
	
	  
	  
	}
	
      }
    }
  }


  if(!Skip_HLT) {

    //print the results
    

    boost::filesystem::create_directory("../data/HLT_virtual");
    for(int iens=0; iens < (signed)Ens_Tag.size() ; iens++ ) {

      boost::filesystem::create_directory("../data/HLT_virtual/"+Ens_Tag[iens]);
      for(int xg=0; xg < (signed)xgamma_list.size(); xg++) {

	Vfloat xk_list_tailored;
	if(xgamma_list[xg]=="0.5") {xk_list_tailored = {0.6,0.65,0.7,0.75} ; }
	else if(xgamma_list[xg]=="0.7") {xk_list_tailored = {0.6} ; }
	else if(xgamma_list[xg]=="0.3") {xk_list_tailored = {0.6,0.65,0.7,0.75,0.8,0.85};}
	else if(xgamma_list[xg]=="0.1") {xk_list_tailored = {0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95} ; }
	else crash("xgamma value: "+xgamma_list[xg]+" not allowed!");
	
	boost::filesystem::create_directory("../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]);

	//print as a function of epsilon for fixed xk

	for(int xk=0; xk<(signed)xk_list_tailored.size(); xk++) {
	  Print_To_File({}, { eps_list,  H1_RE[iens][xg][xk].ave(), H1_RE[iens][xg][xk].err(), H1_IM[iens][xg][xk].ave(), H1_IM[iens][xg][xk].err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/H1_xk_"+xk_list[xk],"","");
	  Print_To_File({}, { eps_list,  H1_RE_odg[iens][xg][xk].ave(), H1_RE_odg[iens][xg][xk].err(), H1_IM_odg[iens][xg][xk].ave(), H1_IM_odg[iens][xg][xk].err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/H1_xk_"+xk_list[xk]+"_odg","","");
	  Print_To_File({}, { eps_list,  H1_RE_naive[iens][xg][xk].ave(), H1_RE_naive[iens][xg][xk].err(), H1_IM_naive[iens][xg][xk].ave(), H1_IM_naive[iens][xg][xk].err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/H1_xk_"+xk_list[xk]+"_naive","","");

	  Print_To_File({}, { eps_list,  H2_RE[iens][xg][xk].ave(), H2_RE[iens][xg][xk].err(), H2_IM[iens][xg][xk].ave(), H2_IM[iens][xg][xk].err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/H2_xk_"+xk_list[xk],"","");
	  Print_To_File({}, { eps_list,  H2_RE_odg[iens][xg][xk].ave(), H2_RE_odg[iens][xg][xk].err(), H2_IM_odg[iens][xg][xk].ave(), H2_IM_odg[iens][xg][xk].err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/H2_xk_"+xk_list[xk]+"_odg","","");
	  Print_To_File({}, { eps_list,  H2_RE_naive[iens][xg][xk].ave(), H2_RE_naive[iens][xg][xk].err(), H2_IM_naive[iens][xg][xk].ave(), H2_IM_naive[iens][xg][xk].err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/H2_xk_"+xk_list[xk]+"_naive","","");

	  Print_To_File({}, { eps_list,  FA_RE[iens][xg][xk].ave(), FA_RE[iens][xg][xk].err(), FA_IM[iens][xg][xk].ave(), FA_IM[iens][xg][xk].err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/FA_xk_"+xk_list[xk],"","");
	  Print_To_File({}, { eps_list,  FA_RE_odg[iens][xg][xk].ave(), FA_RE_odg[iens][xg][xk].err(), FA_IM_odg[iens][xg][xk].ave(), FA_IM_odg[iens][xg][xk].err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/FA_xk_"+xk_list[xk]+"_odg","","");
	  Print_To_File({}, { eps_list,  FA_RE_naive[iens][xg][xk].ave(), FA_RE_naive[iens][xg][xk].err(), FA_IM_naive[iens][xg][xk].ave(), FA_IM_naive[iens][xg][xk].err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/FA_xk_"+xk_list[xk]+"_naive","","");
	  
	  Print_To_File({}, { eps_list,  FV_RE[iens][xg][xk].ave(), FV_RE[iens][xg][xk].err(), FV_IM[iens][xg][xk].ave(), FV_IM[iens][xg][xk].err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/FV_xk_"+xk_list[xk],"","");
	  Print_To_File({}, { eps_list,  FV_RE_odg[iens][xg][xk].ave(), FV_RE_odg[iens][xg][xk].err(), FV_IM_odg[iens][xg][xk].ave(), FV_IM_odg[iens][xg][xk].err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/FV_xk_"+xk_list[xk]+"_odg","","");
	  Print_To_File({}, { eps_list,  FV_RE_naive[iens][xg][xk].ave(), FV_RE_naive[iens][xg][xk].err(), FV_IM_naive[iens][xg][xk].ave(), FV_IM_naive[iens][xg][xk].err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/FV_xk_"+xk_list[xk]+"_naive","","");

	  

	  
	}
	
	for(int eps=0;eps<(signed)eps_list.size(); eps++) {
	  //reshuffle the results
	  
	  distr_t_list H1_RE_reshuffled(UseJack);
	  distr_t_list H1_IM_reshuffled(UseJack);
	  distr_t_list H1_RE_reshuffled_odg(UseJack);
	  distr_t_list H1_IM_reshuffled_odg(UseJack);
	  distr_t_list H1_RE_reshuffled_naive(UseJack);
	  distr_t_list H1_IM_reshuffled_naive(UseJack);

	  distr_t_list H2_RE_reshuffled(UseJack);
	  distr_t_list H2_IM_reshuffled(UseJack);
	  distr_t_list H2_RE_reshuffled_odg(UseJack);
	  distr_t_list H2_IM_reshuffled_odg(UseJack);
	  distr_t_list H2_RE_reshuffled_naive(UseJack);
	  distr_t_list H2_IM_reshuffled_naive(UseJack);

	  distr_t_list FA_RE_reshuffled(UseJack);
	  distr_t_list FA_IM_reshuffled(UseJack);
	  distr_t_list FA_RE_reshuffled_odg(UseJack);
	  distr_t_list FA_IM_reshuffled_odg(UseJack);
	  distr_t_list FA_RE_reshuffled_naive(UseJack);
	  distr_t_list FA_IM_reshuffled_naive(UseJack);

	  distr_t_list FV_RE_reshuffled(UseJack);
	  distr_t_list FV_IM_reshuffled(UseJack);
	  distr_t_list FV_RE_reshuffled_odg(UseJack);
	  distr_t_list FV_IM_reshuffled_odg(UseJack);
	  distr_t_list FV_RE_reshuffled_naive(UseJack);
	  distr_t_list FV_IM_reshuffled_naive(UseJack);

	 
	  

	  for(int xk=0; xk < (signed)xk_list_tailored.size(); xk++) {

	    H1_RE_reshuffled.distr_list.push_back( H1_RE[iens][xg][xk].distr_list[eps] );
	    H1_IM_reshuffled.distr_list.push_back( H1_IM[iens][xg][xk].distr_list[eps] );
	    H1_RE_reshuffled_odg.distr_list.push_back( H1_RE_odg[iens][xg][xk].distr_list[eps] );
	    H1_IM_reshuffled_odg.distr_list.push_back( H1_IM_odg[iens][xg][xk].distr_list[eps] );
	    H1_RE_reshuffled_naive.distr_list.push_back( H1_RE_naive[iens][xg][xk].distr_list[eps] );
	    H1_IM_reshuffled_naive.distr_list.push_back( H1_IM_naive[iens][xg][xk].distr_list[eps] );

	    H2_RE_reshuffled.distr_list.push_back( H2_RE[iens][xg][xk].distr_list[eps] );
	    H2_IM_reshuffled.distr_list.push_back( H2_IM[iens][xg][xk].distr_list[eps] );
	    H2_RE_reshuffled_odg.distr_list.push_back( H2_RE_odg[iens][xg][xk].distr_list[eps] );
	    H2_IM_reshuffled_odg.distr_list.push_back( H2_IM_odg[iens][xg][xk].distr_list[eps] );
	    H2_RE_reshuffled_naive.distr_list.push_back( H2_RE_naive[iens][xg][xk].distr_list[eps] );
	    H2_IM_reshuffled_naive.distr_list.push_back( H2_IM_naive[iens][xg][xk].distr_list[eps] );

	    FA_RE_reshuffled.distr_list.push_back( FA_RE[iens][xg][xk].distr_list[eps] );
	    FA_IM_reshuffled.distr_list.push_back( FA_IM[iens][xg][xk].distr_list[eps] );
	    FA_RE_reshuffled_odg.distr_list.push_back( FA_RE_odg[iens][xg][xk].distr_list[eps] );
	    FA_IM_reshuffled_odg.distr_list.push_back( FA_IM_odg[iens][xg][xk].distr_list[eps] );
	    FA_RE_reshuffled_naive.distr_list.push_back( FA_RE_naive[iens][xg][xk].distr_list[eps] );
	    FA_IM_reshuffled_naive.distr_list.push_back( FA_IM_naive[iens][xg][xk].distr_list[eps] );

	    FV_RE_reshuffled.distr_list.push_back( FV_RE[iens][xg][xk].distr_list[eps] );
	    FV_IM_reshuffled.distr_list.push_back( FV_IM[iens][xg][xk].distr_list[eps] );
	    FV_RE_reshuffled_odg.distr_list.push_back( FV_RE_odg[iens][xg][xk].distr_list[eps] );
	    FV_IM_reshuffled_odg.distr_list.push_back( FV_IM_odg[iens][xg][xk].distr_list[eps] );
	    FV_RE_reshuffled_naive.distr_list.push_back( FV_RE_naive[iens][xg][xk].distr_list[eps] );
	    FV_IM_reshuffled_naive.distr_list.push_back( FV_IM_naive[iens][xg][xk].distr_list[eps] );

	  }


	  //H1
	  Print_To_File({}, { xk_list_tailored , H1_RE_reshuffled.ave(), H1_RE_reshuffled.err(), H1_IM_reshuffled.ave(), H1_IM_reshuffled.err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/H1_eps_"+to_string_with_precision(eps,3),"","");
	  Print_To_File({}, { xk_list_tailored , H1_RE_reshuffled_odg.ave(), H1_RE_reshuffled_odg.err(), H1_IM_reshuffled_odg.ave(), H1_IM_reshuffled_odg.err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/H1_eps_"+to_string_with_precision(eps,3)+"_odg","","");
	  Print_To_File({}, { xk_list_tailored , H1_RE_reshuffled_naive.ave(), H1_RE_reshuffled_naive.err(), H1_IM_reshuffled_naive.ave(), H1_IM_reshuffled_naive.err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/H1_eps_"+to_string_with_precision(eps,3)+"_naive","","");

	  //H2
	  Print_To_File({}, { xk_list_tailored , H2_RE_reshuffled.ave(), H2_RE_reshuffled.err(), H2_IM_reshuffled.ave(), H2_IM_reshuffled.err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/H2_eps_"+to_string_with_precision(eps,3),"","");
	  Print_To_File({}, { xk_list_tailored , H2_RE_reshuffled_odg.ave(), H2_RE_reshuffled_odg.err(), H2_IM_reshuffled_odg.ave(), H2_IM_reshuffled_odg.err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/H2_eps_"+to_string_with_precision(eps,3)+"_odg","","");
	  Print_To_File({}, { xk_list_tailored , H2_RE_reshuffled_naive.ave(), H2_RE_reshuffled_naive.err(), H2_IM_reshuffled_naive.ave(), H2_IM_reshuffled_naive.err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/H2_eps_"+to_string_with_precision(eps,3)+"_naive","","");

	  //FA
	  Print_To_File({}, { xk_list_tailored , FA_RE_reshuffled.ave(), FA_RE_reshuffled.err(), FA_IM_reshuffled.ave(), FA_IM_reshuffled.err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/FA_eps_"+to_string_with_precision(eps,3),"","");
	  Print_To_File({}, { xk_list_tailored , FA_RE_reshuffled_odg.ave(), FA_RE_reshuffled_odg.err(), FA_IM_reshuffled_odg.ave(), FA_IM_reshuffled_odg.err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/FA_eps_"+to_string_with_precision(eps,3)+"_odg","","");
	  Print_To_File({}, { xk_list_tailored , FA_RE_reshuffled_naive.ave(), FA_RE_reshuffled_naive.err(), FA_IM_reshuffled_naive.ave(), FA_IM_reshuffled_naive.err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/FA_eps_"+to_string_with_precision(eps,3)+"_naive","","");


	  //FV
	  Print_To_File({}, { xk_list_tailored , FV_RE_reshuffled.ave(), FV_RE_reshuffled.err(), FV_IM_reshuffled.ave(), FV_IM_reshuffled.err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/FV_eps_"+to_string_with_precision(eps,3),"","");
	  Print_To_File({}, { xk_list_tailored , FV_RE_reshuffled_odg.ave(), FV_RE_reshuffled_odg.err(), FV_IM_reshuffled_odg.ave(), FV_IM_reshuffled_odg.err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/FV_eps_"+to_string_with_precision(eps,3)+"_odg","","");
	  Print_To_File({}, { xk_list_tailored , FV_RE_reshuffled_naive.ave(), FV_RE_reshuffled_naive.err(), FV_IM_reshuffled_naive.ave(), FV_IM_reshuffled_naive.err() }, "../data/HLT_virtual/"+Ens_Tag[iens]+"/xg_"+xgamma_list[xg]+"/FV_eps_"+to_string_with_precision(eps,3)+"_naive","","");
	
	  
	}


	//print jackknives for all epsilon and xk
	boost::filesystem::create_directory("../data/HLT_virtual/jackknives");
	boost::filesystem::create_directory("../data/HLT_virtual/jackknives/"+Ens_Tag[iens]);
	boost::filesystem::create_directory("../data/HLT_virtual/jackknives/"+Ens_Tag[iens]+"/FA");
	boost::filesystem::create_directory("../data/HLT_virtual/jackknives/"+Ens_Tag[iens]+"/FV");
	boost::filesystem::create_directory("../data/HLT_virtual/jackknives/"+Ens_Tag[iens]+"/H1");
	boost::filesystem::create_directory("../data/HLT_virtual/jackknives/"+Ens_Tag[iens]+"/H2");
	boost::filesystem::create_directory("../data/HLT_virtual/jackknives/"+Ens_Tag[iens]+"/FA/xg_"+xgamma_list[xg]);
	boost::filesystem::create_directory("../data/HLT_virtual/jackknives/"+Ens_Tag[iens]+"/FV/xg_"+xgamma_list[xg]);
	boost::filesystem::create_directory("../data/HLT_virtual/jackknives/"+Ens_Tag[iens]+"/H1/xg_"+xgamma_list[xg]);
	boost::filesystem::create_directory("../data/HLT_virtual/jackknives/"+Ens_Tag[iens]+"/H2/xg_"+xgamma_list[xg]);
	for(int ixk=0;ixk<(signed)xk_list_tailored.size();ixk++) {
	  for(int ieps=0;ieps<(signed)eps_list.size();ieps++) {
	    //print jack FA
	    Print_To_File({},{FA_RE[iens][xg][ixk].distr_list[ieps].distr, FA_IM[iens][xg][ixk].distr_list[ieps].distr}, "../data/HLT_virtual/jackknives/"+Ens_Tag[iens]+"/FA/xg_"+xgamma_list[xg]+"/xk_"+to_string_with_precision(xk_list_tailored[ixk],2)+"_eps_"+to_string_with_precision(eps_list[ieps],3)+".jack","","");
	    //print jack FV
	    Print_To_File({},{FV_RE[iens][xg][ixk].distr_list[ieps].distr, FV_IM[iens][xg][ixk].distr_list[ieps].distr}, "../data/HLT_virtual/jackknives/"+Ens_Tag[iens]+"/FV/xg_"+xgamma_list[xg]+"/xk_"+to_string_with_precision(xk_list_tailored[ixk],2)+"_eps_"+to_string_with_precision(eps_list[ieps],3)+".jack","","");
	    //print jack H1
	    Print_To_File({},{H1_RE[iens][xg][ixk].distr_list[ieps].distr, H1_IM[iens][xg][ixk].distr_list[ieps].distr}, "../data/HLT_virtual/jackknives/"+Ens_Tag[iens]+"/H1/xg_"+xgamma_list[xg]+"/xk_"+to_string_with_precision(xk_list_tailored[ixk],2)+"_eps_"+to_string_with_precision(eps_list[ieps],3)+".jack","","");
	    //print jack H2
	    Print_To_File({},{H2_RE[iens][xg][ixk].distr_list[ieps].distr, H2_IM[iens][xg][ixk].distr_list[ieps].distr}, "../data/HLT_virtual/jackknives/"+Ens_Tag[iens]+"/H2/xg_"+xgamma_list[xg]+"/xk_"+to_string_with_precision(xk_list_tailored[ixk],2)+"_eps_"+to_string_with_precision(eps_list[ieps],3)+".jack","","");
	  }
	}  
	
      }
    }
  }


  


 
  //Gounaris Sakurai model

  int npts_spline= 300;
  int Luscher_num_zeroes= 4;
  int Nresonances=2;
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
    cout<<"Computing n(Lusch): "<<L_zero<<endl;
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

    for(int istep=1;istep<=npts_spline-1;istep++) { double pt= sx+dz*istep;
      phi_data[L_zero].push_back( phi(sqrt(pt)));}

    phi_data[L_zero].push_back(M_PI/2.0);
    double sx_der_loc =  phi_der_for_back(sqrt(sx)+1e-14, 1);
    double dx_der_loc =  phi_der_for_back(sqrt(dx)-1e-14, -1);
    sx_der.push_back(sx_der_loc);
    dx_der.push_back(dx_der_loc);

    phi_der_data[L_zero].push_back(sx_der_loc);
    for(int istep=1;istep<=npts_spline-1;istep++) { double pt= sx+dz*istep;
      phi_der_data[L_zero].push_back( phi_der(sqrt(pt)));}
    phi_der_data[L_zero].push_back(dx_der_loc);

  }


  LL_functions LL(phi_data,phi_der_data,sx_der, dx_der, sx_int, Dz, Nresonances, Luscher_zeroes);

  //###########################################END INTERPOLATION PHI FUNCTION AND DERIVATIVES################################
  cout<<"####Spline for phi(x) and phi'(x) successfully generated!"<<endl;
  

  double vol1= 3.8160*fm_to_inv_Gev;
  double vol2= 5.09*fm_to_inv_Gev;
  double vol3= 7.63*fm_to_inv_Gev;
  
  Vfloat En1;
  LL.Find_pipi_energy_lev(vol1, 0.775, 5.0, 0.135, 0.0, En1);
  Vfloat En2;
  LL.Find_pipi_energy_lev(vol2, 0.775, 5.0, 0.135, 0.0, En2);
  Vfloat En3;
  LL.Find_pipi_energy_lev(vol3, 0.775, 5.0, 0.135, 0.0, En3);


 


  vector<double> eps_list_GS { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.125, 0.150, 0.175, 0.20, 0.25, 0.3, 0.35}; //GeV


  // Kl4-enhanced Gounaris-Sakurai model


  for(int ixg=0;ixg<(signed)xgamma_list.size() ; ixg++) {


    bool no_R=false;

    double xg= stod(xgamma_list[ixg]);
    cout<<"Computing model for xg: "<<xg<<endl;
    double k= MK_FLAG*xg/2;
    double t=1;
    double m= MK_FLAG;
    double th = 2*sqrt( pow(0.135,2) + pow(k/2,2));
    double Emax= 4*MK_FLAG; //sqrt( pow(1.0,2) + pow(k,2));
    double ek = MK_FLAG*sqrt( 1 + pow(xg/2,2) );
    
    
    auto C_rho_12_V = [&xg, &k, &LL, &t, &m](double E) {
      double s= E*E - k*k;
      double R=1; //LL.F_pi_GS_mod(2*0.135,0.775,5.0,0.135,0.0);
      double LEC= -0.4587*5.0/R;
      double Fem= pow(LL.F_pi_GS_mod(sqrt(s), 0.775, 5.0,0.135,0.0),2);
      double std_f= pow(s- pow(2*0.135,2) ,3.0/2)/(12.0*pow(2*M_PI,2));
      double res= k*std_f*Fem/(m*m*sqrt(s));
      
      
      return LEC*res*exp(-E*t);
    };
    
    auto C_rho_11_A = [&xg, &k, &LL, &t, &m](double E) {
      double s= E*E - k*k;
      double q2= s/pow(2*0.135,2) -1;
      double R=1.0; //1.0/LL.F_pi_GS_mod(2*0.135,0.775,5.0,0.135,0.0);
      double Gs= 5.0/R; //5.0;
      double Fem= pow(LL.F_pi_GS_mod(sqrt(s), 0.775, 5.0,0.135,0.0),2);
      double std_f= pow(s- pow(2*0.135,2) ,3.0/2)/(12.0*pow(2*M_PI,2));
      double res= -std_f*Gs*Fem/(m*sqrt(s));
      return -res*exp(-E*t);
    };

    
    auto C_rho_33_A = [&xg, &k, &LL, &t, &m, &no_R](double E) {
      double s= E*E - k*k;
      double R=1.0; //LL.F_pi_GS_mod(2*0.135,0.775,5.0,0.135,0.0);
      double X= m*k; 
      double fovG= -0.055;
      double hatfovG = fovG*m*m/(X);
      double q2= s/pow(2*0.135,2) -1;
      double Gs= 5.0/R;
      double RovG =(no_R==true)?0.0:-1.4;
      //RovG= (m*m/(s -2*m*sqrt( s + k*k)));
      double Fem= pow(LL.F_pi_GS_mod(sqrt(s), 0.775, 5.0,0.135,0.0),2);
      double std_f= pow(s- pow(2*0.135,2) ,3.0/2)/(12.0*pow(2*M_PI,2));
      double std_ft= pow(s,3.0/2)*(1.0 -pow(2*0.135,2)/s)/(12.0*pow(2*M_PI,2));
      double res_1= -Gs*std_f*Fem*E*E/(m*pow(s,3.0/2));
      double res_2= -Gs*std_f*k*k*E*(RovG*Fem)/(m*m*pow(s,3.0/2));
      double res_3= -Gs*std_ft*k*k*E*(hatfovG*Fem)/(m*m*pow(s,3.0/2));
      return -(res_1 + res_2 + 0*res_3)*exp(-E*t);
    };
    
    auto C_rho_03_A  = [&xg, &k, &LL, &t, &m, &no_R](double E) {
      
      double s= E*E - k*k;
      double X=m*k;
      double fovG= -0.055;
      double hatfovG = fovG*m*m/(X);
      double RovG =(no_R==true)?0.0:-1.4;
      double q2= s/pow(2*0.135,2) -1;
      double Gs= 5.0+ 0.0*q2;
      //RovG= (m*m/(s -2*m*sqrt( s + k*k)));

      double Fem= pow(LL.F_pi_GS_mod(sqrt(s), 0.775, 5.0,0.135,0.0),2);
      double std_f= pow(s- pow(2*0.135,2) ,3.0/2)/(12.0*pow(2*M_PI,2));
      double std_ft= pow(s,3.0/2)*(1.0 -pow(2*0.135,2)/s)/(12.0*pow(2*M_PI,2));
      double res_1= -Gs*std_f*Fem*k*E/(m*pow(s,3.0/2));
      double res_2= -Gs*std_f*k*k*k*(RovG*Fem)/(m*m*pow(s,3.0/2));
      double res_3= -Gs*std_ft*k*k*k*(hatfovG*Fem)/(m*m*pow(s,3.0/2));
      return -(res_1 + res_2 + 0*res_3)*exp(-E*t);
    };
    
    
    auto C_rho_30_A  = [&xg, &k, &LL, &t, &m, &no_R](double E) {
      
      double s= E*E - k*k;
          
      double X= m*k;
      double fovG= -0.055;
      double hatfovG = fovG*m*m/(X);
      double RovG =(no_R==true)?0.0:-1.4;
      //RovG= (m*m/(s -2*m*sqrt( s + k*k)));
    
      double Fem= pow(LL.F_pi_GS_mod(sqrt(s), 0.775, 5.0,0.135,0.0),2);
      double q2= s/pow(2*0.135,2) -1;
      double Gs= 5.0*(1.0+ 0.0*q2);

      double std_f= pow(s- pow(2*0.135,2) ,3.0/2)/(12.0*pow(2*M_PI,2));
      double std_ft= pow(s,3.0/2)*(1.0 -pow(2*0.135,2)/s)/(12.0*pow(2*M_PI,2));
      double res_1= -Gs*std_f*Fem*k*E/(m*pow(s,3.0/2));
      double res_2= -Gs*std_f*k*E*(RovG*Fem)*(E/m -1.0)/(m*pow(s,3.0/2));
      double res_3= -Gs*std_ft*k*E*(hatfovG*Fem)/(m*pow(s,3.0/2));
      return -(res_1 + res_2+ 0*res_3)*exp(-E*t);
    };

    double dt=0.01*fm_to_inv_Gev; //fm 
    int N=400;

    Vfloat times;
    Vfloat rat_33_11;
    Vfloat rat_30_03;
    Vfloat rat_12_11;

    Vfloat corr_30_A, corr_03_A, corr_33_A, corr_11_A, corr_12_V, corr_30_A_kaon, corr_33_A_kaon, corr_03_A_kaon;
  
    for(int n=1; n<N;n++) {

      //cout<<"Computing time: "<<n*dt/fm_to_inv_Gev<<" fm"<<endl;
    
      t= n*dt;
      double tolerance =1e-5;

      times.push_back( t/fm_to_inv_Gev);

      //##### V12 ######
      double valV,errV;
      gsl_function_pp<decltype(C_rho_12_V)> F_V(C_rho_12_V);
      gsl_integration_workspace * wV = gsl_integration_workspace_alloc (30000);
      gsl_function *GV = static_cast<gsl_function*>(&F_V);
      gsl_integration_qags(GV, th, Emax,  0.0, tolerance, 30000, wV, &valV, &errV);
      gsl_integration_workspace_free (wV);
      if(fabs(errV)/valV > 1e-4) cout<<"Warning precision reached: "<<errV/fabs(valV)<<endl;
      //##### A11 #####
      double valA11,errA11;
      gsl_function_pp<decltype(C_rho_11_A)> F_A11(C_rho_11_A);
      gsl_integration_workspace * wA11 = gsl_integration_workspace_alloc (30000);
      gsl_function *GA11 = static_cast<gsl_function*>(&F_A11);
      gsl_integration_qags(GA11, th, Emax, 0.0, tolerance, 30000, wA11, &valA11, &errA11);
      gsl_integration_workspace_free (wA11);
      if(fabs(errA11)/valA11 > 1e-4) cout<<"Warning precision reached: "<<errA11/fabs(valA11)<<endl;
      //##### A33 #####
      double valA33,errA33;
      gsl_function_pp<decltype(C_rho_33_A)> F_A33(C_rho_33_A);
      gsl_integration_workspace * wA33 = gsl_integration_workspace_alloc (30000);
      gsl_function *GA33 = static_cast<gsl_function*>(&F_A33);
      gsl_integration_qags(GA33, th, Emax,0.0, tolerance, 30000, wA33, &valA33, &errA33);
      gsl_integration_workspace_free (wA33);
      if(fabs(errA33)/valA33 > 1e-4) cout<<"Warning precision reached: "<<errA33/fabs(valA33)<<endl;
      //##### A30 #####
      double valA30,errA30;
      gsl_function_pp<decltype(C_rho_30_A)> F_A30(C_rho_30_A);
      gsl_integration_workspace * wA30 = gsl_integration_workspace_alloc (30000);
      gsl_function *GA30 = static_cast<gsl_function*>(&F_A30);
      gsl_integration_qags(GA30, th, Emax, 0.0, tolerance, 30000, wA30, &valA30, &errA30);
      gsl_integration_workspace_free (wA30);
      if(fabs(errA30)/valA30 > 1e-4) cout<<"Warning precision reached: "<<errA30/fabs(valA30)<<endl;
      //##### A03 #####
      double valA03,errA03;
      gsl_function_pp<decltype(C_rho_03_A)> F_A03(C_rho_03_A);
      gsl_integration_workspace * wA03 = gsl_integration_workspace_alloc (30000);
      gsl_function *GA03 = static_cast<gsl_function*>(&F_A03);
      gsl_integration_qags(GA03, th, Emax,  0.0, tolerance, 30000, wA03, &valA03, &errA03); //qagiu
      gsl_integration_workspace_free (wA03);
      if(fabs(errA03)/valA03 > 1e-4) cout<<"Warning precision reached: "<<errA03/fabs(valA03)<<endl;

  

      double fk=0.156;
      double Fem= 1.3;


      double kaon_30 = Fem*exp( -(m +ek)*t)*k*ek*(fk/(2*ek));
      double kaon_33 = Fem*exp( -(m +ek)*t)*k*k*(fk/(2*ek));
      double kaon_03 = Fem*exp( -(m +ek)*t)*(ek-m)*k*(fk/(2*ek));

      rat_30_03.push_back( (valA30+kaon_30)/(valA03+kaon_03));
      rat_33_11.push_back( (valA33+kaon_33)/valA11);
      rat_12_11.push_back( valV/valA11);
    

      corr_30_A.push_back(valA30);
      corr_03_A.push_back(valA03);
      corr_11_A.push_back(valA11);
      corr_33_A.push_back(valA33);
      corr_12_V.push_back(valV);
      corr_30_A_kaon.push_back(kaon_30);
      corr_33_A_kaon.push_back(kaon_33);
      corr_03_A_kaon.push_back(kaon_03);
    
    
    
    }


    //print model prediction to file

    boost::filesystem::create_directory("../data/HLT_virtual/GS_Kl4");

    Print_To_File({ }, {times, rat_30_03, rat_33_11, rat_12_11} ,  "../data/HLT_virtual/GS_Kl4/xg_"+to_string_with_precision(xg,3),"","");

    Print_To_File({ }, {times, corr_12_V, corr_11_A, corr_33_A, corr_30_A, corr_03_A, corr_30_A_kaon, corr_33_A_kaon} ,  "../data/HLT_virtual/GS_Kl4/corr_xg_"+to_string_with_precision(xg,3),"","");


    corr_30_A = summ_master(corr_30_A, corr_30_A_kaon);
    corr_33_A = summ_master(corr_33_A, corr_33_A_kaon);
    corr_03_A = summ_master(corr_03_A, corr_03_A_kaon);

    for(int ixk=0;ixk<(signed)xk_list_double.size();ixk++) {

      double MP=MK_FLAG;
      double kz= k;
      double xk = xk_list_double[ixk];
      double k2= pow(xk*MP,2);
      double Eg_v = sqrt( k2 + k*k);
      
    
      Eigen::MatrixXd A(3,3); //coefficient matrix

      
      //   [   H1    ]     =    [           ]^-1       [       H30 -H03*(MP-Eg)/(2MP-Eg)          ]
      //   [   H2    ]     =    [     A     ]     *    [              H33                         ]
      //   [   FA    ]     =    [           ]          [              H11                         ]
       

      //                   [  -Eg*kz/(2MP-Eg)            -kz*(MP-Eg)/(2MP -Eg)           kz*MP/(2MP-Eg)                 ]
      //         A  =      [  -Eg*Eg/MP                   Eg*kz^2/(2MP*Eg -k^2)         -Eg*(MP-Eg)/MP                  ]
      //                   [  -k^2/MP                             0                     -(MP*Eg -k^2)/MP                ]


      int sign=1;

      //A(0,0) = -sign*Eg_v*kz/MP ; A(0,1) = -sign*Eg_v*(MP-Eg_v)/(2*MP*Eg_v -k2); A(0,2) = sign*Eg_v*kz/MP; 
    
      A(0,0) =  -sign*Eg_v*kz/(2*MP-Eg_v);     A(0,1) = -sign*kz*(MP-Eg_v)/(2*MP-Eg_v);         A(0,2) = sign*kz*MP/(2*MP-Eg_v);
    
      A(1,0) = -Eg_v*Eg_v/MP;                  A(1,1) = Eg_v*kz*kz/(2*MP*Eg_v - k2);            A(1,2) = -Eg_v*(MP-Eg_v)/MP;
    
      A(2,0) = -k2/MP;                         A(2,1) = 0.0;                                    A(2,2) = -(MP*Eg_v - k2)/MP;

      Eigen::MatrixXd B = A.inverse();

      //cout<<"B-matrix for xg: "<<xg<<", xk: "<<xk<<endl;
      //cout<<B<<endl<<endl;

      //from the correlator to the form factors
    
      Vfloat H1_corr, H2_corr, FV_corr, FA_corr;
      for(int t=0;t<(signed)times.size();t++) {
      
      
	Eigen::VectorXd X(3);
      
	X(0) = corr_30_A[t] - corr_03_A[t]*(MP-Eg_v)/(2*MP-Eg_v);
	X(1) = corr_33_A[t];
	X(2) = corr_11_A[t];
      
	Eigen::VectorXd Y = B*X;
      
	H1_corr.push_back( Y(0));
	H2_corr.push_back( Y(1));
	FA_corr.push_back( Y(2));
	FV_corr.push_back( -corr_12_V[t]/k);
      
      }
    
    
      Print_To_File({ }, {times, H1_corr,FA_corr, H2_corr, FV_corr} ,  "../data/HLT_virtual/GS_Kl4/F_corr_xg_"+to_string_with_precision(xg,3)+"_xk_"+to_string_with_precision(xk,2),"","");

      //compute the smeared amplitude

      for(int iff=0;iff<4;iff++) {


	//KK pole contribution

	double fk=0.156;
	double Fem= 1.3;
	Eigen::VectorXd X_F_Kaon(3);	
	X_F_Kaon(2) = 0.0;
	X_F_Kaon(1) = Fem*k*k*(fk/(2*ek));
	X_F_Kaon(0) = Fem*k*ek*(fk/(2*ek)) -  Fem*(ek-m)*k*(fk/(2*ek))*(MP-Eg_v)/(2*MP-Eg_v);


	double amp= 0.0;
	Eigen::VectorXd Y_F_Kaon= B*X_F_Kaon;
	if(iff<3) amp = Y_F_Kaon(iff);
      
	Vfloat RE_F_mod, IM_F_mod, RE_F_mod_noK, RE_F_mod_noR, IM_F_mod_noR;
      
	for(int ie=0; ie<(signed)eps_list_GS.size(); ie++) {

	  double eps= eps_list_GS[ie];
	  t=0;

	  bool is_eps_0 = false;
	  if(fabs(eps) < 1e-10) is_eps_0=true;

	  auto IM_rho_F = [&](double E) {

	    Eigen::VectorXd X_F(3);
	
	    X_F(2)=  C_rho_11_A(E);
	    X_F(1)=  C_rho_33_A(E);
	    X_F(0) = C_rho_30_A(E) - C_rho_03_A(E)*(MP-Eg_v)/(2*MP-Eg_v);	
	
	    double X_V= -C_rho_12_V(E)/k;

	    Eigen::VectorXd Y_F = B*X_F;

	    double X= (iff==3)?X_V:Y_F(iff);

	    if(is_eps_0) return X*M_PI;
	
	    return X*eps/( (E-Eg_v)*(E-Eg_v) + eps*eps);
	  };

	  auto RE_rho_F = [&](double E) {

	    Eigen::VectorXd X_F(3);
	
	    X_F(2)=  C_rho_11_A(E);
	    X_F(1)=  C_rho_33_A(E);
	    X_F(0) = C_rho_30_A(E) - C_rho_03_A(E)*(MP-Eg_v)/(2*MP-Eg_v);

	    Eigen::VectorXd Y_F = B*X_F;

	    double X_V= -C_rho_12_V(E)/k;
	
	    double X= (iff==3)?X_V:Y_F(iff);

	    if(is_eps_0) { return X*( Eg_v/E); }

	    if(is_eps_0) crash("you shold not be here");

	    return X*( (E-Eg_v)/( (E-Eg_v)*(E-Eg_v) + eps*eps)   - (1.0/(E)) );
	  };


	  //integrate

	  

	  if(fabs(eps) > 1e-10) {

	    double val,err;
	    double tolerance = 1e-5;
	    gsl_function_pp<decltype(RE_rho_F)> F_RE_corr(RE_rho_F);
	    gsl_integration_workspace * w_RE = gsl_integration_workspace_alloc (30000);
	    gsl_function *G_RE = static_cast<gsl_function*>(&F_RE_corr);
	    gsl_integration_qags(G_RE, th, Emax, 0.0, tolerance, 30000, w_RE, &val, &err);
	    gsl_integration_workspace_free (w_RE);
	    if( err/fabs(val) > 5*tolerance) crash("In RE-mod inf-vol, for ixk: "+to_string(xk)+" eps: "+to_string(eps)+", gls integration not able to achieve target precision. Precision reached: "+to_string_with_precision(err/fabs(val),7));
	    double val_kaon_RE=amp*( (m+ek-Eg_v)/( (m+ek-Eg_v)*(m+ek-Eg_v) + eps*eps)   - (1.0/(m+ek)) );
	    RE_F_mod.push_back(val+ val_kaon_RE);
	    RE_F_mod_noK.push_back(val);	    
	    gsl_function_pp<decltype(IM_rho_F)> F_IM_corr(IM_rho_F);
	    gsl_integration_workspace * w_IM = gsl_integration_workspace_alloc (30000);
	    gsl_function *G_IM = static_cast<gsl_function*>(&F_IM_corr);
	    gsl_integration_qags(G_IM, th, Emax, 0.0, tolerance, 30000, w_IM, &val, &err);
	    gsl_integration_workspace_free (w_IM);
	    if( err/fabs(val) > 5*tolerance) crash("In IM-mod inf-vol, for ixk: "+to_string(xk)+" eps: "+to_string(eps)+", gls integration not able to achieve target precision. Precision reached: "+to_string_with_precision(err/fabs(val),7));
	    double val_kaon_IM= amp*eps/( (m+ek-Eg_v)*(m+ek-Eg_v) + eps*eps);
	    IM_F_mod.push_back(val+val_kaon_IM);

	    no_R=true;
	    gsl_function_pp<decltype(RE_rho_F)> F_noR_RE_corr(RE_rho_F);
	    gsl_integration_workspace * w_RE_noR = gsl_integration_workspace_alloc (30000);
	    gsl_function *G_RE_noR = static_cast<gsl_function*>(&F_noR_RE_corr);
	    gsl_integration_qags(G_RE_noR, th, Emax, 0.0, tolerance, 30000, w_RE_noR, &val, &err);
	    gsl_integration_workspace_free (w_RE_noR);
	    if( err/fabs(val) > 5*tolerance) crash("In RE-mod inf-vol, for ixk: "+to_string(xk)+" eps: "+to_string(eps)+", gls integration not able to achieve target precision. Precision reached: "+to_string_with_precision(err/fabs(val),7));
	    RE_F_mod_noR.push_back(val+ val_kaon_RE);
	    gsl_function_pp<decltype(IM_rho_F)> F_noR_IM_corr(IM_rho_F);
	    gsl_integration_workspace * w_IM_noR = gsl_integration_workspace_alloc (30000);
	    gsl_function *G_IM_noR = static_cast<gsl_function*>(&F_noR_IM_corr);
	    gsl_integration_qags(G_IM_noR, th, Emax, 0.0, tolerance, 30000, w_IM_noR, &val, &err);
	    gsl_integration_workspace_free (w_IM_noR);
	    if( err/fabs(val) > 5*tolerance) crash("In IM-mod inf-vol, for ixk: "+to_string(xk)+" eps: "+to_string(eps)+", gls integration not able to achieve target precision. Precision reached: "+to_string_with_precision(err/fabs(val),7));
	    IM_F_mod_noR.push_back(val+ val_kaon_IM);

	    no_R=false;

	  }

	  else {

	    double val,err;
	    double tolerance = 1e-5;
	    is_eps_0=true;
	    gsl_function_pp<decltype(RE_rho_F)> F_RE_corr(RE_rho_F);
	    gsl_integration_workspace * w_RE = gsl_integration_workspace_alloc (30000);
	    gsl_function *G_RE = static_cast<gsl_function*>(&F_RE_corr);
	    gsl_integration_qawc(G_RE, th+1e-7, Emax, Eg_v, 0.0, tolerance, 30000, w_RE, &val, &err);
	    gsl_integration_workspace_free (w_RE);
	    if( err/fabs(val) > 5*tolerance) crash("In RE-mod inf-vol, for ixk: "+to_string(xk)+" eps: "+to_string(eps)+", gls integration not able to achieve target precision. Precision reached: "+to_string_with_precision(err/fabs(val),7));
	    double val_kaon_RE = amp*( 1.0/(m+ek-Eg_v)   - (1.0/(m+ek)) );
	    RE_F_mod.push_back(val+val_kaon_RE);
	    RE_F_mod_noK.push_back(val);
	
	    IM_F_mod.push_back( IM_rho_F(Eg_v));

	    no_R=true;
	    gsl_function_pp<decltype(RE_rho_F)> F_noR_RE_corr(RE_rho_F);
	    gsl_integration_workspace * w_RE_noR = gsl_integration_workspace_alloc (30000);
	    gsl_function *G_RE_noR = static_cast<gsl_function*>(&F_noR_RE_corr);
	    gsl_integration_qawc(G_RE_noR, th+1e-7, Emax, Eg_v, 0.0, tolerance, 30000, w_RE_noR, &val, &err);
	    gsl_integration_workspace_free (w_RE_noR);
	    if( err/fabs(val) > 5*tolerance) crash("In RE-mod inf-vol, for ixk: "+to_string(xk)+" eps: "+to_string(eps)+", gls integration not able to achieve target precision. Precision reached: "+to_string_with_precision(err/fabs(val),7));
	    RE_F_mod_noR.push_back(val+val_kaon_RE);
	    IM_F_mod_noR.push_back(IM_rho_F(Eg_v));
	    no_R=false;

	  }
      
	}

	//print
	string tag_FF="";
	if(iff==0) tag_FF="H1";
	else if(iff==1) tag_FF="H2";
	else if(iff==2) tag_FF="FA";
	else if(iff==3) tag_FF="FV";
	else crash("tag not found");
      
	Print_To_File({ }, {eps_list_GS, RE_F_mod,IM_F_mod, RE_F_mod_noK, RE_F_mod_noR, IM_F_mod_noR} ,  "../data/HLT_virtual/GS_Kl4/"+tag_FF+"_mod_xg_"+to_string_with_precision(xg,3)+"_xk_"+to_string_with_precision(xk,2),"","");
    
      }
    }

  }
  
  
  
  cout<<"Computation done! exiting!"<<endl;
  exit(-1);
  

  

  




};

