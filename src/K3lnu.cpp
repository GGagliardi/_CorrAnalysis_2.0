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
