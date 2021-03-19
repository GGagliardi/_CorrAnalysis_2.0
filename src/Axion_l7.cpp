#include "../include/Axion_l7.h"


using namespace std;
namespace plt = matplotlibcpp;

const double MPiPhys=0.135;
const bool Use_JB_distribution= false;
const bool UseJack=1;
const int Njacks=15;
const int Nboots=200;


void Neutral_pi_TM() {

  data_t pi_conn, pi_disc, pi_bubble;

   
  pi_conn.Read("../datasets_pi_0_TM", "conn_contr", "P5P5");
  pi_disc.Read("../datasets_pi_0_TM", "disc_contr", "P5");
  pi_bubble.Read("../datasets_pi_0_TM", "bubble", "P5");

  int Nens= pi_conn.size;
  cout<<"N_ens: "<<Nens<<endl;

 

  for(int i=0; i<Nens;i++) {
    distr_t_list pi_conn_distr(UseJack), pi_disc_distr(UseJack);
    CorrAnalysis Corr(UseJack, Njacks, Nboots);
    Corr.Perform_Nt_t_average=1;
    Corr.Reflection_sign=1;
    if(pi_conn.Tag[i].substr(0,1)=="A") Corr.Tmin=12;
    else if(pi_conn.Tag[i].substr(0,1)=="B") Corr.Tmin=16;
    else if(pi_conn.Tag[i].substr(0,1)=="D") Corr.Tmin=21;
    else if(pi_conn.Tag[i].substr(0,1)=="a") Corr.Tmin=12;
    else Corr.Tmin=21;
    Corr.Tmax = pi_conn.nrows[i]/2-6;
    Corr.Nt = pi_conn.nrows[i];
    //LatticeInfo L_info("LOCAL");

    //L_info.LatInfo(pi_conn.Tag[i].substr(0,1));
    double L, NNT, m;
    Read_pars_from_ensemble_tag(pi_conn.Tag[i], m, L, NNT);
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/pi_0_TM");
    pi_conn_distr=Corr.corr_t(pi_conn.col(0)[i], "../data/pi_0_TM/pi_conn."+pi_conn.Tag[i]);
    
    distr_t_list Pi_plus_eff_mass = Corr.effective_mass_t(pi_conn_distr, "../data/pi_0_TM/pi_plus_eff_mass."+pi_conn.Tag[i]);

    auto F= [&](const Vfloat& par) { if((signed)par.size() != 2) crash("Lambda function expects par[2], but par["+to_string((signed)par.size())+"] provided"); return par[0] -par[1]*par[1];};
    Vfloat bubble_average(pi_bubble.col(0)[i][0].size(),0);
    for(int t=0;t<Corr.Nt;t++) {
      if(bubble_average.size() != pi_bubble.col(0)[i][t].size()) crash("Number of confs not costants over time in Neutral_pi_TM");
      for(unsigned int iconf=0; iconf < bubble_average.size();iconf++) {
	bubble_average[iconf] += pi_bubble.col(0)[i][t][iconf]/Corr.Nt;
      }
    }
    for(int t=0;t<Corr.Nt;t++) {
	   if(UseJack==1) {
	     Jackknife J(100,Njacks);
	     pi_disc_distr.distr_list.push_back(J.DoJack(F, 2, pi_disc.col(0)[i][t], bubble_average));
	   }
	   else {
	     Bootstrap B(Nboots, 432423, Corr.Nt);
	     pi_disc_distr.distr_list.push_back(B.DoBoot(F, 2, pi_disc.col(0)[i][t], bubble_average));
	   }
    }
    Print_To_File({}, {pi_disc_distr.ave(), pi_disc_distr.err()},"../data/pi_0_TM/pi_disc."+pi_disc.Tag[i]+".t"  ,"OUT", "");

  
    distr_t_list Pi_0_eff_mass= Corr.effective_mass_t(pi_conn_distr+2.0*pi_disc_distr, "../data/pi_0_TM/pi_0_eff_mass."+pi_disc.Tag[i]);

    for(int t=0; t < Corr.Nt;t++) printV(Pi_0_eff_mass.distr_list[t].distr, "t: "+to_string(t), 1);

 
    Pfloat res_charged_pion= Corr.Fit_(Pi_plus_eff_mass);
    Pfloat res_neutral_pion= Corr.Fit_(Pi_0_eff_mass);
    cout<<"M_pi+: "<<res_charged_pion.first<<"("<<res_charged_pion.second<<")"<<endl;
    cout<<"M_pi0: "<<res_neutral_pion.first<<"("<<res_neutral_pion.second<<")"<<endl;
  }
  return;
}

void Eta_TM() {


 data_t eta_conn, eta_disc, eta_bubble;

   
  eta_conn.Read("../datasets_eta_TM", "conn_contr", "P5P5");
  eta_disc.Read("../datasets_eta_TM", "disc_contr", "P5");
  eta_bubble.Read("../datasets_eta_TM", "bubble", "P5");

  int Nens= eta_conn.size;
  cout<<"N_ens: "<<Nens<<endl;

 

  for(int i=0; i<Nens;i++) {
    distr_t_list eta_conn_distr(UseJack), eta_disc_distr(UseJack);
    CorrAnalysis Corr(UseJack, Njacks, Nboots);
    Corr.Perform_Nt_t_average=1;
    Corr.Reflection_sign=1;
    if(eta_conn.Tag[i].substr(0,1)=="A") Corr.Tmin=12;
    else if(eta_conn.Tag[i].substr(0,1)=="B") Corr.Tmin=16;
    else if(eta_conn.Tag[i].substr(0,1)=="D") Corr.Tmin=21;
    else if(eta_conn.Tag[i].substr(0,2)=="a") Corr.Tmin=12;
    else Corr.Tmin=21;
    Corr.Tmax = eta_conn.nrows[i]/2-6;
    Corr.Nt = eta_conn.nrows[i];
    //LatticeInfo L_info("LOCAL");

    //L_info.LatInfo(eta_conn.Tag[i].substr(0,1));
    double L, NNT, m;
    Read_pars_from_ensemble_tag(eta_conn.Tag[i], m, L, NNT);
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/eta_TM");
    eta_conn_distr=Corr.corr_t(eta_conn.col(0)[i], "../data/eta_TM/eta_conn."+eta_conn.Tag[i]);
    
    distr_t_list Eta_conn_eff_mass = Corr.effective_mass_t(eta_conn_distr, "../data/eta_TM/eta_conn_eff_mass."+eta_conn.Tag[i]);

    auto F= [&](const Vfloat& par) { if((signed)par.size() != 2) crash("Lambda function expects par[2], but par["+to_string((signed)par.size())+"] provided"); return par[0] -par[1]*par[1];};
    Vfloat bubble_average(eta_bubble.col(0)[i][0].size(),0);
    for(int t=0;t<Corr.Nt;t++) {
      if(bubble_average.size() != eta_bubble.col(0)[i][t].size()) crash("Number of confs not costants over time in eta_TM");
      for(unsigned int iconf=0; iconf < bubble_average.size();iconf++) {
	bubble_average[iconf] += eta_bubble.col(0)[i][t][iconf]/Corr.Nt;
      }
    }
    for(int t=0;t<Corr.Nt;t++) {
	   if(UseJack==1) {
	     Jackknife J(100,Njacks);
	     eta_disc_distr.distr_list.push_back(J.DoJack(F, 2, eta_disc.col(0)[i][t], bubble_average));
	   }
	   else {
	     Bootstrap B(Nboots, 432423, Corr.Nt);
	     eta_disc_distr.distr_list.push_back(B.DoBoot(F, 2, eta_disc.col(0)[i][t], bubble_average));
	   }
    }
    Print_To_File({}, {eta_disc_distr.ave(), eta_disc_distr.err()},"../data/eta_TM/eta_disc."+eta_disc.Tag[i]+".t"  ,"OUT", "");

    
 
    distr_t_list Eta_eff_mass= Corr.effective_mass_t(eta_conn_distr+2.0*eta_disc_distr, "../data/eta_TM/eta_eff_mass."+eta_disc.Tag[i]);
    Pfloat res_eta_conn= Corr.Fit_(Eta_conn_eff_mass);
    Pfloat res_eta= Corr.Fit_(Eta_eff_mass);
    cout<<"M_eta_conn: "<<res_eta_conn.first<<"("<<res_eta_conn.second<<")"<<endl;
    cout<<"M_eta: "<<res_eta.first<<"("<<res_eta.second<<")"<<endl;
  }
  return;
}




void Axion_l7_analysis() {
  
  data_t pi_data_r0,pi_data_r1, IE1_data, IE2_data, IH1_data, IH2_data, disc_pi0;
  data_t IB1_data, IB2_data, IS1_data;
  pi_data_r0.Read("../datasets_axion", "mes_contr_00_r0", "P5P5");
  //  pi_data_r1.Read("../datasets_axion", "mes_contr_00_r1", "P5P5");
  // disc_pi0.Read("../datasets_axion", "mes_contr_disc_pi0", "S0P5");
  IE1_data.Read("../datasets_axion", "mes_contr_IE1", "P5P5");
  //IE2_data.Read("../datasets_axion", "mes_contr_IE2", "P5P5");
  IH1_data.Read("../datasets_axion", "mes_contr_IH1", "S0P5");
  //IH2_data.Read("../datasets_axion", "mes_contr_IH2", "S0P5");
  IB1_data.Read("../datasets_axion", "mes_contr_IBub1", "S0P5");
  // IB2_data.Read("../datasets_axion", "mes_contr_IBub2", "S0P5");
  IS1_data.Read("../datasets_axion", "mes_contr_IS1", "P5P5");
  //IS2_data.Read("../datasets_axion", "mes_contr_IS2", "P5P5");
  //IS3_data.Read("../datasets_axion", "mes_contr_IS3", "P5P5");
  //IS4_data.Read("../datasets_axion", "mes_contr_IS4", "P5P5");
 

  
    
    
  
  int Nens = pi_data_r0.size;  //== dm_exch_data.size .....
  cout<<"N_ens: "<<Nens<<endl;

  distr_t_list Pi_r0_corr(UseJack), Pi_r1_corr(UseJack), Exch_r0_corr(UseJack), Exch_r1_corr(UseJack), Hand_r0_corr(UseJack), Hand_r1_corr(UseJack);
  distr_t_list Bubble_r0_R_corr(UseJack), Bubble_r1_R_corr(UseJack), Bubble_r0_I_corr(UseJack), Bubble_r1_I_corr(UseJack);
  distr_t_list IS1_corr(UseJack), IS2_corr(UseJack), IS3_corr(UseJack), IS4_corr(UseJack);
  
  distr_t_list Pi_plus_corr(UseJack), Pi_plus_eff_mass(UseJack);

 
  distr_t_list Exch_r0_eff_slope(UseJack), Exch_r1_eff_slope(UseJack), Hand_r0_eff_slope(UseJack), Hand_r1_eff_slope(UseJack), Global_r0_eff_slope(UseJack), Global_r1_eff_slope(UseJack);
  
  for(int i=0; i < Nens; i++) {
   
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    if(pi_data_r0.Tag[i].substr(0,1) == "A")  Corr.Tmin = 5;
      else if(pi_data_r0.Tag[i].substr(0,1) =="B") Corr.Tmin = 8;
      else Corr.Tmin= 10;
    Corr.Tmax = pi_data_r0.nrows[i]/2 -12;
    Corr.Nt = pi_data_r0.nrows[i];
    LatticeInfo L_info("LOCAL");
    
    L_info.LatInfo(pi_data_r0.Tag[i].substr(0,1));
    double L, NNT, m;
    Read_pars_from_ensemble_tag(pi_data_r0.Tag[i], m, L, NNT);
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/axion_l7");
    

    //compute correlators

    Pi_r0_corr = Corr.corr_t(pi_data_r0.col(0)[i], "../data/axion_l7/pi_corr_r0."+pi_data_r0.Tag[i]);
    //Pi_r1_corr = Corr.corr_t(pi_data_r1.col(0)[i], "../data/axion_l7/pi_corr_r1."+pi_data_r0.Tag[i]);
    //Pi_0_disc_corr = Corr.corr_t(disc_pi0.col(1)[i], "../data/axion_l7/disc_pi0."+pi_data_r0.Tag[i]);
    Exch_r0_corr = Corr.corr_t(IE1_data.col(0)[i], "../data/axion_l7/Exch_r0."+pi_data_r0.Tag[i]);
    //Exch_r1_corr = Corr.corr_t(IE2_data.col(0)[i], "../data/axion_l7/Exch_r1."+pi_data_r0.Tag[i]);
    Hand_r0_corr = Corr.corr_t(IH1_data.col(0)[i], "../data/axion_l7/Hand_r0."+pi_data_r0.Tag[i]);
    //Hand_r1_corr = Corr.corr_t(IH2_data.col(0)[i], "../data/axion_l7/Hand_r1."+pi_data_r0.Tag[i]);
    //Bubble_r0_R_corr = Corr.corr_t(IB1_data.col(0)[i], "");
    IS1_corr = Corr.corr_t(IS1_data.col(0)[i], "../data/axion_l7/IS_corr_r0r0."+pi_data_r0.Tag[i]);
    

      
    //erase gluon-disconnected contribution from Hand_r0 and Hand_r1

    /*
    for(int ts=0; ts<Corr.Nt;ts++) {
      for(int tt=0;tt<Corr.Nt;tt++) {
	Hand_r0_corr.distr_list[ts] = Hand_r0_corr.distr_list[ts]- 0*(1.0/(double)Corr.Nt)*Bubble_r0_R_corr.distr_list[tt]*Bubble_r0_R_corr.distr_list[(tt+ts)%Corr.Nt] - 0*(1.0/(double)Corr.Nt)*Bubble_r0_I_corr.distr_list[tt]*Bubble_r0_I_corr.distr_list[(tt+ts)%Corr.Nt];
	Hand_r1_corr.distr_list[ts] = Hand_r1_corr.distr_list[ts]- 0*(1.0/(double)Corr.Nt)*Bubble_r1_R_corr.distr_list[tt]*Bubble_r1_R_corr.distr_list[(tt+ts)%Corr.Nt] - 0*(1.0/(double)Corr.Nt)*Bubble_r1_I_corr.distr_list[tt]*Bubble_r1_I_corr.distr_list[(tt+ts)%Corr.Nt];
    }
    }
    

    */


    

    //compute effective masses of neutral and charged pion

    Pi_plus_corr=Pi_r0_corr;
    Pi_plus_eff_mass = Corr.effective_mass_t(Pi_plus_corr, "../data/axion_l7/pi_plus_eff_mass."+pi_data_r0.Tag[i]);
    
    distr_t_list f0_distr = 2.0*m*Corr.decay_constant_t(Pi_plus_corr, "../data/axion_l7/pi_plus_decay_const."+pi_data_r0.Tag[i]);
    distr_t f0 = Corr.Fit_distr(f0_distr);
    distr_t Pi_plus_fit = Corr.Fit_distr(Pi_plus_eff_mass);
       
    
     
    distr_t_list ratio_corr_exch_r0 = Exch_r0_corr/Pi_plus_corr;
    //distr_t_list ratio_corr_exch_r1 = Exch_r1_corr/Pi_0_corr;

    Print_To_File({}, {ratio_corr_exch_r0.ave(), ratio_corr_exch_r0.err()},"../data/axion_l7/ratio_corr_exch_r0."+pi_data_r0.Tag[i]  ,"OUT", "");
    //Print_To_File({}, {ratio_corr_exch_r1.ave(), ratio_corr_exch_r1.err()},"../data/axion_l7/ratio_corr_exch_r1."+pi_data_r0.Tag[i]  ,"OUT", "");

    distr_t_list ratio_corr_hand_r0 = Hand_r0_corr/Pi_plus_corr;
    //distr_t_list ratio_corr_hand_r1 = Hand_r1_corr/Pi_0_corr;

    Print_To_File({}, {ratio_corr_hand_r0.ave(), ratio_corr_hand_r0.err()},"../data/axion_l7/ratio_corr_hand_r0."+pi_data_r0.Tag[i]  ,"OUT", "");
    //Print_To_File({}, {ratio_corr_hand_r1.ave(), ratio_corr_hand_r1.err()},"../data/axion_l7/ratio_corr_hand_r1."+pi_data_r0.Tag[i]  ,"OUT", "");

    
    distr_t_list ratio_IS1;
    ratio_IS1 = IS1_corr/Pi_plus_corr;
   
    distr_t_list eff_mass_IS1;

    eff_mass_IS1 = Corr.effective_slope_t(IS1_corr, Pi_plus_corr, "../data/axion_l7/IS_r0r0_eff_slope."+pi_data_r0.Tag[i]);
   
    Print_To_File({}, {ratio_IS1.ave(), ratio_IS1.err()},"../data/axion_l7/ratio_corr_IS_r0r0."+pi_data_r0.Tag[i]  ,"OUT", "");
    //Print_To_File({}, {ratio_IS2.ave(), ratio_IS2.err()},"../data/axion_l7/ratio_corr_IS_r0r1."+pi_data_r0.Tag[i]  ,"OUT", "");
    //Print_To_File({}, {ratio_IS3.ave(), ratio_IS3.err()},"../data/axion_l7/ratio_corr_IS_r1r0."+pi_data_r0.Tag[i]  ,"OUT", "");
    //Print_To_File({}, {ratio_IS4.ave(), ratio_IS4.err()},"../data/axion_l7/ratio_corr_IS_r1r1."+pi_data_r0.Tag[i]  ,"OUT", "");
     
  
    
   
    distr_t ZP_jack(UseJack), ZS_jack(UseJack);
    GaussianMersenne Gauss_ZP(43423,L_info.ZP_err/sqrt(Njacks-1));
    GaussianMersenne Gauss_ZS(65443,L_info.ZS_err/sqrt(Njacks-1));
    for(int nj=0;nj<Njacks;nj++) ZP_jack.distr.push_back( L_info.ZP +Gauss_ZP());
    for(int nj=0;nj<Njacks;nj++) ZS_jack.distr.push_back(L_info.ZS +Gauss_ZS()); 
    
    //compute effective slopes
    //distr_t_list Exch_r0_eff_slope_2nd= ZP_jack*ZP_jack*Corr.effective_slope_t_2nd_ord(Exch_r0_corr, Pi_plus_corr, IS1_corr,  "../data/axion_l7/Exch_r0_eff_slope_2nd."+pi_data_r0.Tag[i]);

    //compute A'M' term

    
    
    distr_t_list A_prime_M_prime = eff_mass_IS1*Corr.residue_correction_t(IS1_corr, Pi_plus_corr, "../data/axion_l7/residue_correction."+pi_data_r0.Tag[i]);

     
    Exch_r0_eff_slope = Corr.effective_slope_t(Exch_r0_corr, Pi_plus_corr, "../data/axion_l7/Exch_r0_eff_slope."+pi_data_r0.Tag[i]);

    Hand_r0_eff_slope = Corr.effective_slope_t(Hand_r0_corr, Pi_plus_corr, "../data/axion_l7/Hand_r0_eff_slope."+pi_data_r0.Tag[i]);

    Global_r0_eff_slope = Corr.effective_slope_t(Exch_r0_corr-Hand_r0_corr, Pi_plus_corr, "../data/axion_l7/Global_r0_eff_slope."+pi_data_r0.Tag[i]);

    //fit Global_r0_eff_slope
    distr_t Fit_result = Corr.Fit_distr(Global_r0_eff_slope);

    distr_t_list ell_7= Global_r0_eff_slope*(ZS_jack*ZS_jack/(ZP_jack*ZP_jack))*f0*f0*(1.0/(Pi_plus_fit*Pi_plus_fit*Pi_plus_fit))*(m*m);

    Print_To_File({}, {ell_7.ave(), ell_7.err()}, "../data/axion_l7/ell_7."+pi_data_r0.Tag[i]+".t", "OUT", "");
    
    //plots
    plt::clf();
    plt::figure_size(1000,1000);
    plt::xlim(1.0, (double)Corr.Nt/2.0 +1);
    plt::xlabel("t");
    plt::ylabel("effective slope");
    Vfloat times;
    for(int tt=0;tt<Corr.Nt;tt++) times.push_back(tt);
    plt::errorbar(times, Exch_r0_eff_slope.ave(),Exch_r0_eff_slope.err(), { {"c", "black"}, {"marker", "."} , {"ls" , "-"}, {"label", "connected"}});
    plt::errorbar(times, Global_r0_eff_slope.ave(), Global_r0_eff_slope.err(), { {"c", "red"}, {"marker", "."} , {"ls" , "-"}, {"label", "global"}});
    plt::errorbar(times, Hand_r0_eff_slope.ave(), Hand_r0_eff_slope.err(), { {"c", "green"}, {"marker", "."}, {"ls", "-"}, {"label", "disconnected"}});
    plt::legend();
    string path="../plots/axion_l7/eff_slopes_"+pi_data_r0.Tag[i]+".png";
    plt::save(path.c_str());
    plt::close();
  
    
    cout<<"#####################"<<endl;
    cout<<"Ensemble: "<<pi_data_r0.Tag[i]<<endl;
    //try to reconstruct l7
    cout<<"m: "<<m<<endl;
    cout<<"Zp: "<<L_info.ZP<<endl;
    cout<<"Zs: "<<L_info.ZS<<endl;
    cout<<"pi mass: "<<Pi_plus_fit.ave()<<endl;
    cout<<"f0: "<<f0.ave()<<endl;
    cout<<"dm: "<<Fit_result.ave()<<"("<<Fit_result.err()<<")"<<endl;
    distr_t l7 =  (ZS_jack*ZS_jack/(ZP_jack*ZP_jack))*f0*f0*Fit_result*(1.0/(Pi_plus_fit*Pi_plus_fit*Pi_plus_fit))*(m*m);
   
    cout<<"l7: "<<l7.ave()<<"("<<l7.err()<<")"<<endl;
    cout<<"#####################"<<endl;
     
  }


  
     
  return;
}
 
					
  






















