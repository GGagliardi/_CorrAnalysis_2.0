#include "../include/Axion_l7.h"


using namespace std;
namespace plt = matplotlibcpp;

const double MPiPhys=0.135;
const bool Use_JB_distribution= false;
const bool UseJack=1;
const int Njacks=20;
const int Nboots=200;









void Axion_l7_analysis() {
  
  data_t pi_data_r0,pi_data_r1, IE1_data, IE2_data, IH1_data, IH2_data, disc_pi0;
  data_t IB1_data, IB2_data, IS1_data, IS2_data, IS3_data, IS4_data;
  pi_data_r0.Read("../datasets_axion", "mes_contr_00_r0", "P5P5");
  //  pi_data_r1.Read("../datasets_axion", "mes_contr_00_r1", "P5P5");
  disc_pi0.Read("../datasets_axion", "mes_contr_disc_pi0", "S0P5");
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
  
  distr_t_list Pi_0_disc_corr(UseJack), Pi_0_corr(UseJack), Pi_plus_corr(UseJack);

  distr_t_list Pi_0_eff_mass(UseJack), Pi_plus_eff_mass(UseJack);

  distr_t_list Exch_r0_eff_slope(UseJack), Exch_r1_eff_slope(UseJack), Hand_r0_eff_slope(UseJack), Hand_r1_eff_slope(UseJack), Global_r0_eff_slope(UseJack), Global_r1_eff_slope(UseJack);
  
  for(int i=0; i < Nens; i++) {
   
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    if(pi_data_r0.Tag[i].substr(0,1) == "A")  Corr.Tmin = 15;
      else if(pi_data_r0.Tag[i].substr(0,1) =="B") Corr.Tmin = 13;
      else Corr.Tmin= 19;
    Corr.Tmax = pi_data_r0.nrows[i]/2 -3;
    Corr.Nt = pi_data_r0.nrows[i];
    LatticeInfo L_info("LOCAL");
    
    L_info.LatInfo(pi_data_r0.Tag[i].substr(0,1));
    double Zp= L_info.ZP;
    double L, NNT, m;
    Read_pars_from_ensemble_tag(pi_data_r0.Tag[i], m, L, NNT);
    boost::filesystem::create_directory("../data");
    boost::filesystem::create_directory("../data/axion_l7");
    


    //compute correlators

    Pi_r0_corr = Corr.corr_t(pi_data_r0.col(0)[i], "../data/axion_l7/pi_corr_r0."+pi_data_r0.Tag[i]);
    //Pi_r1_corr = Corr.corr_t(pi_data_r1.col(0)[i], "../data/axion_l7/pi_corr_r1."+pi_data_r0.Tag[i]);
    //Pi_0_disc_corr = Corr.corr_t(disc_pi0.col(0)[i], "../data/axion_l7/disc_pi0."+pi_data_r0.Tag[i]);
    Exch_r0_corr = Corr.corr_t(IE1_data.col(0)[i], "../data/axion_l7/Exch_r0."+pi_data_r0.Tag[i]);
    //Exch_r1_corr = Corr.corr_t(IE2_data.col(0)[i], "../data/axion_l7/Exch_r1."+pi_data_r0.Tag[i]);
    Hand_r0_corr = Corr.corr_t(IH1_data.col(0)[i], "../data/axion_l7/Hand_r0."+pi_data_r0.Tag[i]);
    //Hand_r1_corr = Corr.corr_t(IH2_data.col(0)[i], "../data/axion_l7/Hand_r1."+pi_data_r0.Tag[i]);
    //Bubble_r0_R_corr = Corr.corr_t(IB1_data.col(0)[i], "");
    //Bubble_r1_R_corr = Corr.corr_t(IB2_data.col(0)[i], "");
    //Bubble_r0_I_corr = Corr.corr_t(IB1_data.col(1)[i], "");
    //Bubble_r1_I_corr = Corr.corr_t(IB2_data.col(1)[i], "");
    IS1_corr = Corr.corr_t(IS1_data.col(0)[i], "../data/axion_l7/IS_corr_r0r0."+pi_data_r0.Tag[i]);
    //IS2_corr = Corr.corr_t(IS2_data.col(0)[i], "../data/axion_l7/IS_corr_r0r1."+pi_data_r0.Tag[i]);
    //IS3_corr = Corr.corr_t(IS3_data.col(0)[i], "../data/axion_l7/IS_corr_r1r0."+pi_data_r0.Tag[i]);
    //IS4_corr = Corr.corr_t(IS4_data.col(0)[i], "../data/axion_l7/IS_corr_r1r1."+pi_data_r0.Tag[i]);


   
    //erase gluon-disconnected contribution from Hand_r0 and Hand_r1

    /*
    for(int ts=0; ts<Corr.Nt;ts++) {
      for(int tt=0;tt<Corr.Nt;tt++) {
	Hand_r0_corr.distr_list[ts] = Hand_r0_corr.distr_list[ts]- 0*(1.0/(double)Corr.Nt)*Bubble_r0_R_corr.distr_list[tt]*Bubble_r0_R_corr.distr_list[(tt+ts)%Corr.Nt] - 0*(1.0/(double)Corr.Nt)*Bubble_r0_I_corr.distr_list[tt]*Bubble_r0_I_corr.distr_list[(tt+ts)%Corr.Nt];
	Hand_r1_corr.distr_list[ts] = Hand_r1_corr.distr_list[ts]- 0*(1.0/(double)Corr.Nt)*Bubble_r1_R_corr.distr_list[tt]*Bubble_r1_R_corr.distr_list[(tt+ts)%Corr.Nt] - 0*(1.0/(double)Corr.Nt)*Bubble_r1_I_corr.distr_list[tt]*Bubble_r1_I_corr.distr_list[(tt+ts)%Corr.Nt];
    }
    }
    
    Pi_0_corr = Pi_r1_corr - Pi_0_disc_corr;

    */


    

    //compute effective masses of neutral and charged pion

    Pi_plus_corr=Pi_r0_corr;
      
    Pi_plus_eff_mass = Corr.effective_mass_t(Pi_plus_corr, "../data/axion_l7/pi_plus_eff_mass."+pi_data_r0.Tag[i]);
    distr_t_list f0_distr = 2.0*m*Corr.decay_constant_t(Pi_plus_corr, "../data/axion_l7/pi_plus_decay_const."+pi_data_r0.Tag[i]);
    distr_t f0 = Corr.Fit_distr(f0_distr);
    distr_t Pi_plus_fit = Corr.Fit_distr(Pi_plus_eff_mass);
       
    //Pi_0_eff_mass = Corr.effective_mass_t(Pi_0_corr, "../data/axion_l7/pi_zero_eff_mass."+pi_data_r0.Tag[i]);

    
    distr_t_list ratio_corr_exch_r0 = Exch_r0_corr/Pi_plus_corr;
    //distr_t_list ratio_corr_exch_r1 = Exch_r1_corr/Pi_0_corr;

    Print_To_File({}, {ratio_corr_exch_r0.ave(), ratio_corr_exch_r0.err()},"../data/axion_l7/ratio_corr_exch_r0."+pi_data_r0.Tag[i]  ,"OUT", "");
    //Print_To_File({}, {ratio_corr_exch_r1.ave(), ratio_corr_exch_r1.err()},"../data/axion_l7/ratio_corr_exch_r1."+pi_data_r0.Tag[i]  ,"OUT", "");

    distr_t_list ratio_corr_hand_r0 = Hand_r0_corr/Pi_plus_corr;
    //distr_t_list ratio_corr_hand_r1 = Hand_r1_corr/Pi_0_corr;

    Print_To_File({}, {ratio_corr_hand_r0.ave(), ratio_corr_hand_r0.err()},"../data/axion_l7/ratio_corr_hand_r0."+pi_data_r0.Tag[i]  ,"OUT", "");
    //Print_To_File({}, {ratio_corr_hand_r1.ave(), ratio_corr_hand_r1.err()},"../data/axion_l7/ratio_corr_hand_r1."+pi_data_r0.Tag[i]  ,"OUT", "");




    distr_t_list ratio_IS1, ratio_IS2, ratio_IS3, ratio_IS4;
    ratio_IS1 = IS1_corr/Pi_plus_corr;
    //ratio_IS2 = IS2_corr/Pi_plus_corr;
    //ratio_IS3 = IS3_corr/Pi_plus_corr;
    //ratio_IS4 = IS4_corr/Pi_plus_corr;

    distr_t_list eff_mass_IS1, eff_mass_IS2, eff_mass_IS3, eff_mass_IS4;

    eff_mass_IS1 = Corr.effective_slope_t(IS1_corr, Pi_plus_corr, "../data/axion_l7/IS_r0r0_eff_slope."+pi_data_r0.Tag[i]);
    //eff_mass_IS2 = Corr.effective_slope_t(IS2_corr, Pi_plus_corr, "../data/axion_l7/IS_r0r1_eff_slope."+pi_data_r0.Tag[i]);
    //eff_mass_IS3 = Corr.effective_slope_t(IS3_corr, Pi_plus_corr, "../data/axion_l7/IS_r1r0_eff_slope."+pi_data_r0.Tag[i]);
    //eff_mass_IS4 = Corr.effective_slope_t(IS4_corr, Pi_plus_corr, "../data/axion_l7/IS_r1r1_eff_slope."+pi_data_r0.Tag[i]);

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
   
    
    Exch_r0_eff_slope = Corr.effective_slope_t(ZP_jack*ZP_jack*Exch_r0_corr, Pi_plus_corr, "../data/axion_l7/Exch_r0_eff_slope."+pi_data_r0.Tag[i]);
    //Exch_r1_eff_slope = Corr.effective_slope_t(Exch_r1_corr, Pi_0_corr, "../data/axion_l7/Exch_r1_eff_slope."+pi_data_r0.Tag[i]);
    //distr_t_list Hand_r0_eff_slope_2nd = ZS_jack*ZS_jack*Corr.effective_slope_t_2nd_ord( Hand_r0_corr, Pi_plus_corr, IS1_corr, "../data/axion_l7/Hand_r0_eff_slope_2nd."+pi_data_r0.Tag[i]);
    Hand_r0_eff_slope = Corr.effective_slope_t(Hand_r0_corr, Pi_plus_corr, "../data/axion_l7/Hand_r0_eff_slope."+pi_data_r0.Tag[i]);
    //Hand_r1_eff_slope = Corr.effective_slope_t(Hand_r1_corr, Pi_0_corr, "../data/axion_l7/Hand_r1_eff_slope."+pi_data_r0.Tag[i]);
    Global_r0_eff_slope = Exch_r0_eff_slope -1.0*(ZS_jack*ZS_jack)*Hand_r0_eff_slope;
    distr_t_list Global_r0_eff_slope_2nd = Corr.effective_slope_t_2nd_ord(ZP_jack*ZP_jack*Exch_r0_corr, Pi_plus_corr, Hand_r0_corr, "../data/axion_l7/subtracted_eff_slope."+pi_data_r0.Tag[i]);
    //Print_To_File({}, {Exch_r0_eff_slope_2nd.ave(), Exch_r0_eff_slope_2nd.err()},"../data/axion_l7/Exch_r0_eff_slope_renorm_2nd."+pi_data_r0.Tag[i]  ,"OUT", "");
    Print_To_File({}, {Global_r0_eff_slope.ave(), Global_r0_eff_slope.err()},"../data/axion_l7/Global_r0_eff_slope."+pi_data_r0.Tag[i]  ,"OUT", "");
    //Print_To_File({}, {Global_r0_eff_slope_2nd.ave(), Global_r0_eff_slope_2nd.err()},"../data/axion_l7/Global_r0_eff_slope_2nd."+pi_data_r0.Tag[i]  ,"OUT", "");
    //Global_r1_eff_slope = Corr.effective_slope_t(2.0*(Exch_r1_corr - Hand_r1_corr), Pi_0_corr, "../data_axion_l7/l7_r1_eff_slope."+pi_data_r0.Tag[i]);
    //fit Global_r0_eff_slope
    distr_t Fit_result = Corr.Fit_distr(Global_r0_eff_slope);
    distr_t Fit_result_sub = Corr.Fit_distr(Global_r0_eff_slope_2nd);

    //plots
    plt::clf();
    plt::figure_size(1000,1000);
    plt::xlim(1.0, (double)Corr.Nt/2.0 +1);
    plt::xlabel("t");
    //plt::ylim(0.0,3000.0); 
    plt::ylabel("effective slope");
    Vfloat times;
    for(int tt=0;tt<Corr.Nt;tt++) times.push_back(tt);
    plt::errorbar(times, Exch_r0_eff_slope.ave(),Exch_r0_eff_slope.err(), { {"c", "black"}, {"marker", "."} , {"ls" , "-"}, {"label", "exchange"}});
    plt::errorbar(times, Global_r0_eff_slope.ave(), Global_r0_eff_slope.err(), { {"c", "red"}, {"marker", "."} , {"ls" , "-"}, {"label", "exchange-disconnected"}});
    //plt::errorbar(times, Exch_r0_eff_slope_2nd.ave(),Exch_r0_eff_slope_2nd.err(), { {"c", "blue"}, {"marker", "."} , {"ls" , "-"}, {"label", "exchange 2nd ord"}});
    plt::errorbar(times, Global_r0_eff_slope_2nd.ave(), Global_r0_eff_slope_2nd.err(), { {"c", "green"}, {"marker", "."} , {"ls" , "-"}, {"label", "exchange-disconnected subtracted"}});
    plt::legend();
    string path="../plots/axion_l7/eff_slope_tanh_"+pi_data_r0.Tag[i]+".png";
    plt::save(path.c_str());
    /*
    plt::clf();
    plt::figure_size(1000,1000);
    plt::xlim(1.0, (double)Corr.Nt/2.0 +1);
    plt::ylim(0.0,2050.0); 
    plt::xlabel("t");
    plt::ylabel("effective slope");
    plt::errorbar(times, Exch_r0_eff_slope_2nd.ave(),Exch_r0_eff_slope_2nd.err(), { {"c", "black"}, {"marker", "."} , {"ls" , "-"}, {"label", "exchange"}});
    plt::errorbar(times, Global_r0_eff_slope_2nd.ave(), Global_r0_eff_slope_2nd.err(), { {"c", "red"}, {"marker", "."} , {"ls" , "-"}, {"label", "exchange-disconnected"}});
    plt::legend();
    path="../plots/axion_l7/eff_slope_subtracted"+pi_data_r0.Tag[i]+ ".png";
    plt::save(path.c_str());
    */
    
    
    cout<<"Ensemble: "<<pi_data_r0.Tag[i]<<" dm: "<<Fit_result.ave()<<"("<<Fit_result.err()<<")"<<endl;
    //try to reconstruct l7
    cout<<"m: "<<m<<endl;
    cout<<"Zp: "<<Zp<<endl;
    cout<<"pi mass: "<<Pi_plus_fit.ave()<<endl;
    cout<<"f0: "<<f0.ave()<<endl;
    distr_t l7 = (2.0*f0*f0*Fit_result*(1.0/(Pi_plus_fit*Pi_plus_fit*Pi_plus_fit)))*(m*m);
    distr_t l7_v2 = l7/(ZP_jack*ZP_jack);
    distr_t l7_sub = (2.0*f0*f0*Fit_result_sub*(1.0/(Pi_plus_fit*Pi_plus_fit*Pi_plus_fit)))*(m*m);
    distr_t l7_sub_v2 = l7_sub/(ZP_jack*ZP_jack);

    cout<<"l7: "<<l7.ave()<<"("<<l7.err()<<")"<<endl;
    cout<<"l7_sub: "<<l7_sub.ave()<<"("<<l7_sub.err()<<")"<<endl;
    cout<<"l7_v2: "<<l7_v2.ave()<<"("<<l7_v2.err()<<")"<<endl;
    cout<<"l7_sub_v2: "<<l7_sub_v2.ave()<<"("<<l7_sub_v2.err()<<")"<<endl;
     
  }


  
     
  return;
}
 
					
  






















