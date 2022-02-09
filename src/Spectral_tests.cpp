#include "../include/Spectral_tests.h"



const double MPiPhys=0.135;
const double alpha = 1.0/137.04;
const double e2 = alpha*4.0*M_PI;
const bool UseJack=1;
const int Njacks=50;
const int Nboots=200;
const double qu= 2.0/3.0;
const double qd= -1.0/3.0;
const double qs= qd;
const double qc= qu;
const double fm_to_inv_Gev= 1.0/0.1975;
const int prec = 128;
const double rho_R= 12*M_PI*M_PI;
const double Nc=3;
const double Rpert= Nc*( qu*qu + qd*qd);
const double m_mu= 0.10565837; //GeV
const string SM_TYPE= "GAUSSIAN";

using namespace std;


void Spectral_tests() {

 
  //create output directories
  boost::filesystem::create_directory("../data/R_ratio");
  boost::filesystem::create_directory("../data/R_ratio/corr");

  //LOAD GM2 DATA
  GaussianMersenne GM(981832);

  
  data_t  V_light_1, V_light_2, V_light_3;
  data_t  V_light_OS_1, V_light_OS_2, V_light_OS_3;

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
  //#################################END READING########################




  //###########################################################################################


  //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D
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


  //############################################################################################

  int Nens_light =  V_light_1.size;

  for(int i_ens=0;i_ens<Nens_light;i_ens++) {

 
    
    CorrAnalysis Corr(UseJack, Njacks,Nboots);
    Corr.Nt = V_light_1.nrows[i_ens];
    int T = Corr.Nt;

    cout<<"Analyzing Ensemble: "<<V_light_1.Tag[i_ens]<<endl;
    cout<<"NT: "<<T<<endl;

    //get lattice spacing
    distr_t a_distr(UseJack);
    if(V_light_1.Tag[i_ens].substr(1,1)=="A") {a_distr=a_A;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="B") {a_distr=a_B;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="C") {a_distr=a_C;}
    else if(V_light_1.Tag[i_ens].substr(1,1)=="D") {a_distr=a_D;}
    else crash("lattice spacing distribution for Ens: "+V_light_1.Tag[i_ens]+" not found");
  
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(V_light_1.Tag[i_ens]);
 
   
    distr_t_list  V_light_1_distr, V_light_2_distr, V_light_3_distr, disco_distr;
    distr_t_list  V_light_OS_1_distr, V_light_OS_2_distr, V_light_OS_3_distr, V_light_OS_distr;  
    distr_t_list  V_light_distr; 
    distr_t Zv(UseJack), Za(UseJack);
  

    //vector light sector
    V_light_1_distr = Corr.corr_t(V_light_1.col(0)[i_ens], "");
    V_light_2_distr = Corr.corr_t(V_light_2.col(0)[i_ens], "");
    V_light_3_distr = Corr.corr_t(V_light_3.col(0)[i_ens], "");
    V_light_OS_1_distr = Corr.corr_t(V_light_OS_1.col(0)[i_ens], "");
    V_light_OS_2_distr = Corr.corr_t(V_light_OS_2.col(0)[i_ens], "");
    V_light_OS_3_distr = Corr.corr_t(V_light_OS_3.col(0)[i_ens], "");
    //sum over the Lorentz indices of the e.m. current
    V_light_distr= ((pow(qu,2)+pow(qd,2))/3.0)*(V_light_1_distr+ V_light_2_distr + V_light_3_distr);
    V_light_OS_distr =  ((pow(qu,2)+pow(qd,2))/3.0)*(V_light_OS_1_distr+ V_light_OS_2_distr + V_light_OS_3_distr);


    //############################  GET RENORMALIZATION CONSTANT FROM HADRONIC METHOD ################################################




    //resample lattice spacing
    //generate jackknife sample of input parameters
  
    if(UseJack)  { for(int ijack=0;ijack<Njacks;ijack++) {
	Za.distr.push_back(  L_info.Za_WI + GM()*L_info.Za_WI_err*(1.0/sqrt(Njacks-1.0)));
	Zv.distr.push_back(  L_info.Zv_WI + GM()*L_info.Zv_WI_err*(1.0/sqrt(Njacks-1.0)));
      }
    }
    else {
      for (int iboot=0; iboot<Nboots;iboot++) {
	Za.distr.push_back(  L_info.Za_WI + GM()*L_info.Za_WI_err);
	Zv.distr.push_back(  L_info.Zv_WI + GM()*L_info.Zv_WI_err);
      }
    }



    //############################################################################################


    //leading order perturbative substraction of a^2/t^2 lattice artifacts
  
    //free corr LO artifacts
    Vfloat free_corr_log_art(Corr.Nt);
    for(int t=0;t<Corr.Nt;t++) { if( t*a_distr.ave()/fm_to_inv_Gev < 2.0 && t != 0) {   free_corr_log_art[t] = -1.0*(qu*qu + qd*qd)*(1.0/(2.0*M_PI*M_PI*pow(t,5)));} else free_corr_log_art[t] = 0.0;}


    //multiply corr using Zv and Za
    V_light_distr = (V_light_distr)*Za*Za;
    V_light_OS_distr = (V_light_OS_distr)*Zv*Zv;


    //print correlator
  
    Print_To_File({}, {V_light_distr.ave(), V_light_distr.err(), V_light_OS_distr.ave(), V_light_OS_distr.err()}, "../data/R_ratio/corr/corr_"+V_light_1.Tag[i_ens]+".dat", "", "# t  tm  OS");


    cout<<"Perturbative corrections added"<<endl;





    //#########################   RECONSTRUCT THE SMEARED R-RATIO ################################


    //set tmax to the value where the error on V(t) is larger than x%
    //#############################################################################################################################

    bool Found_error_less_x_percent=false;
    double x=5;
    double tmax=1;
    while(!Found_error_less_x_percent && tmax < Corr.Nt/2) {
   
      if( (V_light_distr.distr_list[tmax]).err()/fabs( (V_light_distr.distr_list[tmax]).ave()) <  0.01*x) tmax++;
      else Found_error_less_x_percent=true;
    }

    Found_error_less_x_percent=false;
    double tmax_OS=1;
    while(!Found_error_less_x_percent && tmax < Corr.Nt/2) {
   
      if( (V_light_OS_distr.distr_list[tmax_OS]).err()/fabs( (V_light_OS_distr.distr_list[tmax_OS]).ave()) <  0.01*x) tmax_OS++;
      else Found_error_less_x_percent=true;
    }
    
    //#############################################################################################################################

  

    //##################################################
    //#############################SET GENERAL PARAMETERS######################################
    double Emax= 4.0*a_distr.ave();  //[x GeV]*a
    double E0 =  2*MPiPhys*a_distr.ave(); // two pion treshold ~ 280 MeV
    double Estart= 2*MPiPhys*a_distr.ave(); // start from 300 MeV
    double sigma = 1.0*a_distr.ave(); // smearing
    double step_size= 0.25; //step size in energy in units of sigma
    int Npoints = (int)((Emax - Estart)/(sigma*step_size));
    tmax= (T/2 -4);
    tmax_OS = T/2 -4;
    //#########################################################################################


    cout<<"Reconstructing R(E)....."<<endl;
    cout<<"aEmax: "<<Emax<<" -> Emax: "<<Emax/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"aE0: "<<E0<<" -> E0: "<<E0/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"asigma: "<<sigma<<" -> sigma: "<<sigma/a_distr.ave()<<" [GeV] "<<endl;
    cout<<"Npoints: "<<Npoints<<endl;
    cout<<"Corr(0) tm: "<<(rho_R*V_light_distr).ave(0)<<endl;
    cout<<"Corr(0) OS: "<<(rho_R*V_light_OS_distr).ave(0)<<endl;
    cout<<"Zv: "<<L_info.Zv_WI<< " +- "<<L_info.Zv_WI_err<<endl;
    cout<<"Za: "<<L_info.Za_WI<< " +- "<<L_info.Za_WI_err<<endl;
    cout<<"SM_TYPE: "<<SM_TYPE<<endl;
    

    // lambda function to be used as a smearing func.

    const auto f = [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0) -> PrecFloat {

		     if(SM_TYPE == "WIN")  return (1.0/(1 + exp(-2*(E-(m-s))/s)) - 1.0/(1 + exp(-2*(E-(m+s))/s)))/(2*s*sqr(E));

		     else if (SM_TYPE == "GAUSSIAN") return Get_exact_gauss(E, m, s, E0)/sqr(E);

		     else crash("SM_TYPE: "+SM_TYPE+" not yet implemented");

		     return PrecFloat(0);
		   
		   };  

 

    distr_t_list Spectral_dens;
    distr_t_list Spectral_dens_OS;
    distr_t_list Spectral_dens_SANF;
    distr_t_list Spectral_dens_OS_SANF;
    Vfloat Ergs;
    Vfloat Ergs_GeV;
    Vfloat lambdas_tm;
    Vfloat lambdas_OS;
    Vfloat syst_tm(Npoints);
    Vfloat syst_OS(Npoints);
    Vfloat syst_tm_SANF(Npoints);
    Vfloat syst_OS_SANF(Npoints);


    //COMPUTE THE SMEARED R-RATIO

    for(int ip=0; ip<Npoints;ip++) {

      double mean = Estart + ip*((Emax-Estart)/((double)Npoints));
      Ergs.push_back(mean);
      Ergs_GeV.push_back(mean/a_distr.ave());
      double sigma_resc=sigma;
      double lambda_Estar;
      double lambda_Estar_SANF;
      //tmax = (int)((T/2-2)*(Estart/mean));
      //tmax_OS = (int)((T/2-2)*(Estart/mean));
      //double sigma_resc= (mean/a_distr.ave() > 0.8)?sigma*sqrt(mean/guess_resonance_pos):sigma;

 
      cout<<"Evaluating R(E*) at aE*: "<<mean<<" -> "<<"E*: "<<mean/a_distr.ave()<<" [GeV] ..."<<endl;
      cout<<"Rescaled asigma(E*): "<<sigma_resc<<" -> sigma(E*): "<<sigma_resc/a_distr.ave()<<" [GeV] "<<endl;
      cout<<"[tmin, tmax] (tm): ["<<1<<", "<<tmax<<"]"<<endl;
      cout<<"[tmin, tmax] (OS): ["<<1<<", "<<tmax_OS<<"]"<<endl;
      Spectral_dens.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax, prec, SM_TYPE+"_ov_E2",f, rho_R*V_light_distr, syst_tm[ip], Rpert, lambda_Estar, "TANT", "tm", "VV_"+V_light_1.Tag[i_ens] )) ;
      cout<<"R(E*) tm (TANTALO): "<<Spectral_dens.ave(ip)<<" +- ("<<Spectral_dens.err(ip)<<")_stat ("<<syst_tm[ip]<<")_syst "<<endl;
      Spectral_dens_SANF.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax, prec, SM_TYPE+"_ov_E2",f, rho_R*V_light_distr, syst_tm_SANF[ip], Rpert, lambda_Estar_SANF, "SANF", "tm", "VV_"+V_light_1.Tag[i_ens] )) ;
      cout<<"R(E*) tm (SANFILIPPO): "<<Spectral_dens_SANF.ave(ip)<<" +- ("<<Spectral_dens_SANF.err(ip)<<")_stat ("<<syst_tm_SANF[ip]<<")_syst "<<endl;
      lambdas_tm.push_back(lambda_Estar);
    
      Spectral_dens_OS.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS, prec, SM_TYPE+"_ov_E2",f, rho_R*V_light_OS_distr, syst_OS[ip], Rpert, lambda_Estar, "TANT", "OS", "VV_"+V_light_1.Tag[i_ens] )) ;
      cout<<"R(E*) OS (TANTALO): "<<Spectral_dens_OS.ave(ip)<<" +- ("<<Spectral_dens_OS.err(ip)<<")_stat ("<<syst_OS[ip]<<")_syst "<<endl;
      Spectral_dens_OS_SANF.distr_list.push_back(Get_Laplace_transfo(  mean,  sigma_resc, E0,  T, tmax_OS, prec, SM_TYPE+"_ov_E2",f, rho_R*V_light_OS_distr, syst_OS_SANF[ip], Rpert, lambda_Estar_SANF, "SANF", "OS", "VV_"+V_light_1.Tag[i_ens] )) ;
      cout<<"R(E*) OS (SANFILIPPO): "<<Spectral_dens_OS_SANF.ave(ip)<<" +- ("<<Spectral_dens_OS_SANF.err(ip)<<")_stat ("<<syst_OS_SANF[ip]<<")_syst "<<endl;
      lambdas_OS.push_back(lambda_Estar);
     
      cout<<"...done"<<endl;
      //############################################################################################



    }

  
    cout<<"printing output..."<<endl;


    //print to file
    Print_To_File({}, {Ergs, Ergs_GeV, lambdas_tm,  Spectral_dens.ave(), Spectral_dens.err(), syst_tm, Spectral_dens_SANF.ave(), Spectral_dens_SANF.err(), syst_tm_SANF,  lambdas_OS, Spectral_dens_OS.ave(), Spectral_dens_OS.err(),  syst_OS, Spectral_dens_OS_SANF.ave(), Spectral_dens_OS_SANF.err(), syst_OS_SANF}, "../data/R_ratio/"+SM_TYPE+"_"+V_light_1.Tag[i_ens]+"_sigma_"+to_string_with_precision(sigma,3)+".dat", "", "# aE* E*(GeV)   R(E)_tm (lambda <> stat. syst.) [TANTALO]  R(E)_tm (<> stat. syst.) [SANFILIPPO]   R(E)_OS (lambda <> stat. syst.) [TANTALO]  R(E)_OS (<> stat. syst.) [SANFILIPPO]");

    cout<<"done!"<<endl;


  }
  
  return;
}
