#include "../include/Kl4_HLT.h"
#include "input.h"
#include "numerics.h"


using namespace std;


const bool UseJack=true;
const int Njacks=35;
const string SM_FF = "FF_Exp";
//const vector<double> sigmas({0.5, 0.4, 0.3, 0.2}); 
//const vector<double> xks({0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0});
const vector<double> sigmas({1e-3});
const vector<double> xks({0.2});
const int Nsigmas= sigmas.size();
const int Nergs = xks.size();
const int prec = 128;
const int fmtIgev = 1.0 / 0.197327;
const double MK= 0.493677;



void Kl4_HLT() {


  //###########################################################
  auto K_RE= [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {
    if(SM_FF=="FF_Gauss") {
      PrecFloat x= (E-m);
      PrecFloat cosh_ov_sinh_half= (exp(x) + exp(-3*x))/(1-exp(-2*x)); 
      return  2*cos(s)/(cosh_ov_sinh_half - cos(2*s)/sinh(x));
    }
    if(SM_FF=="FF_Exp") {
      PrecFloat x= (E-m);
      //PrecFloat cosh_ov_sinh_half= (exp(x) + exp(-3*x))/(1-exp(-2*x)); 
      //return  2*cos(s)/(cosh_ov_sinh_half - cos(2*s)/sinh(x));
      return ( cos(s)*exp(-x)-exp(-2*x))/( 1 + exp(-2*x) -2*cos(s)*exp(-x));
    }
    exit(-1);
    return PrecFloat(0);
  };
  
  
  auto K_IM = [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {
    if(SM_FF=="FF_Gauss") return precPi()*Get_exact_gauss(E, m, s, E0);
    else if(SM_FF=="FF_Exp") {
      PrecFloat x= (E-m);
      return (exp(-x)*sin(s))/( 1 + exp(-2*x) -2*cos(s)*exp(-x));
      //PrecFloat cosh_ov_cosh_half= (exp(x) + exp(-3*x))/(1+exp(-2*x)); 
      //return 2*sin(s)/(cosh_ov_cosh_half - cos(2*s)/cosh(x));
    }
    return PrecFloat(0);
  };
  //###########################################################################
  

  vector<string> EnsList({"cB211b.072.64"});


  for(int iens=0; iens < (signed)EnsList.size(); iens++) {


    string Ens= EnsList[iens];

    //read from file

    int tw=25;
    //Lattice info
    LatticeInfo L_info;
    L_info.LatInfo_new_ens(Ens);
    int T= L_info.T;
    double a= L_info.a_from_afp_FLAG*fmtIgev;
    double xg=0.5;

    distr_t_list C(UseJack);
    distr_t_list C_bt(UseJack);
    
    for(int t=tw;t<T;t++) C.distr_list.emplace_back( UseJack,  Read_From_File("../data_Kl4_HLT/"+Ens+"/xg_0.5/xg_0.5/t_"+to_string(t), 0, 1) );
    for(int t=tw;t<T;t++) C_bt.distr_list.emplace_back( UseJack, Read_From_File("../data_Kl4_HLT/"+Ens+"/xg_0.5/xg_0.5_bs1/t_"+to_string(t), 0, 1) );

    boost::filesystem::create_directory("../data/Kl4_HLT");
    boost::filesystem::create_directory("../data/Kl4_HLT/"+Ens);
    boost::filesystem::create_directory("../data/Kl4_HLT/"+Ens+"/xg_0.5");

    Print_To_File({}, { C.ave(), C.err()}, "../data/Kl4_HLT/"+Ens+"/xg_0.5/C.dat", "", "");

    //Generate Covariance matrix
    Vfloat cov_T;
    for(int tt=0;tt< T-tw;tt++)
	for(int rr=0;rr< T-tw;rr++) {
	  cov_T.push_back( (C_bt.distr_list[tt]%C_bt.distr_list[rr])*(C.err(tt)*C.err(rr)/(C_bt.err(tt)*C_bt.err(rr))));
	}
																     

    //we now use HLT method
    vector<distr_t_list> RE_F , IM_F;
    for(int is=0;is<Nsigmas;is++) {
      RE_F.emplace_back( UseJack, Nergs) ;
      IM_F.emplace_back( UseJack, Nergs) ;
    }

    for(int ierg=0; ierg< Nergs;ierg++) {

      #pragma omp parallel for schedule(dynamic)
      for(int is=0; is<Nsigmas; is++) {

	double syst_RE, syst_IM;
	double l_re, l_im;

	distr_t ID = Get_id_jack_distr(Njacks);

	double E0 = 2*0.135*a;

	double sigma = sigmas[is]*a;

	double Erg= a*MK*sqrt( pow(xks[ierg],2) + pow(xg/2.0,2) );

	double mult_t_RE= 1e-4; double mult_t_IM= 1e-4;

	double kz= xg*a*MK/2;

	cout<<"kz: "<<kz<<endl;


	if(is == 0) {
	  if(ierg==0) {mult_t_RE=0.004; mult_t_IM=0.1;}
	  else if(ierg==1) {mult_t_RE=0.004; mult_t_IM=0.1;}
	  else if(ierg==2) {mult_t_RE=0.005; mult_t_IM=0.1;}
	  else if(ierg==3) {mult_t_RE=0.005; mult_t_IM=0.04;}
	  else if(ierg==4) {mult_t_RE=0.006; mult_t_IM=0.05;}
	  else if(ierg==5) {mult_t_RE=0.005; mult_t_IM=0.03;}
	  else if(ierg==6) {mult_t_RE=0.004; mult_t_IM=0.03;}
	  else if(ierg==7) {mult_t_RE=0.003; mult_t_IM=0.02;}
	  else if(ierg==8) {mult_t_RE=0.003; mult_t_IM=0.01;}
	  else if(ierg==9) {mult_t_RE=0.002; mult_t_IM=0.01;}
	  else crash("ierg: "+to_string(ierg)+" not yet implemented");
	}
	else if( is == 1) {
	  if(ierg==0) {mult_t_RE=0.02; mult_t_IM=0.1;}
	  else if(ierg==1) {mult_t_RE=0.03; mult_t_IM=0.1;}
	  else if(ierg==2) {mult_t_RE=0.02; mult_t_IM=0.2;}
	  else if(ierg==3) {mult_t_RE=0.02; mult_t_IM=0.1;}
	  else if(ierg==4) {mult_t_RE=0.02; mult_t_IM=0.2;}
	  else if(ierg==5) {mult_t_RE=0.01; mult_t_IM=0.1;}
	  else if(ierg==6) {mult_t_RE=0.01; mult_t_IM=0.1;}
	  else if(ierg==7) {mult_t_RE=0.01; mult_t_IM=0.06;}
	  else if(ierg==8) {mult_t_RE=0.01; mult_t_IM=0.06;}
	  else if(ierg==9) {mult_t_RE=0.01; mult_t_IM=0.06;}
	  else crash("ierg: "+to_string(ierg)+" not yet implemented");
	}
	else if( is == 2) {
	  if(ierg==0) {mult_t_RE=0.03; mult_t_IM=0.1;}
	  else if(ierg==1) {mult_t_RE=0.03; mult_t_IM=0.2;}
	  else if(ierg==2) {mult_t_RE=0.03; mult_t_IM=0.2;}
	  else if(ierg==3) {mult_t_RE=0.035; mult_t_IM=0.2;}
	  else if(ierg==4) {mult_t_RE=0.05; mult_t_IM=0.3;}
	  else if(ierg==5) {mult_t_RE=0.03; mult_t_IM=0.3;}
	  else if(ierg==6) {mult_t_RE=0.02; mult_t_IM=0.3;}
	  else if(ierg==7) {mult_t_RE=0.04; mult_t_IM=0.1;}
	  else if(ierg==8) {mult_t_RE=0.035; mult_t_IM=0.15;}
	  else if(ierg==9) {mult_t_RE=0.02; mult_t_IM=0.1;}
	  else crash("ierg: "+to_string(ierg)+" not yet implemented");
	}

	else if( is == 3) {
	  if(ierg==0) {mult_t_RE=0.03; mult_t_IM=0.1;}
	  else if(ierg==1) {mult_t_RE=0.03; mult_t_IM=0.2;}
	  else if(ierg==2) {mult_t_RE=0.03; mult_t_IM=0.2;}
	  else if(ierg==3) {mult_t_RE=0.035; mult_t_IM=0.2;}
	  else if(ierg==4) {mult_t_RE=0.05; mult_t_IM=0.3;}
	  else if(ierg==5) {mult_t_RE=0.03; mult_t_IM=0.3;}
	  else if(ierg==6) {mult_t_RE=0.02; mult_t_IM=0.3;}
	  else if(ierg==7) {mult_t_RE=0.04; mult_t_IM=0.1;}
	  else if(ierg==8) {mult_t_RE=0.035; mult_t_IM=0.15;}
	  else if(ierg==9) {mult_t_RE=0.03; mult_t_IM=0.25;}
	  else crash("ierg: "+to_string(ierg)+" not yet implemented");
	}

	else crash("isigma: "+to_string(is)+" not yet implemented!");

	//last four arguments are IntUpToEmax?, Emax, alpha, UseBW?
	
	RE_F[is].distr_list[ierg] = Get_Laplace_transfo_tmin( Erg,  sigma, E0,  128, 1,  30, prec, "Kl4_RE_xg_0.5", K_RE, C, syst_RE, mult_t_RE ,  l_re , "TANT", "E0_2PI", Ens , 1e-3,0, ID, ID*0.0,  "Kl4", cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.00,1)/kz;
	
	IM_F[is].distr_list[ierg] = Get_Laplace_transfo_tmin( Erg,  sigma, E0,  128, 1,  30, prec, "Kl4_IM_xg_0.5", K_IM, C, syst_IM, mult_t_IM ,  l_im , "TANT", "E0_2PI", Ens , 1e-3,0, ID, ID*0.0,  "Kl4", cov_T, fake_func,0, fake_func_d ,  1 , 4.0, 0.00,1)/kz;
	
      }
      
    }



    


    //Print as a function of  energy for various sigmas

 

    for(int is=0;is<Nsigmas;is++) {
      Print_To_File({} , { xks, RE_F[is].ave(), RE_F[is].err(), IM_F[is].ave(), IM_F[is].err()} , "../data/Kl4_HLT/"+Ens+"/xg_0.5/sigma_"+to_string_with_precision(sigmas[is],3)+".dat", "", "# xk  RE  IM");
    }




  }



  
  return;

}
