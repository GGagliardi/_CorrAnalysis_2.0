#include "../include/inclusive_Ds_CKM_analysis.h"
#include "numerics.h"
#include "random.h"


using namespace std;

const bool UseJack = 0;
const int Nboots=1000;





void Get_inclusive_Ds_CKM_analysis() {

  //###########    EXP  INPUT #################

  double Gamma_exp= 8.31*1e-14;  //comb BES CLEO  [ 10^-14 GeV]
  double Gamma_exp_err = 0.20*1e-14; //comb BES CLEO [ 10^-14 GeV]

  double M1_exp = 0.446;  //BES [GeV]
  double M1_exp_err = 0.007 ; //BES [GeV]

  double M2_exp = 0.2245; //BES [GeV^2]
  double M2_exp_err = 0.0046; //BES [GeV^2]

  //###########################################




  //########## LATT INPUT #####################

  //double Vcs_den=0.975; double Vcs_den_err =0.006;
  //double Vcd_den=0.221; double Vcd_den_err = 0.004;

  //double GCS= 8.59; double GCS_err= 0.56;
  //double GCD =12.69; double GCD_err=0.93;

  //double M1CS=0.443; double M1CS_err=0.023;
  //double M1CD=0.714; double M1CD_err=0.057;

  //double M2CS=0.217; double M2CS_err=0.010;
  //double M2CD=0.406; double M2CD_err=0.040;


  //###########################################



  //distr_t Gexp_D(UseJack), M1exp_D(UseJack), M2exp_D(UseJack);

  //distr_t G_cs_D(UseJack), G_cd_D(UseJack);
  //distr_t M1_cs_D(UseJack), M1_cd_D(UseJack);
  //distr_t M2_cs_D(UseJack), M2_cd_D(UseJack);

  //distr_t Vcs_den_D(UseJack), Vcd_den_D(UseJack);


  distr_t Vcs2_den_D(UseJack, Read_From_File("../data/inclusive_CKM_analysis/boot.dat",0,18));
  distr_t Vcd2_den_D(UseJack, Read_From_File("../data/inclusive_CKM_analysis/boot.dat",1,18));
  distr_t G_cs_D(UseJack, Read_From_File("../data/inclusive_CKM_analysis/boot.dat",2,18));
  distr_t G_cd_D(UseJack, Read_From_File("../data/inclusive_CKM_analysis/boot.dat",3,18));
  distr_t M1_hat_cs_D(UseJack, Read_From_File("../data/inclusive_CKM_analysis/boot.dat",5,18));
  distr_t M1_hat_cd_D(UseJack, Read_From_File("../data/inclusive_CKM_analysis/boot.dat",6,18));
  distr_t M2_hat_cs_D(UseJack, Read_From_File("../data/inclusive_CKM_analysis/boot.dat",10,18));
  distr_t M2_hat_cd_D(UseJack, Read_From_File("../data/inclusive_CKM_analysis/boot.dat",11,18));
  distr_t Gamma(UseJack, Read_From_File("../data/inclusive_CKM_analysis/boot.dat",4,18));
  distr_t Gexp_D(UseJack, Read_From_File("../data/inclusive_CKM_analysis/boot.dat",15,18));
  distr_t M1exp_D(UseJack, Read_From_File("../data/inclusive_CKM_analysis/boot.dat",16,18));
  distr_t M2exp_D(UseJack, Read_From_File("../data/inclusive_CKM_analysis/boot.dat",17,18));

    

  for(int i=0;i<Nboots;i++) {
    cout<<"iboot: "<<i<<" "<<M1_hat_cs_D.distr[i]<<" "<<M1_hat_cd_D.distr[i]<<endl;
  }

  distr_t M1_cs_D= M1_hat_cs_D/Gamma;
  distr_t M1_cd_D= M1_hat_cd_D/Gamma;

  distr_t M2_cs_D= M2_hat_cs_D/Gamma;
  distr_t M2_cd_D= M2_hat_cd_D/Gamma;

  distr_t Vcs_den_D = SQRT_D(Vcs2_den_D);
  distr_t Vcd_den_D = SQRT_D(Vcd2_den_D);
  

  GaussianMersenne GM(48593);

  for(int ib=0;ib<Nboots;ib++) {

    //Gexp_D.distr.push_back( Gamma_exp + Gamma_exp_err*GM());
    //M1exp_D.distr.push_back( M1_exp + M1_exp_err*GM());
    //M2exp_D.distr.push_back( M2_exp + M2_exp_err*GM());

    //G_cs_D.distr.push_back( GCS + GCS_err*GM());
    //M1_cs_D.distr.push_back( M1CS + M1CS_err*GM());
    //M2_cs_D.distr.push_back( M2CS + M2CS_err*GM());

    //G_cd_D.distr.push_back( GCD + GCD_err*GM());
    //M1_cd_D.distr.push_back( M1CD + M1CD_err*GM());
    //M2_cd_D.distr.push_back( M2CD + M2CD_err*GM());

    //Vcs_den_D.distr.push_back( Vcs_den + Vcs_den_err*GM());
    //Vcd_den_D.distr.push_back( Vcd_den + Vcd_den_err*GM());
    

  }




  //normalize

  //G_cs_D= G_cs_D/(Vcs_den_D*Vcs_den_D);
  //G_cd_D= G_cd_D/(Vcd_den_D*Vcd_den_D);

  //M1_cs_D = M1_cs_D/(Vcs_den_D*Vcs_den_D);
  //M1_cd_D = M1_cd_D/(Vcd_den_D*Vcd_den_D);

  //M2_cs_D = M2_cs_D/(Vcs_den_D*Vcs_den_D);
  //M2_cd_D = M2_cd_D/(Vcd_den_D*Vcd_den_D);



  //get lines


  double dx=0.0001;
  double Vcd_min=0.1; double Vcd_max=0.50;

  int Nsteps= (int)( (Vcd_max-Vcd_min)/dx );

  Vfloat Vcd2_vals;
  distr_t_list Vcs2_vals_G(UseJack);
  distr_t_list Vcs2_vals_M1(UseJack);
  distr_t_list Vcs2_vals_M2(UseJack);
  distr_t_list Vcs2_vals_M1_hat(UseJack);
  distr_t_list Vcs2_vals_M2_hat(UseJack);
  distr_t_list Vcs2_vals_ave(UseJack);

  distr_t_list Vcs2_vals_G_exp(UseJack);
  distr_t_list Vcs2_vals_M1_exp(UseJack);
  distr_t_list Vcs2_vals_M2_exp(UseJack);
  distr_t_list Vcs2_vals_ave_exp(UseJack);

  distr_t Vcs_uni_G(UseJack), Vcs_uni_M1(UseJack), Vcs_uni_M2(UseJack), Vcs_uni_M1_hat(UseJack), Vcs_uni_M2_hat(UseJack);

  Vcs_uni_G =  SQRT_D( (Gexp_D -1.5*Vcd2_den_D*G_cd_D)/(G_cs_D - 0.0*G_cd_D));
  Vcs_uni_M1 = SQRT_D( (M1exp_D - M1_cd_D)/(M1_cs_D - M1_cd_D));
  Vcs_uni_M2 = SQRT_D( (M2exp_D - M2_cd_D)/(M2_cs_D - M2_cd_D));
  Vcs_uni_M1_hat = SQRT_D(   (M1exp_D*Gexp_D -  Vcd2_den_D*M1_hat_cd_D)/(M1_hat_cs_D - 0.0*M1_hat_cd_D));
  Vcs_uni_M2_hat = SQRT_D(   (M2exp_D*Gexp_D -  Vcd2_den_D*M2_hat_cd_D)/(M2_hat_cs_D - 0.0*M2_hat_cd_D));

  double w1_uni= 1.0/pow(Vcs_uni_G.err(),2); double  w2_uni = 1.0/pow( Vcs_uni_M1_hat.err(),2); double w3_uni = 1.0/pow(Vcs_uni_M2_hat.err(),2);
  double sum_uni = w1_uni + w2_uni+w3_uni;
  w1_uni /= sum_uni; w2_uni /= sum_uni; w3_uni /= sum_uni;

  distr_t Vcs_uni = w1_uni*Vcs_uni_G + w2_uni*Vcs_uni_M1_hat + w3_uni*Vcs_uni_M2_hat;


  cout<<"Vcs(uni,G): "<<Vcs_uni_G.ave()<<" "<<Vcs_uni_G.err()<<endl;
  cout<<"Vcs(uni,M1_hat): "<<Vcs_uni_M1_hat.ave()<<" "<<Vcs_uni_M1_hat.err()<<endl;
  cout<<"Vcs(uni,M2_hat): "<<Vcs_uni_M2_hat.ave()<<" "<<Vcs_uni_M2_hat.err()<<endl;
  cout<<"Vcs(uni): "<<Vcs_uni.ave()<<" "<<Vcs_uni.err()<<endl;
  
  
  

  for(int istep=0;istep<Nsteps;istep++) {

    double Vcd_t= Vcd_min + istep*dx;
    double Vcd2_t= Vcd_t*Vcd_t;
    Vcd2_vals.push_back( Vcd_t*Vcd_t);
    Vcs2_vals_G.distr_list.push_back(  ( Gexp_D - Vcd2_t*G_cd_D)/G_cs_D );
    Vcs2_vals_M1.distr_list.push_back( ( M1exp_D - Vcd2_t*M1_cd_D)/M1_cs_D );
    Vcs2_vals_M2.distr_list.push_back( ( M2exp_D - Vcd2_t*M2_cd_D)/M2_cs_D );

    distr_t b1 = M1exp_D - (M1_hat_cs_D/G_cs_D);
    distr_t b2 = M2exp_D - (M2_hat_cs_D/G_cs_D);
    distr_t a1 = (M1_hat_cd_D/G_cs_D) - (M1exp_D*G_cd_D/G_cs_D);
    distr_t a2 = (M2_hat_cd_D/G_cs_D) - (M2exp_D*G_cd_D/G_cs_D);

    Vcs2_vals_M1_hat.distr_list.push_back( Vcd2_t*a1/b1);
    Vcs2_vals_M2_hat.distr_list.push_back( Vcd2_t*a2/b2);

    double w1= 1.0/pow(Vcs2_vals_G.err(istep),2);
    double w2= 1.0/pow(Vcs2_vals_M1.err(istep),2);
    double w3= 1.0/pow(Vcs2_vals_M2.err(istep),2);

    double sum=w1+w2+w3;
    w1 /= sum;   w2 /= sum; w3 /= sum;

    Vcs2_vals_ave.distr_list.push_back( w1*Vcs2_vals_G.distr_list[istep] + w2*Vcs2_vals_M1.distr_list[istep] + w3*Vcs2_vals_M2.distr_list[istep] );


    Vcs2_vals_G_exp.distr_list.push_back(  ( Gexp_D - Vcd2_t*G_cd_D.ave())/G_cs_D.ave());
    Vcs2_vals_M1_exp.distr_list.push_back( ( M1exp_D - Vcd2_t*M1_cd_D.ave())/M1_cs_D.ave() );
    Vcs2_vals_M2_exp.distr_list.push_back( ( M2exp_D - Vcd2_t*M2_cd_D.ave())/M2_cs_D.ave() );


    double w1_exp= 1.0/pow(Vcs2_vals_G_exp.err(istep),2);
    double w2_exp= 1.0/pow(Vcs2_vals_M1_exp.err(istep),2);
    double w3_exp= 1.0/pow(Vcs2_vals_M2_exp.err(istep),2);

    double sum_exp=w1_exp+w2_exp+w3_exp;
    w1_exp /= sum_exp;   w2_exp /= sum_exp; w3_exp /= sum_exp;

    Vcs2_vals_ave_exp.distr_list.push_back( w1_exp*Vcs2_vals_G_exp.distr_list[istep] + w2_exp*Vcs2_vals_M1_exp.distr_list[istep] + w3_exp*Vcs2_vals_M2_exp.distr_list[istep] );

   
    
    
  }


  //solve three-equations

  //  exp1 =    x*a1  + y*a2     ->   x*(a1*b2 - b1*a2) =  exp1*b2 - exp2*a2 
  //  exp2 =    x*b1  + y*b2     ->   y*(a2*b1 - b2*a1) =  exp1*b1 - exp2*a1

  //G-M1

  distr_t  Vcs2_GM1 = ( Gexp_D*M1_cd_D - M1exp_D*G_cd_D)/( G_cs_D*M1_cd_D - M1_cs_D*G_cd_D);
  distr_t  Vcd2_GM1 = ( Gexp_D*M1_cs_D - M1exp_D*G_cs_D)/( G_cd_D*M1_cs_D - M1_cd_D*G_cs_D);

  distr_t  Vcs2_GM2 = ( Gexp_D*M2_cd_D - M2exp_D*G_cd_D)/( G_cs_D*M2_cd_D - M2_cs_D*G_cd_D);
  distr_t  Vcd2_GM2 = ( Gexp_D*M2_cs_D - M2exp_D*G_cs_D)/( G_cd_D*M2_cs_D - M2_cd_D*G_cs_D);

  distr_t Vcd2_M1M2 = ( M1exp_D*M2_cd_D - M2exp_D*M1_cd_D)/(M1_cs_D*M2_cd_D - M2_cs_D*M1_cd_D);
  distr_t Vcs2_M1M2 = ( M1exp_D*M2_cs_D - M2exp_D*M1_cs_D)/(M1_cd_D*M2_cs_D - M2_cd_D*M1_cs_D);


  //print to file

  boost::filesystem::create_directory("../data/inclusive_CKM_analysis");

  Print_To_File({}, { Vcd2_vals, Vcs2_vals_ave.ave(), Vcs2_vals_ave.err(),  Vcs2_vals_G.ave(), Vcs2_vals_G.err(), Vcs2_vals_M1.ave(), Vcs2_vals_M1.err(), Vcs2_vals_M2.ave(), Vcs2_vals_M2.err() }, "../data/inclusive_CKM_analysis/traj", "", "");

  Print_To_File({}, { Vcd2_vals, Vcs2_vals_M1_hat.ave(), Vcs2_vals_M1_hat.err(), Vcs2_vals_M2_hat.ave(), Vcs2_vals_M2_hat.err() }, "../data/inclusive_CKM_analysis/traj_hat", "", "");

  Print_To_File({}, { Vcd2_vals, Vcs2_vals_ave_exp.ave(), Vcs2_vals_ave_exp.err(),  Vcs2_vals_G_exp.ave(), Vcs2_vals_G_exp.err(), Vcs2_vals_M1_exp.ave(), Vcs2_vals_M1_exp.err(), Vcs2_vals_M2_exp.ave(), Vcs2_vals_M2_exp.err() }, "../data/inclusive_CKM_analysis/traj_exp", "", "");

  ofstream Print_CKM("../data/inclusive_CKM_analysis/Vus_Vud");

  Print_CKM<<"1 "<<Vcs2_GM1.ave()<<" "<<Vcs2_GM1.err()<<" "<<SQRT_D(Vcd2_GM1).ave()<<" "<<SQRT_D(Vcd2_GM1).err()<<endl;
  Print_CKM<<"2 "<<SQRT_D(Vcs2_GM2).ave()<<" "<<SQRT_D(Vcs2_GM2).err()<<" "<<SQRT_D(Vcd2_GM2).ave()<<" "<<SQRT_D(Vcd2_GM2).err()<<endl;
  Print_CKM<<"3 "<<SQRT_D(Vcs2_M1M2).ave()<<" "<<SQRT_D(Vcs2_M1M2).err()<<" "<<SQRT_D(Vcd2_M1M2).ave()<<" "<<SQRT_D(Vcd2_M1M2).err()<<endl;

  Print_CKM.close();

  //weighted average of the three contributions

  double w1_cs = 1.0/(pow(Vcs2_GM1.err(),2)); double w2_cs = 1.0/pow(Vcs2_GM2.err(),2); double w3_cs = 1.0/(pow(Vcs2_M1M2.err(),2));
  double w1_cd = 1.0/(pow(Vcd2_GM1.err(),2)); double w2_cd = 1.0/pow(Vcd2_GM2.err(),2); double w3_cd = 1.0/(pow(Vcd2_M1M2.err(),2));
  double sum_cs= w1_cs + w2_cs + w3_cs;
  double sum_cd= w1_cd + w2_cd + w3_cd;
  //normalize
  w1_cs/= sum_cs; w2_cs/=sum_cs; w3_cs/=sum_cs;
  w1_cd/= sum_cd; w2_cd/=sum_cd; w3_cd/=sum_cd; 


  distr_t Vcs2_weighted= w1_cs*Vcs2_GM1 + w2_cs*Vcs2_GM2 + w3_cs*Vcs2_M1M2;
  distr_t Vcd2_weighted= w1_cd*Vcd2_GM1 + w2_cd*Vcd2_GM2 + w3_cd*Vcd2_M1M2;


  cout<<"|Vcs| : "<<SQRT_D(Vcs2_weighted).ave()<<" "<<SQRT_D(Vcs2_weighted).err()<<endl;
  cout<<"|Vcd| : "<<SQRT_D(Vcd2_weighted).ave()<<" "<<SQRT_D(Vcd2_weighted).err()<<endl;


  

  
  
  

  
  
  



  return;
}
