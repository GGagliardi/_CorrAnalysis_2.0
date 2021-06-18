#include "../include/LatInfo.h"


using namespace std;

const string Path_to_covariate_input_pars = "../Nf211_ETM_Ensemble_Info/CovMatrix_Branch";



void LatticeInfo::LatInfo(string S) {

  Za=0.0;
  Za_err=1.0;
  Za_M2=0.0;
  Za_M2_err= 1.0;


  if(S=="A") {
    Beta= 1.90;
    ZP= 0.529;
    ZP_err=0.007;
    ZS=0.747;
    ZS_err=0.012;
    a = 0.0885;
    a_err = 0.0036;
    ainv= 2.224;
    ainv_err= 0.068;
    Zm_fact=1.629;
    Zm_fact_err= 0.041;
    Za_fact=0.859;
    Za_fact_err=0.015;
    Zm_fact_M2= 1.637;
    Zm_fact_M2_err=0.014;
    Za_fact_M2= 0.990;
    Za_fact_M2_err=0.009;
    
    if(CURRENT_TYPE=="LOCAL") {
      Zv = 0.587;
      Zv_err=0.004;
      Zv_M2=0.608;
      Zv_M2_err=0.003;
      Za = 0.731;
      Za_err = 0.008;
      Za_M2 = 0.703;
      Za_M2_err = 0.002;
      
    }
    else if (CURRENT_TYPE=="CONSERVED") {Zv= 1.0; Zv_err=1.0; Zv_M2=1.0; Zv_M2_err=1.0;}
    else crash("CURRENT_TYPE: "+CURRENT_TYPE+" not yet implemented");
  }
  else if(S=="B") {
    Beta= 1.95;
    ZP= 0.509;
    ZP_err=0.004;
    ZS=0.713;
    ZS_err=0.009;
    a = 0.0815;
    a_err = 0.0030;
    ainv= 2.416;
    ainv_err= 0.063;
    Zm_fact=1.514;
    Zm_fact_err= 0.033;
    Za_fact=0.873;
    Za_fact_err=0.013;
    Zm_fact_M2= 1.585;
    Zm_fact_M2_err=0.012;
    Za_fact_M2= 0.980;
    Za_fact_M2_err=0.008;
    if(CURRENT_TYPE=="LOCAL") {
      Zv=0.603;
      Zv_err=0.003;
      Zv_M2 = 0.614;
      Zv_M2_err = 0.002;
      Za = 0.737;
      Za_err = 0.005;
      Za_M2 = 0.714;
      Za_M2_err = 0.002;
      
    }
    else if(CURRENT_TYPE=="CONSERVED") {Zv= 1.0; Zv_err=1.0;  Zv_M2=1.0; Zv_M2_err=1.0;}
    else crash("CURRENT_TYPE: "+CURRENT_TYPE+" not yet implemented");
    
  }
  
  else if(S=="D") {
    Beta= 2.10;
    ZP= 0.516;
    ZP_err=0.02;
    ZS=0.700;
    ZS_err=0.006;
    a = 0.0619;
    a_err = 0.0018;
    ainv= 3.184;
    ainv_err= 0.059;
    Zm_fact=1.459;
    Zm_fact_err= 0.017;
    Za_fact=0.909;
    Za_fact_err=0.006;
    Zm_fact_M2= 1.462;
    Zm_fact_M2_err=0.006;
    Za_fact_M2= 0.958;
    Za_fact_M2_err=0.003;
    if(CURRENT_TYPE=="LOCAL") {
      Zv= 0.655 ;
      Zv_err=0.003;
      Zv_M2 = 0.657;
      Zv_M2_err = 0.002;
      Za = 0.762;
      Za_err = 0.004;
      Za_M2 = 0.752;
      Za_M2_err = 0.002;
    }
    else if(CURRENT_TYPE=="CONSERVED") {Zv=1.0; Zv_err=1.0;  Zv_M2=1.0; Zv_M2_err=1.0;}
    else crash("CURRENT_TYPE: "+CURRENT_TYPE+" not yet implemented");
  }

   else if(S=="a" || S=="b") {
    Beta= 2.10;
    ZP= 0.733207;
    ZP_err=0.00152823;
    //ZP=1.0;
    //ZP_err= 0.0000001;
    ZS=1.0;
    ZS_err=0.00000;
    a = 0.0619;
    a_err = 0.0018;
    if(CURRENT_TYPE=="LOCAL") {
      Zv= 0.655 ;
      Zv_err=0.003;
      Zv_M2 = 0.657;
      Zv_M2_err = 0.002;
      Za = 0.762;
      Za_err = 0.004;
      Za_M2 = 0.752;
      Za_M2_err = 0.002;
    }
    else if(CURRENT_TYPE=="CONSERVED") {Zv=1.0; Zv_err=1.0;  Zv_M2=1.0; Zv_M2_err=1.0;}
    else crash("CURRENT_TYPE: "+CURRENT_TYPE+" not yet implemented");
  }

  
  else crash("Simulation point: "+S+" is not defined");

  return;


}


void LatticeInfo::LatInfo_new_ens(string Tag) {

  if(Tag.substr(1,1)=="A") {
    a= 0.09471;
    a_err= 0.00039;
    a_nucleon=0.09295;
    a_nucleon_err=0.00047;
  }

  else if(Tag.substr(1,1)=="B") {
    a= 0.08161;
    a_err= 0.00030;
    a_nucleon=0.07975;
    a_nucleon_err=0.00032;
  }

  else if(Tag.substr(1,1)=="C") {
    a= 0.06942;
    a_err= 0.00026;
    a_nucleon=0.06860;
    a_nucleon_err=0.00020;
  }

  else crash("In LatticeInfo::LatInfo_new_ens Ensemble: "+Tag+" not found");

  if(Tag=="cA211.53.24") {
    L=24; T=48; ml=0.00530;
  }
  else if(Tag=="cA211.53.24") {
    L=24; T=48; ml=0.00400;
  }
  else if(Tag=="cA211.53.24") {
    L=32; T=64; ml=0.00300;
  }
  else if(Tag=="cA211.53.24") {
    L=48; T=96; ml=0.00120;
  }
  else if(Tag=="cA211.53.24") {
    L=32; T=64; ml=0.00250;
  }
  else if(Tag=="cA211.53.24") {
    L=48; T=96; ml=0.00250;
  }
  else if(Tag=="cA211.53.24") {
    L=64; T=128; ml=0.00140;
  }
  else if(Tag=="cA211.53.24") {
    L=64; T=128; ml=0.00072;
  }
  else if(Tag=="cA211.53.24") {
    L=48; T=96; ml=0.00200;
  }
  else if(Tag=="cA211.53.24") {
    L=80; T=160; ml=0.00060;
  }
  else crash("In LatticeInfo::LatInfo_new_ens Ensemble: "+Tag+" not found");
  
   



  return;
}




void Read_pars_from_ensemble_tag(vector<string>& Ens_vec, Vfloat& m_lat, Vfloat& L, Vfloat& T) {


  for (auto &Ens: Ens_vec) {

    m_lat.push_back(0.0001*stod(Ens.substr(1, Ens.find(".") -1)));
    L.push_back((double)stoi(Ens.substr(Ens.find(".")+1,Ens.find("_") -1 - Ens.find("."))));
    T.push_back((double)stoi(Ens.substr(Ens.find("_")+1, string::npos)));
  }

  return;
}
void Read_pars_from_ensemble_tag(string Ens, double &m_lat, double &L, double &T) {


   m_lat =0.0001*stod(Ens.substr(1, Ens.find(".") -1));
   L =  (double)stoi(Ens.substr(Ens.find(".")+1,Ens.find("_") -1 - Ens.find(".")));
   T = (double)stoi(Ens.substr(Ens.find("_")+1, string::npos));


   return;


}


void ReadBranch(int k, Eigen::MatrixXd& CovMatrixInput, Eigen::VectorXd& Ave_input_parameters ) {

  //read averages from file
  CovMatrixInput.resize(9,9);
  Ave_input_parameters.resize(9);
  ifstream ReadFromFile(Path_to_covariate_input_pars+to_string(k)+".txt");
  double m_l, B0, f0,ainv_temp1, ainv_temp2, ainv_temp3, Zp_temp1, Zp_temp2, Zp_temp3;
  ReadFromFile>>m_l>>B0>>f0>>ainv_temp1>>ainv_temp2>>ainv_temp3>>Zp_temp1>>Zp_temp2>>Zp_temp3;
  Ave_input_parameters<<m_l,B0,f0,ainv_temp1,ainv_temp2,ainv_temp3,Zp_temp1,Zp_temp2,Zp_temp3;
 


  string CommentLine="";
  while(CommentLine=="") getline(ReadFromFile, CommentLine);
  string Expected= "# Cov Matrix";
  if(CommentLine != Expected) crash("Error while reading CovMatrix, branch nr. "+to_string(k)+"\n Expected: '# Cov Matrix' \t Obtained: "+CommentLine);
  //read Covariance Matrix
  for(int i=0; i<9;i++) {
    for(int j=0;j<9;j++) {
      double temp;
      ReadFromFile >> temp;
      CovMatrixInput(i,j) = temp;
    }
  }

  double stop;
  ReadFromFile >> stop;
  if(!ReadFromFile.eof()) crash("Error: after reading branch nr. "+to_string(k)+" File is !eof");


  ReadFromFile.close();
  cout<<"LOADED COVARIANCE MATRIX FOR BRANCH nr. "<<k<<endl;
  cout<<Ave_input_parameters<<endl;
  cout<<CovMatrixInput<<endl;
 
  
  return;
}
