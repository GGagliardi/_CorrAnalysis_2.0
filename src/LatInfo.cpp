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
    a = 0.0885;
    a_err = 0.0036;
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
    a = 0.0815;
    a_err = 0.0030;
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

