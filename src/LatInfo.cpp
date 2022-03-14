#include "../include/LatInfo.h"


using namespace std;

const string Path_to_covariate_input_pars = "../Nf211_ETM_Ensemble_Info/CovMatrix_Branch";



void LatticeInfo::LatInfo(string S) {

  Za=0.0;
  Za_err=1.0;
  Za_M2=0.0;
  Za_M2_err= 1.0;
  double resc=1.0;


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
      Za = 0.731*resc;
      Za_err = 0.008*resc;
      Za_M2 = 0.703*resc;
      Za_M2_err = 0.002*resc;
      
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
      Za = 0.737*resc;
      Za_err = 0.005*resc;
      Za_M2 = 0.714*resc;
      Za_M2_err = 0.002*resc;
      
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
      Za = 0.762*resc;
      Za_err = 0.004*resc;
      Za_M2 = 0.752*resc;
      Za_M2_err = 0.002*resc;
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
    a= 0.09471;          //real is 0.09471 let us change it by 1.5 sigma
    a_err= 0.00039;
    a_nucleon=0.09295;
    a_nucleon_err=0.00047;
    Za=0.724;
    Za_err=0.009;
    Za_WI = 0.695347;
    Za_WI_err= 0.00216267;
    Zv=0.6960;
    Zv_err=0.0007;
    Zv_WI =0.690716;
    Zv_WI_err = 8.78188e-05;
    Zv_WI_charm_extr= 0.601735;
    Zv_WI_charm_extr_err=0.000634;
    Za_WI_charm_extr = 0.66016;
    Za_WI_charm_extr_err = 0.0014;
    if(Tag == "cA211a.53.24") {
      Za_WI_strange= 0.728838 ;
      Za_WI_strange_err = 0.00160579;
      Zv_WI_strange = 0.68670;
      Zv_WI_strange_err = 0.00014;
    }
    else if(Tag == "cA211a.40.24") {
      Za_WI_strange= 0.728838 ;
      Za_WI_strange_err = 0.00160579; 
      Zv_WI_strange = 0.68670;
      Zv_WI_strange_err = 0.00014;
    }
    else if(Tag == "cA211ab.30.32") {
       Za_WI_strange= 0.728838 ;
       Za_WI_strange_err = 0.00160579;
       Zv_WI_strange = 0.68670;
       Zv_WI_strange_err = 0.00014;
    }
    
    ms_L=0.02050; ms_M=0.02300;
    mc_L= 0.26500; mc_M=0.29000; mc_H=0.30000;
  }

  else if(Tag.substr(1,1)=="B") {
    a= 0.08161;
    a_err= 0.00015; //0.00030;   RIMETTI CORRETTO ERRORE
    a_nucleon=0.07975;
    a_nucleon_err=0.00032;
    Za=0.746;
    Za_err=0.008;
    if(Tag == "cB211b.072.64") {
      Za_WI = 0.71403 ;
      Za_WI_err = 0.00077 ;
      Zv_WI = 0.709932 ;
      Zv_WI_err = 0.000007;

      //
      Za_WI_strange= 0.74298;
      Za_WI_strange_err =0.00018;
      Zv_WI_strange = 0.70619;
      Zv_WI_strange_err = 3e-05;
    }
    else if (Tag == "cB211b.072.96") {
      Za_WI = 0.71577 ;
      Za_WI_err = 0.00035 ;
      Zv_WI = 0.709950 ;
      Zv_WI_err = 0.000005;

      //
      Za_WI_strange= 0.742803;
      Za_WI_strange_err = 0.00022 ;
      Zv_WI_strange = 0.706211;
      Zv_WI_strange_err = 3e-05;
    }
    Zv=0.7131;
    Zv_err=0.0006;
    ms_L=0.019; ms_M=0.021;
    mc_L= 0.21000; mc_M=0.23000; mc_H=0.25000;
  }

  else if(Tag.substr(1,1)=="C") {
    a= 0.06942; 
    a_err= 0.00013; // 0.00026; RIMETTI CORRETTO ERRORE
    a_nucleon=0.06860;
    a_nucleon_err=0.00020;
    Za=0.761;
    Za_err=0.008;
    Za_WI = 0.73803 ;
    Za_WI_err = 0.00047 ;
    Zv=0.7310;
    Zv_err=0.0005;
    Zv_WI = 0.728477 ;
    Zv_WI_err = 0.000005;

    //
    Za_WI_strange= 0.7584361736;
    Za_WI_strange_err =0.00020;
    Zv_WI_strange = 0.725291;
    Zv_WI_strange_err = 2.5e-05;
    
    ms_L= 0.01600; ms_M= 0.01800;
    mc_L=0.17500; mc_M=0.19500; mc_H=0.21500;
  }

  else if(Tag.substr(1,1)=="D") {
    a= 0.0577;
    a_err=  0.0001; //0.0002; RIMETTI CORRETTO ERRORE 
    a_nucleon= 0.05777; //fake
    a_nucleon_err= 0.0002; //fake
    Za= 0.76217;
    Za_err=0.00024;
    Za_WI = 0.76134 ;
    Za_WI_err = 0.00027 ;
    Zv_WI = 0.746595 ;
    Zv_WI_err = 0.000005;
    Zv= 0.746599;
    Zv_err = 0.000005;

    //
    Za_WI_strange= 0.77395;
    Za_WI_strange_err = 0.00011;
    Zv_WI_strange = 0.744037;
    Zv_WI_strange_err = 2e-05;


    
    ms_L= 0.014; ms_M= 0.015;
    mc_L=0.165; mc_M=0.175; mc_H=0.175; //mc_H is fake 



  }

  else crash("In LatticeInfo::LatInfo_new_ens Ensemble: "+Tag+" not found while reading lattice spacing");

  if(Tag=="cA211a.53.24") {
    L=24; T=48; ml=0.00530; 
  }
  else if(Tag=="cA211a.40.24") {
    L=24; T=48; ml=0.00400; 
  }
  else if(Tag=="cA211ab.30.32") {
    L=32; T=64; ml=0.00300;
  }
  else if(Tag=="cA211a.12.48") {
    L=48; T=96; ml=0.00120;
  }
  else if(Tag=="cB211a.25.24") {
    L=24; T=48; ml=0.00250;
  }
  else if(Tag=="cB211a.25.32") {
    L=32; T=64; ml=0.00250;
  }
  else if(Tag=="cB211a.25.48") {
    L=48; T=96; ml=0.00250;
  }
  else if(Tag=="cB211a.14.64") {
    L=64; T=128; ml=0.00140; 
  }
  else if(Tag=="cB211b.072.64") {
    L=64; T=128; ml=0.00072;
  }
  else if(Tag=="cB211b.072.96") {
    L=96; T=192; ml=0.00072;
  }
  else if(Tag=="cC211a.20.48") {
    L=48; T=96; ml=0.00200;
  }
  else if(Tag=="cC211a.06.80") {
    L=80; T=160; ml=0.00060; 
  }
  else if(Tag=="cD211a.054.96") {
    L=96; T=192; ml=0.00054;
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

