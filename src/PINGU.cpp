#include "../include/PINGU.h"
#include "Corr_analysis.h"
#include "numerics.h"
#include <ostream>


using namespace std;

const double alpha = 1.0/137.035999;
const bool UseJack = 1;
const int Njacks_JJ=150;
const double qu = 2.0/3.0;
const double qd = -qu/2.0;
const double qc = qu;
const double qs = qd;
const int Nhits = 80;
const int Nc = 3;
const int Ndirac = 4;
const bool verbose = false;
const double fm_to_inv_Gev = 1.0 / 0.197327;
const double t_O4 = 12.0;
const int int_t_O4 = 12;
const bool Use_V3=true;
const int GLB_SIGN_REV = (Use_V3) ? -1 : 1;
const bool Re_Im = (Use_V3) ? 1 : 0;
const bool r_average_charm = true;
const bool Is_Bs = true;
const string Tg=(Use_V3) ? "V3P5": "V0P5";
const string HS = (Is_Bs) ? "HSS" : "HSL";
const string PA = (Is_Bs) ? "Ds" : "D";
const string DT = (Is_Bs) ? "phi": "K";
const bool Perform_spectral_reco = true;
const bool Use_conserved_current=true;
const bool Skip_HLT = true;
const double mom = 0.244;
const double mhL = 2.81;
const double mhS = 2.91;
const double mk =0.4946;
const double meta=0.6872;

enum { COL1, COL2, DIR1, DIR2 };
enum {
  COL1_Tcc,
  COL2_Tcc,
  DIR1_Tcc,
  DIR2_Tcc,
  COL1_Tbs,
  COL2_Tbs,
  DIR1_Tbs,
  DIR2_Tbs
};

//1 = fermion, 2= antifermion


vector<string> GG({"G0", "G1", "G2", "G3", "G5"});


void FERM4_T_gprod(const FERM4_T& in, FERM4_T& out, const C_MATRIX& g, int str_id) {

  int S= g.size();
  if(S != in.strp[str_id]) crash("FERM4_T_prod called with g of size: "+to_string(S)+" in.strp["+to_string(str_id)+"] = "+to_string(in.strp[str_id]));
  int N=1;
  for(int i=0;i<str_id;i++) N *= in.strp[i];

 
  
  out = in;
#pragma omp parallel for
  for(int pos=0;pos<in.size();pos++) {
    out.T_list[pos] = 0.0*out.T_list[pos]; //sets to 0
    int l= in.coord(pos,str_id);
    for(int k=0;k<S;k++) {     out.T_list[pos] = out.T_list[pos] + in.T_list[ pos + (k-l)*N ]*g[k][l] ;  }
  }
  

  return;
}

void summ_FERM4_T(FERM4_T &in, const FERM4_T& b) { // sums b and in and assigns it to in


  //check if dimensions match

  assert(in.size() == b.size());
  assert(in.strp.size() == b.strp.size()) ;
  for(int i=0;i<(signed)in.strp.size();i++) assert( in.strp[i] == b.strp[i] );

#pragma omp parallel for
  for(int i=0; i< in.size();i++ ) in.T_list[i] = in.T_list[i] + b.T_list[i];


  return;
}


complex_distr_t_list glb_red_FERM4_T( const FERM4_T &A, const FERM4_T &B, const vector<pair<int,int>>& rules ) {

  int N=1;

  vector<int> DIMS(rules.size());
  map<int, int> NT_MAP_A;
  map<int, int> NT_MAP_B;
  //check
  for(int i=0;i<(signed)rules.size();i++) {
    int r0=rules[i].first;
    int r1=rules[i].second;
    int Ni=0;
   

    bool isA_0 = r0<A.dims;
    bool isA_1 = r1<A.dims;

    if(isA_0 && isA_1) { assert( A.strp[r0] == A.strp[r1] ); Ni = A.strp[r0];}
    if(isA_0 &&  !isA_1) { assert( A.strp[r0] == B.strp[r1-A.dims] ); Ni =A.strp[r0];  }
    if(!isA_0 && isA_1) { assert( B.strp[r0-A.dims] == A.strp[r1]); Ni= A.strp[r1]; }
    if(!isA_0 && !isA_1) { assert( B.strp[r0-A.dims] == B.strp[r1-A.dims]); Ni= B.strp[r0-A.dims];}

    assert(Ni != 0);

    if(isA_0) NT_MAP_A.insert( make_pair( r0, i));
    else NT_MAP_B.insert( make_pair(r0-A.dims,i));

    if(isA_1) NT_MAP_A.insert( make_pair( r1, i));
    else NT_MAP_B.insert( make_pair(r1-A.dims,i));
    
    DIMS[i] =Ni;
    N*=Ni;
  }


  if(verbose) {
    //check that every element of A and B belongs to the map
    //print A and B map
    cout<<"PRINTING A-MAP"<<endl;
    for(auto it = NT_MAP_A.cbegin(); it != NT_MAP_A.cend(); ++it)
      {
	cout << it->first << " " << it->second<< "\n";
      }
    cout<<"PRINTING B-MAP"<<endl;
    for(auto it = NT_MAP_B.cbegin(); it != NT_MAP_B.cend(); ++it)
      {
	cout << it->first << " " << it->second<< "\n";
      }
  
  }

  
  int T1= A.t_size();
  int T2= B.t_size();
  int TT= T1*T2;
  int Njacks = A.T_list[0].distr_list[0].size();
  int UseJack = A.T_list[0].UseJack;

  if(verbose) {
    cout<<"TT: "<<TT<<endl;
    cout<<"Njacks: "<<Njacks<<endl;
    cout<<"Usejack?: "<<UseJack<<endl;
  }

  
  assert(Njacks == B.T_list[0].distr_list[0].size());
  assert(UseJack == B.T_list[0].UseJack);

  complex_distr_t_list ret( 0.0*Get_id_distr_list(TT,Njacks,UseJack), 0.0*Get_id_distr_list(TT,Njacks,UseJack))  ;

 


    for(int n=0;n<N;n++) {

      int loc_N=n;

      //get indices in summation
      vector<int> a;
      for(int i=0;i<(signed)rules.size();i++) { a.push_back( loc_N % DIMS[i]) ; loc_N /= DIMS[i]; }

      int posA = 0; int posB=0; int Na=1; int Nb=1;
      for(int d=0;d<A.dims;d++) { posA += a[NT_MAP_A.find(d)->second]*Na ; Na *= A.strp[d];}
      for(int d=0;d<B.dims;d++) { posB += a[NT_MAP_B.find(d)->second]*Nb ; Nb *= B.strp[d];}

      const complex_distr_t_list& M1 = A.T_list[posA];
      const complex_distr_t_list& M2 = B.T_list[posB];

#pragma omp parallel for
      for(int t1=0;t1<T1;t1++) {
	const complex_distr_t &m = M1.distr_list[t1];
	for(int t2=0;t2<T2;t2++) {

	  int t= t1 + T1*t2;

	  const complex_distr_t &m2 = M2.distr_list[t2];
	  complex_distr_t &m_ret = ret.distr_list[t];
	  
          //ret.distr_list[t] = ret.distr_list[t] + m*M2.distr_list[t2];
	
	  for(int iconf=0;iconf<Njacks;iconf++) m_ret.RE.distr[iconf]  += m.RE.distr[iconf]*m2.RE.distr[iconf] ;
	  for(int iconf=0;iconf<Njacks;iconf++) m_ret.IM.distr[iconf]  += m.RE.distr[iconf]*m2.IM.distr[iconf] ;
	  for(int iconf=0;iconf<Njacks;iconf++) m_ret.IM.distr[iconf]  += m.IM.distr[iconf]*m2.RE.distr[iconf] ;
	  for(int iconf=0;iconf<Njacks;iconf++) m_ret.RE.distr[iconf]  -= m.IM.distr[iconf]*m2.IM.distr[iconf] ;
	     
	      
      
	
	}
      }
    }

  return ret;
}



complex_distr_t_list glb_red_FERM4_T( const FERM4_T &A,  const vector<pair<int,int>>& rules ) {

  int N=1;

  vector<int> DIMS;
  map<int, int> NT_MAP_A;
  //check
  for(int i=0;i<(signed)rules.size();i++) {
    int r0=rules[i].first;
    int r1=rules[i].second;
    int Ni=0;
   
    assert( A.strp[r0] == A.strp[r1] ); Ni = A.strp[r0];
   
    NT_MAP_A.insert( make_pair( r0, i));
    NT_MAP_A.insert( make_pair( r1, i));
        
    
    N *= Ni;
    DIMS.push_back(Ni);
  }

 
 
  int TT= A.t_size();
  int Njacks = A.T_list[0].distr_list[0].size();
  int UseJack = A.T_list[0].UseJack;

  if(verbose) {
    cout<<"TT: "<<TT<<endl;
    cout<<"Njacks: "<<Njacks<<endl;
    cout<<"Usejack?: "<<UseJack<<endl;
  }

  
  
  complex_distr_t_list ret( 0.0*Get_id_distr_list(TT,Njacks,UseJack), 0.0*Get_id_distr_list(TT,Njacks,UseJack))  ;



   

  for(int n=0;n<N;n++) {
    
      int loc_N=n;
      
      //get indices in summation
      vector<int> a;
      for(int i=0;i<(signed)rules.size();i++) { a.push_back( loc_N % DIMS[i]) ; loc_N /= DIMS[i]; }
      
      int posA = 0; int Na=1;
      for(int d=0;d<A.dims;d++) { posA += a[NT_MAP_A.find(d)->second]*Na ; Na *= A.strp[d];}

      complex_distr_t_list M1=A.T_list[posA];
      
#pragma omp parallel for
      for(int t=0;t<TT;t++) {
        
	ret.distr_list[t] = ret.distr_list[t] + M1.distr_list[t];
      
	
      }
  }
  
  return ret;
  
  
}







void Get_PINGU() {

  //TEST_PI_Q2();


  omp_set_num_threads(4);

  if(Perform_spectral_reco) {
    penguin_spectral_reco();
    exit(-1);
  }


  
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





  boost::filesystem::create_directory("../data/PINGU");


 

  //init gamma
  NISSA_GAMMA NG;

  //print nissa gamma
  for(int r=0;r<5;r++) {
    cout<<GG[r]<<endl;
    for(int i=0;i<4;i++) {
      for(int j=0;j<4;j++) {
	cout<<NG.G[r][i][j]<<" ";
      }
      cout<<endl;
    }
  }

  //print B1
  C_MATRIX B1 = complex_prod_matr( NG.G[2], NG.G[3]);
  cout<<"B1"<<endl;
  for(int i=0;i<4;i++) {
    for(int j=0;j<4;j++) {
      cout<<B1[i][j]<<" ";
    }
    cout<<endl;
  }

  //read confs list
  vector<string> Confs;
  ifstream read_conf("../PINGU/cB211b.072.64/confs_list_conserved");
  while(! read_conf.eof()) {
    string a;
    read_conf >> a;
    if(!read_conf.eof()) Confs.push_back(a);
    
  }
  read_conf.close();
  int Nconfs=Confs.size();
  cout<<"Nconfs: "<<Nconfs<<endl;

  
  
  //load data
  CorrAnalysis Corr(UseJack, Nconfs,100);
  CorrAnalysis Corr_JJ(UseJack, Njacks_JJ,100);
  Corr.Nt=128; //B64
  Corr_JJ.Nt=128;
  int T=Corr.Nt;
  double V = (double)pow(T/2,3.0);
  double V4 = V*((double)T);
  Corr.Perform_Nt_t_average =0;
  Corr_JJ.Perform_Nt_t_average=0; 
  

 

  
  
    
  vector<vector<double>> NULL_VECTOR(T*T);
  vector<vector<double>> SNULL_VECTOR(T);
  for(int t=0;t<T*T;t++)
    for(int ic=0;ic<Nconfs;ic++) NULL_VECTOR[t].push_back(0.0);
  for(int t=0;t<T;t++)
    for(int ic=0;ic<Nconfs;ic++) SNULL_VECTOR[t].push_back(0.0);

  
  complex_distr_t_list C_O1(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O2(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O1_odd(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O1_VV(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O1_AA(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O1_VA(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O1_AV(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O2_VV(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O2_AA(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O2_VA(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O2_AV(0, NULL_VECTOR, NULL_VECTOR);

  complex_distr_t_list C_O1_not_ave_VV(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O1_not_ave_AA(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O2_not_ave_VV(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O2_not_ave_AA(0, NULL_VECTOR, NULL_VECTOR);
 
  
  complex_distr_t_list C_O2_odd(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O1_even(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O2_even(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O1_VV_FIERZ(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O1_AA_FIERZ(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O2_VV_FIERZ(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_O2_AA_FIERZ(0, NULL_VECTOR, NULL_VECTOR);
  complex_distr_t_list C_semlep_S0(0, SNULL_VECTOR, SNULL_VECTOR);
  complex_distr_t_list C_semlep_V0(0, SNULL_VECTOR, SNULL_VECTOR);
  complex_distr_t_list C_semlep_V3(0, SNULL_VECTOR, SNULL_VECTOR);
  complex_distr_t_list C_semlep_A3(0, SNULL_VECTOR, SNULL_VECTOR);
  complex_distr_t_list C_semlep_A3_ave(0, SNULL_VECTOR, SNULL_VECTOR);
  complex_distr_t_list C_etacc(0, SNULL_VECTOR, SNULL_VECTOR);

  
  complex_distr_t_list C_cc_VKV3(0, SNULL_VECTOR, SNULL_VECTOR);
  complex_distr_t_list C_cc_VKV0(0, SNULL_VECTOR, SNULL_VECTOR);
  complex_distr_t_list C_cc_VKA3(0, SNULL_VECTOR, SNULL_VECTOR);
  complex_distr_t_list C_cc_VKA0(0, SNULL_VECTOR, SNULL_VECTOR);
  



  //loop over hits

 

  for(int ihit=0;ihit<Nhits;ihit++) {

    
      auto start = chrono::system_clock::now();

        
      FERM4_T S_cc({Nc,Nc,Ndirac,Ndirac}, Nc*Nc*Ndirac*Ndirac);
      FERM4_T T_cc_P({Nc,Nc,Ndirac,Ndirac}, Nc*Nc*Ndirac*Ndirac);
      FERM4_T T_cc_VK({Nc,Nc,Ndirac,Ndirac}, Nc*Nc*Ndirac*Ndirac);
      FERM4_T T_cc_VK_REV({Nc,Nc,Ndirac,Ndirac}, Nc*Nc*Ndirac*Ndirac);
      FERM4_T S_cc_REV({Nc,Nc,Ndirac,Ndirac}, Nc*Nc*Ndirac*Ndirac);

 
#pragma omp parallel for
      for(int alpha=0;alpha<Ndirac;alpha++)
	for(int beta=0;beta<Ndirac;beta++) 
	  for(int a=0;a<Nc;a++)
	    for(int b=0;b<Nc;b++) {

	      string coldir_str = to_string(beta)+"_"+to_string(b)+"_"+to_string(alpha)+"_"+to_string(a);

              string coldir_str_RP =  to_string(alpha)+"_"+to_string(a)+"_"+to_string(beta)+"_"+to_string(b);

	      int pos = b + a*Nc + beta*Nc*Nc + alpha*Nc*Nc*Ndirac;

	      
	      VVfloat READ_S_RE(T), READ_CC_P5P5_RE(T), READ_CC_VK_RE(T);
	      VVfloat READ_CC_VK_REV_RE(T),  READ_S_REV_RE(T);
	      VVfloat READ_CC_VK_RP_RE(T), READ_CC_VK_RP_REV_RE(T);
	      
	      VVfloat READ_S_IM(T), READ_CC_P5P5_IM(T), READ_CC_VK_IM(T);
	      VVfloat READ_CC_VK_REV_IM(T), READ_S_REV_IM(T);
	      VVfloat READ_CC_VK_RP_IM(T), READ_CC_VK_RP_REV_IM(T);

	      //fast read
	      for(int iconf=0;iconf<Nconfs;iconf++) {
		string conf= Confs[iconf];
		string File_T = "../PINGU/cB211b.072.64/"+conf+"/TT_"+coldir_str+"."+to_string(ihit+1);
		string File_T_RP = "../PINGU/cB211b.072.64/"+conf+"/TT_"+coldir_str_RP+"."+to_string(ihit+1);
		string File_T_CONS = "../PINGU_CONSERVED/cB211b.072.64/"+conf+"/TT_CONS_"+coldir_str+"."+to_string(ihit+1);
		string File_T_CONS_RP = "../PINGU_CONSERVED/cB211b.072.64/"+conf+"/TT_CONS_"+coldir_str_RP+"."+to_string(ihit+1);

		VVfloat T_CONS_RE(4), T_CONS_IM(4), T_RP_CONS_RE(4), T_RP_CONS_IM(4);
		VVfloat T_P5P5_RE(6), T_P5P5_IM(6), T_VKP5_RE(6), T_VKP5_IM(6);
		VVfloat T_RP_P5P5_RE(6), T_RP_P5P5_IM(6), T_RP_VKP5_RE(6), T_RP_VKP5_IM(6); 
		string dummy="";
		string file=File_T;
		ifstream read(file);
		getline(read, dummy);
		getline(read,dummy);
		for(int t=0;t<T;t++) {
		  for(int ic=0;ic<(signed)T_P5P5_RE.size();ic++) {
		    double re, im;
		    read>>re>>im;
		    T_P5P5_RE[ic].push_back(re);
		    T_P5P5_IM[ic].push_back(im);
		  }
		}
		getline(read,dummy);
		getline(read,dummy);
		for(int t=0;t<T;t++) {
		  for(int ic=0;ic<(signed)T_P5P5_RE.size();ic++) {
		    double re, im;
		    read>>re>>im;
		    T_VKP5_RE[ic].push_back(re);
		    T_VKP5_IM[ic].push_back(im);
		  }
		}
		read.close();
		file=File_T_RP;
		read.open(file);
		getline(read, dummy);
		getline(read,dummy);
		for(int t=0;t<T;t++) {
		  for(int ic=0;ic<(signed)T_RP_P5P5_RE.size();ic++) {
		    double re, im;
		    read>>re>>im;
		    T_RP_P5P5_RE[ic].push_back(re);
		    T_RP_P5P5_IM[ic].push_back(im);
		  }
		}
		getline(read,dummy);
		getline(read,dummy);
		for(int t=0;t<T;t++) {
		  for(int ic=0;ic<(signed)T_RP_P5P5_RE.size();ic++) {
		    double re, im;
		    read>>re>>im;
		    T_RP_VKP5_RE[ic].push_back(re);
		    T_RP_VKP5_IM[ic].push_back(im);
		  }
		}
		read.close();

		file=File_T_CONS;
		read.open(file);
		getline(read,dummy);
		getline(read,dummy);
		for(int t=0;t<T;t++) {
		  for(int ic=0;ic<(signed)T_CONS_RE.size();ic++) {
		    double re, im;
		    read>>re>>im;
		    T_CONS_RE[ic].push_back(re);
		    T_CONS_IM[ic].push_back(im);
		  }
		}
		read.close();
		file=File_T_CONS_RP;
		read.open(file);
		getline(read,dummy);
		getline(read,dummy);
		for(int t=0;t<T;t++) {
		  for(int ic=0;ic<(signed)T_RP_CONS_RE.size();ic++) {
		    double re, im;
		    read>>re>>im;
		    T_RP_CONS_RE[ic].push_back(re);
		    T_RP_CONS_IM[ic].push_back(im);
		  }
		}
		
	
		//push_back
		//CC
		for(int t=0;t< T;t++) {
		  READ_CC_P5P5_RE[t].push_back( T_P5P5_RE[4][t]);
		  READ_CC_P5P5_IM[t].push_back( T_P5P5_IM[4][t]);
		  READ_CC_VK_RE[t].push_back( (Use_conserved_current)?T_CONS_RE[2][t]:T_VKP5_RE[4][t]);  //cons or RE
		  READ_CC_VK_IM[t].push_back( (Use_conserved_current)?T_CONS_IM[2][t]:T_VKP5_IM[4][t]);  //cons or RE
		  //CC REV
		  READ_CC_VK_REV_RE[t].push_back( (Use_conserved_current)?T_CONS_RE[3][t]:T_VKP5_RE[5][t]); //cons or RE
		  READ_CC_VK_REV_IM[t].push_back( (Use_conserved_current)?T_CONS_IM[3][t]:T_VKP5_IM[5][t]); //cons or RE
		  //CC RP
		  READ_CC_VK_RP_RE[t].push_back( (Use_conserved_current)?T_RP_CONS_RE[2][t]:T_RP_VKP5_RE[4][t]); //cons or RE
		  READ_CC_VK_RP_IM[t].push_back( (Use_conserved_current)?T_RP_CONS_IM[2][t]:T_RP_VKP5_IM[4][t]); //cons or RE
		  //CC RP REV
		  READ_CC_VK_RP_REV_RE[t].push_back( (Use_conserved_current)?T_RP_CONS_RE[3][t]:T_RP_VKP5_RE[5][t]); //cons or RE
		  READ_CC_VK_RP_REV_IM[t].push_back( (Use_conserved_current)?T_RP_CONS_IM[3][t]:T_RP_VKP5_IM[5][t]); //cons or RE
		  //HH
		  int offset= (Is_Bs)?2:0;
		  READ_S_RE[t].push_back( T_P5P5_RE[0+offset][t]);
		  READ_S_IM[t].push_back( T_P5P5_IM[0+offset][t]);
		  //HH REV
		  READ_S_REV_RE[t].push_back( T_P5P5_RE[1+offset][t]);
		  READ_S_REV_IM[t].push_back( T_P5P5_IM[1+offset][t]);
		  
		}
	      }
	      
	      T_cc_P.T_list[pos] = complex_distr_t_list(0, READ_CC_P5P5_RE, READ_CC_P5P5_IM );
	      T_cc_VK.T_list[pos] = complex_distr_t_list( 0 , summ_master(READ_CC_VK_RE, Multiply_Vvector_by_scalar(READ_CC_VK_RP_REV_RE,-1.0*r_average_charm)), summ_master(READ_CC_VK_IM, Multiply_Vvector_by_scalar(READ_CC_VK_RP_REV_IM, 1.0*r_average_charm)));
	      S_cc.T_list[pos] = complex_distr_t_list( 0,  READ_S_RE , READ_S_IM ) ;
	      T_cc_VK_REV.T_list[pos] = complex_distr_t_list( 0 , summ_master(READ_CC_VK_REV_RE, Multiply_Vvector_by_scalar(READ_CC_VK_RP_RE,-1.0*r_average_charm)), summ_master(READ_CC_VK_REV_IM,Multiply_Vvector_by_scalar(READ_CC_VK_RP_IM, 1.0*r_average_charm) ) );
	      S_cc_REV.T_list[pos] = complex_distr_t_list( 0,  READ_S_REV_RE , READ_S_REV_IM ) ;
	      
	      
	      
	      

	      
	      
	      
	      
	    }

    
      //perform all matrix multiplications
      //forward
      FERM4_T T_cc_GAMMA;
      FERM4_T T_cc_GAMMA_AA;
      FERM4_T T_cc_GAMMA_VV;
      FERM4_T T_cc_GAMMA_AV;
      FERM4_T T_cc_GAMMA_VA;
      //backward
      FERM4_T T_cc_GAMMA_REV;
      FERM4_T T_cc_GAMMA_AA_REV;
      FERM4_T T_cc_GAMMA_VV_REV;
      FERM4_T T_cc_GAMMA_AV_REV;
      FERM4_T T_cc_GAMMA_VA_REV;

      FERM4_T S_cc_S0;
      FERM4_T S_cc_V0;
      FERM4_T S_cc_V3;
      FERM4_T S_cc_A3;

      FERM4_T S_cc_S0_REV;
      FERM4_T S_cc_V0_REV;
      FERM4_T S_cc_V3_REV;
      FERM4_T S_cc_A3_REV;
      
      FERM4_T T_cc_TEST;
      FERM4_T T_cc_VKV3, T_cc_VKA3, T_cc_VKA0, T_cc_VKV0;
     
  
 
      //1= fermion, 2=anti-fermion 
      vector<pair<int,int>> RULE_O1;  
      RULE_O1.push_back( make_pair(COL1_Tcc,COL2_Tcc));
      RULE_O1.push_back( make_pair(DIR1_Tcc,DIR2_Tbs));
      RULE_O1.push_back( make_pair(DIR2_Tcc,DIR1_Tbs));
      RULE_O1.push_back( make_pair(COL1_Tbs,COL2_Tbs));

      vector<pair<int,int>> RULE_O2;
      RULE_O2.push_back( make_pair(COL1_Tcc,COL2_Tbs));
      RULE_O2.push_back( make_pair(COL2_Tcc,COL1_Tbs));
      RULE_O2.push_back( make_pair(DIR1_Tcc,DIR2_Tbs));
      RULE_O2.push_back( make_pair(DIR2_Tcc,DIR1_Tbs));

    
      //contract Fierz
      vector<pair<int,int>> RULE_O1_FIERZ;

      RULE_O1_FIERZ.push_back( make_pair(COL1_Tcc,COL2_Tcc));
      RULE_O1_FIERZ.push_back( make_pair(DIR2_Tcc,DIR1_Tcc));
      RULE_O1_FIERZ.push_back( make_pair(DIR2_Tbs,DIR1_Tbs));
      RULE_O1_FIERZ.push_back( make_pair(COL1_Tbs,COL2_Tbs));

      vector<pair<int,int>> RULE_O2_FIERZ;

      RULE_O2_FIERZ.push_back( make_pair(COL1_Tcc,COL2_Tbs));
      RULE_O2_FIERZ.push_back( make_pair(COL2_Tcc,COL1_Tbs));
      RULE_O2_FIERZ.push_back( make_pair(DIR2_Tcc,DIR1_Tcc));
      RULE_O2_FIERZ.push_back( make_pair(DIR2_Tbs,DIR1_Tbs));

      vector<pair<int,int>> RULE_SINGLE;
      RULE_SINGLE.push_back( make_pair(COL1, COL2) );
      RULE_SINGLE.push_back( make_pair(DIR1, DIR2) );

      complex_distr_t_list RES_O1_VV_FIERZ;
      complex_distr_t_list RES_O1_AA_FIERZ;

      complex_distr_t_list RES_O1_VV_FIERZ_REV;
      complex_distr_t_list RES_O1_AA_FIERZ_REV;

      complex_distr_t_list RES_O2_VV_FIERZ;
      complex_distr_t_list RES_O2_AA_FIERZ;

      complex_distr_t_list RES_O2_VV_FIERZ_REV;
      complex_distr_t_list RES_O2_AA_FIERZ_REV;

     
    


      auto end = chrono::system_clock::now();
    
      chrono::duration<double> elapsed_seconds = end-start;
    
      cout<<"elapsed time before applying V-A to ccbar: "<<elapsed_seconds.count()<<" s"<<flush<<endl;

      
    

      for(int D=0; D<Ndirac;D++) {

      
	FERM4_T out;
	FERM4_T out2;
	FERM4_T out_V;
	FERM4_T out_A;
	FERM4_T out_AA;
	FERM4_T out_VV;
	FERM4_T out_AV;
	FERM4_T out_VA;

	FERM4_T out_REV;
	FERM4_T out2_REV;
	FERM4_T out_V_REV;
	FERM4_T out_A_REV;
	FERM4_T out_AA_REV;
	FERM4_T out_VV_REV;
	FERM4_T out_AV_REV;
	FERM4_T out_VA_REV;

	FERM4_T T_cc_GMU;
	FERM4_T T_ss_GMU;

	FERM4_T T_cc_GAMU;
	FERM4_T T_ss_GAMU;

	FERM4_T T_cc_GMU_REV;
	FERM4_T T_ss_GMU_REV;

	FERM4_T T_cc_GAMU_REV;
	FERM4_T T_ss_GAMU_REV;



      
	//MULTIPLY c quark by V-A 
	//FORWARD
	FERM4_T_gprod(T_cc_VK, out_A, complex_prod_matr( NG.G[4], TRANSPOSE( NG.G5PROD(D) )   ) , DIR1);
	FERM4_T_gprod(T_cc_VK, out_V, complex_prod_matr( NG.G[4], TRANSPOSE(NG.G[D] )   ) , DIR1);
	FERM4_T_gprod(T_cc_VK, out, complex_prod_matr( NG.G[4], TRANSPOSE( complex_diff_matr( NG.G[D] , NG.G5PROD(D) ) )) , DIR1);

	FERM4_T_gprod(out, out2, complex_prod_matr( complex_diff_matr( NG.G[D] , NG.G5PROD(D)), NG.G[4]), DIR2);
	out.free_mem();
	FERM4_T_gprod(out_A, out_AA, complex_prod_matr( NG.G5PROD(D), NG.G[4]), DIR2);
	FERM4_T_gprod(out_A, out_AV, complex_prod_matr( NG.G[D], NG.G[4]), DIR2);
	out_A.free_mem();
	FERM4_T_gprod(out_V, out_VV, complex_prod_matr( NG.G[D], NG.G[4]), DIR2);
	FERM4_T_gprod(out_V, out_VA, complex_prod_matr( NG.G5PROD(D), NG.G[4]), DIR2);
	out_V.free_mem();


	

	//ADDITIONAL CONTRACTION FOR FIERZING
       
	FERM4_T_gprod(T_cc_VK, T_cc_GMU, complex_prod_matr( NG.G[D], NG.G[4]), DIR2);
	FERM4_T_gprod(T_cc_VK, T_cc_GAMU, complex_prod_matr( NG.G5PROD(D), NG.G[4]), DIR2);
	FERM4_T_gprod(S_cc, T_ss_GMU, complex_prod_matr( NG.G[D], NG.G[4]), DIR2);
	FERM4_T_gprod(S_cc, T_ss_GAMU, complex_prod_matr( NG.G5PROD(D), NG.G[4]), DIR2);

	FERM4_T_gprod(T_cc_VK_REV, T_cc_GMU_REV, complex_prod_matr( NG.G[D], NG.G[4]), DIR2);
	FERM4_T_gprod(T_cc_VK_REV, T_cc_GAMU_REV, complex_prod_matr( NG.G5PROD(D), NG.G[4]), DIR2);
	FERM4_T_gprod(S_cc_REV, T_ss_GMU_REV, complex_prod_matr( NG.G[D], NG.G[4]), DIR2);
	FERM4_T_gprod(S_cc_REV, T_ss_GAMU_REV, complex_prod_matr( NG.G5PROD(D), NG.G[4]), DIR2);
	
	if(D==0) {
	  RES_O1_VV_FIERZ = glb_red_FERM4_T( T_cc_GMU, T_ss_GMU_REV, RULE_O1_FIERZ );
	  RES_O1_AA_FIERZ= glb_red_FERM4_T( T_cc_GAMU, T_ss_GAMU_REV, RULE_O1_FIERZ );

	  RES_O1_VV_FIERZ_REV = glb_red_FERM4_T( T_cc_GMU_REV, T_ss_GMU, RULE_O1_FIERZ );
	  RES_O1_AA_FIERZ_REV= glb_red_FERM4_T( T_cc_GAMU_REV, T_ss_GAMU, RULE_O1_FIERZ );

	  RES_O2_VV_FIERZ = glb_red_FERM4_T( T_cc_GMU, T_ss_GMU_REV, RULE_O2_FIERZ );
	  RES_O2_AA_FIERZ= glb_red_FERM4_T( T_cc_GAMU, T_ss_GAMU_REV, RULE_O2_FIERZ );

	  RES_O2_VV_FIERZ_REV = glb_red_FERM4_T( T_cc_GMU_REV, T_ss_GMU, RULE_O2_FIERZ );
	  RES_O2_AA_FIERZ_REV= glb_red_FERM4_T( T_cc_GAMU_REV, T_ss_GAMU, RULE_O2_FIERZ );
	  
	}
	else {
	  RES_O1_VV_FIERZ = RES_O1_VV_FIERZ + glb_red_FERM4_T( T_cc_GMU, T_ss_GMU_REV, RULE_O1_FIERZ );
	  RES_O1_AA_FIERZ= RES_O1_AA_FIERZ+ glb_red_FERM4_T( T_cc_GAMU, T_ss_GAMU_REV, RULE_O1_FIERZ );

	  RES_O1_VV_FIERZ_REV = RES_O1_VV_FIERZ_REV + glb_red_FERM4_T( T_cc_GMU_REV, T_ss_GMU, RULE_O1_FIERZ );
	  RES_O1_AA_FIERZ_REV= RES_O1_AA_FIERZ_REV+ glb_red_FERM4_T( T_cc_GAMU_REV, T_ss_GAMU, RULE_O1_FIERZ );

	  RES_O2_VV_FIERZ = RES_O2_VV_FIERZ + glb_red_FERM4_T( T_cc_GMU, T_ss_GMU_REV, RULE_O2_FIERZ );
	  RES_O2_AA_FIERZ= RES_O2_AA_FIERZ+ glb_red_FERM4_T( T_cc_GAMU, T_ss_GAMU_REV, RULE_O2_FIERZ );

	  RES_O2_VV_FIERZ_REV = RES_O2_VV_FIERZ_REV + glb_red_FERM4_T( T_cc_GMU_REV, T_ss_GMU, RULE_O2_FIERZ );
	  RES_O2_AA_FIERZ_REV= RES_O2_AA_FIERZ_REV+ glb_red_FERM4_T( T_cc_GAMU_REV, T_ss_GAMU, RULE_O2_FIERZ );
	}
	
      
	if(D==0) { T_cc_GAMMA = out2; T_cc_GAMMA_AA= out_AA; T_cc_GAMMA_AV= out_AV; T_cc_GAMMA_VA = out_VA; T_cc_GAMMA_VV = out_VV; }
	else  { summ_FERM4_T(T_cc_GAMMA, out2); summ_FERM4_T(T_cc_GAMMA_AA, out_AA); summ_FERM4_T(T_cc_GAMMA_AV, out_AV); summ_FERM4_T(T_cc_GAMMA_VA, out_VA); summ_FERM4_T(T_cc_GAMMA_VV, out_VV); }


	T_cc_GMU.free_mem();
	T_ss_GMU.free_mem();
	T_cc_GAMU.free_mem();
	T_ss_GAMU.free_mem();
	T_cc_GMU_REV.free_mem();
	T_ss_GMU_REV.free_mem();
	T_cc_GAMU_REV.free_mem();
	T_ss_GAMU_REV.free_mem();
	out2.free_mem();
	out_AA.free_mem();
	out_AV.free_mem();
	out_VA.free_mem();
	out_VV.free_mem();

	//BACKWARD
	FERM4_T_gprod(T_cc_VK_REV, out_A_REV, complex_prod_matr( NG.G[4], TRANSPOSE( NG.G5PROD(D) )   ) , DIR1);
	FERM4_T_gprod(T_cc_VK_REV, out_V_REV, complex_prod_matr( NG.G[4], TRANSPOSE(NG.G[D] )   ) , DIR1);
	FERM4_T_gprod(T_cc_VK_REV, out_REV, complex_prod_matr( NG.G[4], TRANSPOSE( complex_diff_matr( NG.G[D] , NG.G5PROD(D) ) )) , DIR1);
	FERM4_T_gprod(out_REV, out2_REV, complex_prod_matr( complex_diff_matr( NG.G[D] , NG.G5PROD(D)), NG.G[4]), DIR2);
	out_REV.free_mem();
	FERM4_T_gprod(out_A_REV, out_AA_REV, complex_prod_matr( NG.G5PROD(D), NG.G[4]), DIR2);
	FERM4_T_gprod(out_A_REV, out_AV_REV, complex_prod_matr( NG.G[D], NG.G[4]), DIR2);
	out_A_REV.free_mem();
	FERM4_T_gprod(out_V_REV, out_VV_REV, complex_prod_matr( NG.G[D], NG.G[4]), DIR2);
	FERM4_T_gprod(out_V_REV, out_VA_REV, complex_prod_matr( NG.G5PROD(D), NG.G[4]), DIR2);
	out_V_REV.free_mem();
	//BACKWARD
	if(D==0) { T_cc_GAMMA_REV = out2_REV; T_cc_GAMMA_AA_REV= out_AA_REV; T_cc_GAMMA_AV_REV= out_AV_REV; T_cc_GAMMA_VA_REV = out_VA_REV; T_cc_GAMMA_VV_REV = out_VV_REV; }
	else  { summ_FERM4_T(T_cc_GAMMA_REV, out2_REV); summ_FERM4_T(T_cc_GAMMA_AA_REV, out_AA_REV); summ_FERM4_T(T_cc_GAMMA_AV_REV, out_AV_REV); summ_FERM4_T(T_cc_GAMMA_VA_REV, out_VA_REV); summ_FERM4_T(T_cc_GAMMA_VV_REV, out_VV_REV); }

      }

      FERM4_T_gprod(S_cc, S_cc_S0, NG.G[4], DIR1);
      FERM4_T_gprod(S_cc, S_cc_V0, complex_prod_matr(NG.G[0],NG.G[4]), DIR2);
      FERM4_T_gprod(S_cc, S_cc_V3, complex_prod_matr(NG.G[3],NG.G[4]), DIR2);
      FERM4_T_gprod(S_cc, S_cc_A3, NG.G[3], DIR2);

      FERM4_T_gprod(S_cc_REV, S_cc_S0_REV, NG.G[4], DIR1);
      FERM4_T_gprod(S_cc_REV, S_cc_V0_REV, complex_prod_matr(NG.G[0],NG.G[4]), DIR2);
      FERM4_T_gprod(S_cc_REV, S_cc_V3_REV, complex_prod_matr(NG.G[3],NG.G[4]), DIR2);
      FERM4_T_gprod(S_cc_REV, S_cc_A3_REV, NG.G[3], DIR2);

      
      FERM4_T_gprod(T_cc_VK, T_cc_VKV3, NG.G5PROD(3), DIR1);
      FERM4_T_gprod(T_cc_VK, T_cc_VKV0, NG.G5PROD(0), DIR1);
      FERM4_T_gprod(T_cc_VK, T_cc_VKA3, NG.G[3], DIR1);
      FERM4_T_gprod(T_cc_VK, T_cc_VKA0, NG.G[0], DIR1);

      T_cc_VK.free_mem();
           

      auto end2 = chrono::system_clock::now();
    
      chrono::duration<double> elapsed_seconds_bC = end2-start;
    
      cout<<"elapsed time before tensor contraction: "<<elapsed_seconds_bC.count()<<" s"<<flush<<endl;
  

    
      //FORWARD
      complex_distr_t_list RES_O1 = glb_red_FERM4_T( T_cc_GAMMA, S_cc_REV, RULE_O1 );
      complex_distr_t_list RES_O2 = glb_red_FERM4_T( T_cc_GAMMA, S_cc_REV, RULE_O2 );
      complex_distr_t_list RES_O1_REV = glb_red_FERM4_T( T_cc_GAMMA_REV, S_cc, RULE_O1 );
      complex_distr_t_list RES_O2_REV = glb_red_FERM4_T( T_cc_GAMMA_REV, S_cc, RULE_O2 );
     
     
      //AA

      complex_distr_t_list RES_O1_AA = glb_red_FERM4_T( T_cc_GAMMA_AA, S_cc_REV, RULE_O1 );
      complex_distr_t_list RES_O2_AA = glb_red_FERM4_T( T_cc_GAMMA_AA, S_cc_REV, RULE_O2 );
      //VV
      complex_distr_t_list RES_O1_VV = glb_red_FERM4_T( T_cc_GAMMA_VV, S_cc_REV, RULE_O1 );
      complex_distr_t_list RES_O2_VV = glb_red_FERM4_T( T_cc_GAMMA_VV, S_cc_REV, RULE_O2 );
      //AV
      complex_distr_t_list RES_O1_AV = glb_red_FERM4_T( T_cc_GAMMA_AV, S_cc_REV, RULE_O1 );
      complex_distr_t_list RES_O2_AV = glb_red_FERM4_T( T_cc_GAMMA_AV, S_cc_REV, RULE_O2 );
      //VA
      complex_distr_t_list RES_O1_VA = glb_red_FERM4_T( T_cc_GAMMA_VA, S_cc_REV, RULE_O1 );
      complex_distr_t_list RES_O2_VA = glb_red_FERM4_T( T_cc_GAMMA_VA, S_cc_REV, RULE_O2 );

      //BACKWARD
      
      //AA
      complex_distr_t_list RES_O1_AA_REV = glb_red_FERM4_T( T_cc_GAMMA_AA_REV, S_cc, RULE_O1 );
      complex_distr_t_list RES_O2_AA_REV = glb_red_FERM4_T( T_cc_GAMMA_AA_REV, S_cc, RULE_O2 );
      //VV
      complex_distr_t_list RES_O1_VV_REV = glb_red_FERM4_T( T_cc_GAMMA_VV_REV, S_cc, RULE_O1 );
      complex_distr_t_list RES_O2_VV_REV = glb_red_FERM4_T( T_cc_GAMMA_VV_REV, S_cc, RULE_O2 );
      //AV
      complex_distr_t_list RES_O1_AV_REV = glb_red_FERM4_T( T_cc_GAMMA_AV_REV, S_cc, RULE_O1 );
      complex_distr_t_list RES_O2_AV_REV = glb_red_FERM4_T( T_cc_GAMMA_AV_REV, S_cc, RULE_O2 );
      //VA
      complex_distr_t_list RES_O1_VA_REV = glb_red_FERM4_T( T_cc_GAMMA_VA_REV, S_cc, RULE_O1 );
      complex_distr_t_list RES_O2_VA_REV = glb_red_FERM4_T( T_cc_GAMMA_VA_REV, S_cc, RULE_O2 );

      complex_distr_t_list RES_S0S = glb_red_FERM4_T(S_cc_S0, RULE_SINGLE);
      complex_distr_t_list RES_V0S = glb_red_FERM4_T(S_cc_V0, RULE_SINGLE);
      complex_distr_t_list RES_V3S = glb_red_FERM4_T(S_cc_V3, RULE_SINGLE);
      complex_distr_t_list RES_A3S = glb_red_FERM4_T(S_cc_A3, RULE_SINGLE);

      complex_distr_t_list RES_S0S_REV = glb_red_FERM4_T(S_cc_S0_REV, RULE_SINGLE);
      complex_distr_t_list RES_V0S_REV = glb_red_FERM4_T(S_cc_V0_REV, RULE_SINGLE);
      complex_distr_t_list RES_V3S_REV = glb_red_FERM4_T(S_cc_V3_REV, RULE_SINGLE);
      complex_distr_t_list RES_A3S_REV = glb_red_FERM4_T(S_cc_A3_REV, RULE_SINGLE);
    
      complex_distr_t_list RES_P5P5_cc = glb_red_FERM4_T(T_cc_P, RULE_SINGLE);

      complex_distr_t_list RES_VKV3_cc = glb_red_FERM4_T(T_cc_VKV3, RULE_SINGLE);
      complex_distr_t_list RES_VKV0_cc = glb_red_FERM4_T(T_cc_VKV0, RULE_SINGLE);
      complex_distr_t_list RES_VKA3_cc = glb_red_FERM4_T(T_cc_VKA3, RULE_SINGLE);
      complex_distr_t_list RES_VKA0_cc = glb_red_FERM4_T(T_cc_VKA0, RULE_SINGLE);

      auto end3 = chrono::system_clock::now();
    
      chrono::duration<double> elapsed_seconds_aC = end3-start;
    
      cout<<"elapsed time after tensor contraction: "<<elapsed_seconds_aC.count()<<" s"<<flush<<endl;

      
    

      cout<<"Adding contraction... "<<flush<<endl;

     
      C_O1 = C_O1 +  0.5*(RES_O1 +GLB_SIGN_REV*RES_O1_REV)/Nhits;
      C_O2 = C_O2 +  0.5*(RES_O2 +GLB_SIGN_REV*RES_O2_REV)/Nhits;

              
      C_O1_even = C_O1_even + 0.5*(RES_O1_AA +RES_O1_VV + GLB_SIGN_REV*RES_O1_AA_REV+ GLB_SIGN_REV*RES_O1_VV_REV)/Nhits;
      C_O2_even = C_O2_even + 0.5*(RES_O2_AA +RES_O2_VV + GLB_SIGN_REV*RES_O2_AA_REV+ GLB_SIGN_REV*RES_O2_VV_REV)/Nhits;

  
      C_O1_odd = C_O1_odd -0.5*(RES_O1_AV +RES_O1_VA + GLB_SIGN_REV*RES_O1_AV_REV + GLB_SIGN_REV*RES_O1_VA_REV)/Nhits;
      C_O2_odd = C_O2_odd -0.5*(RES_O2_AV + RES_O2_VA + GLB_SIGN_REV*RES_O2_AV_REV + GLB_SIGN_REV*RES_O2_VA_REV)/Nhits;

        
      C_O1_VV = C_O1_VV + 0.5*(RES_O1_VV+ GLB_SIGN_REV*RES_O1_VV_REV)/Nhits;
      C_O1_AA = C_O1_AA + 0.5*(RES_O1_AA+ GLB_SIGN_REV*RES_O1_AA_REV)/Nhits;
      C_O1_VA = C_O1_VA + 0.5*(RES_O1_VA+ GLB_SIGN_REV*RES_O1_VA_REV)/Nhits;
      C_O1_AV = C_O1_AV + 0.5*(RES_O1_AV+ GLB_SIGN_REV*RES_O1_VA_REV)/Nhits;

      C_O1_not_ave_VV = C_O1_not_ave_VV + RES_O1_VV/Nhits;
      C_O1_not_ave_AA = C_O1_not_ave_AA + RES_O1_AA/Nhits;

            
      C_O2_VV = C_O2_VV + 0.5*(RES_O2_VV+ GLB_SIGN_REV*RES_O2_VV_REV)/Nhits;
      C_O2_AA = C_O2_AA + 0.5*(RES_O2_AA+ GLB_SIGN_REV*RES_O2_AA_REV)/Nhits;
      C_O2_VA = C_O2_VA + 0.5*(RES_O2_VA+ GLB_SIGN_REV*RES_O2_VA_REV)/Nhits;
      C_O2_AV = C_O2_AV + 0.5*(RES_O2_AV+ GLB_SIGN_REV*RES_O2_VA_REV)/Nhits;

      C_O2_not_ave_VV = C_O2_not_ave_VV + RES_O2_VV/Nhits;
      C_O2_not_ave_AA = C_O2_not_ave_AA + RES_O2_AA/Nhits;

          
      C_semlep_S0 = C_semlep_S0 + 0.5*( RES_S0S + RES_S0S_REV)/Nhits;
      C_semlep_V0 = C_semlep_V0 + 0.5*( RES_V0S + RES_V0S_REV)/Nhits;
      C_semlep_V3 = C_semlep_V3 + 0.5*( RES_V3S - RES_V3S_REV)/Nhits;
      C_semlep_A3 = C_semlep_A3 + RES_A3S/Nhits;
      C_semlep_A3_ave = C_semlep_A3_ave + 0.5*( RES_A3S - RES_A3S_REV)/Nhits;
      
      C_etacc = C_etacc + RES_P5P5_cc/Nhits;
      
      C_cc_VKV3 = C_cc_VKV3 + RES_VKV3_cc/Nhits;
      C_cc_VKV0 = C_cc_VKV0 + RES_VKV0_cc/Nhits;
      C_cc_VKA3 = C_cc_VKA3 + RES_VKA3_cc/Nhits;
      C_cc_VKA0 = C_cc_VKA0 + RES_VKA0_cc/Nhits;

              
      C_O1_VV_FIERZ = C_O1_VV_FIERZ + 0.5*(RES_O1_VV_FIERZ + GLB_SIGN_REV*RES_O1_VV_FIERZ_REV)/Nhits;
      C_O1_AA_FIERZ = C_O1_AA_FIERZ + 0.5*(RES_O1_AA_FIERZ + GLB_SIGN_REV*RES_O1_AA_FIERZ_REV)/Nhits;

      C_O2_VV_FIERZ = C_O2_VV_FIERZ + 0.5*(RES_O2_VV_FIERZ + GLB_SIGN_REV*RES_O2_VV_FIERZ_REV)/Nhits;
      C_O2_AA_FIERZ = C_O2_AA_FIERZ + 0.5*(RES_O2_AA_FIERZ + GLB_SIGN_REV*RES_O2_AA_FIERZ_REV)/Nhits;

      auto end4 = chrono::system_clock::now();
    
      chrono::duration<double> elapsed_seconds_2 = end4-start;
    
      cout<<"Time to compute ihit: "<<ihit<<": "<<elapsed_seconds_2.count()<<" s"<<flush<<endl;
       

    

      cout<<"Done!"<<flush<<endl;

    
    
  }

  cout<<"Contractions computed!"<<flush<<endl;
   

   
  //load 2pt function

  data_t pt2_B, pt2_K;

  Corr.Perform_Nt_t_average=1;
  Corr_JJ.Perform_Nt_t_average=1;

  pt2_B.Read("../PINGU", "mes_contr_PT2_"+PA+"_averaged", "P5P5");
  pt2_K.Read("../PINGU", "mes_contr_PT2_"+DT+"_averaged" , "P5P5");


  
  distr_t_list B_2pt_distr = Corr.corr_t(pt2_B.col(0)[0], "");
  distr_t_list K_2pt_distr = Corr.corr_t(pt2_K.col(0)[0], "");

  distr_t_list B_2pt_distr_JJ = Corr_JJ.corr_t(pt2_B.col(0)[0], "");
  distr_t_list K_2pt_distr_JJ = Corr_JJ.corr_t(pt2_K.col(0)[0], "");
  

  //amputate

  if(Is_Bs) {
    Corr.Tmin=15;
    Corr.Tmax=20; 
  }
  else {
    Corr.Tmin=15;
    Corr.Tmax=20;
  }

  Corr_JJ.Tmin=Corr.Tmin;
  Corr_JJ.Tmax=Corr.Tmax;

  distr_t_list Z_B_distr = Corr.matrix_element_t( B_2pt_distr, "");
  distr_t_list MB_distr = Corr.effective_mass_t(B_2pt_distr, "../data/PINGU/B_eff_mass");

  distr_t_list Z_B_distr_JJ = Corr_JJ.matrix_element_t( B_2pt_distr_JJ, "");
  distr_t_list MB_distr_JJ = Corr_JJ.effective_mass_t(B_2pt_distr_JJ, "");

  distr_t ZB = Corr.Fit_distr(Z_B_distr);
  distr_t MB = Corr.Fit_distr(MB_distr);

  distr_t ZB_JJ = Corr_JJ.Fit_distr(Z_B_distr_JJ);
  distr_t MB_JJ = Corr_JJ.Fit_distr(MB_distr_JJ);

  if(Is_Bs) {
    Corr.Tmin=15;
    Corr.Tmax=45;
  }
  else {
    Corr.Tmin=15;
    Corr.Tmax=21;
  }

  Corr_JJ.Tmin=Corr.Tmin;
  Corr_JJ.Tmax=Corr.Tmax;
  
  distr_t_list Z_K_distr = Corr.matrix_element_t(K_2pt_distr, "");
  distr_t_list MK_distr = Corr.effective_mass_t(K_2pt_distr, "../data/PINGU/K_eff_mass");
  distr_t_list Z_K_distr_JJ = Corr_JJ.matrix_element_t(K_2pt_distr_JJ, "");
  distr_t_list MK_distr_JJ = Corr_JJ.effective_mass_t(K_2pt_distr_JJ, "");

  distr_t ZK = Corr.Fit_distr( Z_K_distr);
  distr_t MK = Corr.Fit_distr( MK_distr);

  distr_t ZK_JJ = Corr_JJ.Fit_distr( Z_K_distr_JJ);
  distr_t MK_JJ = Corr_JJ.Fit_distr( MK_distr_JJ);

  distr_t amp_B = (ZB/(2.0*MB))*EXP_D(-t_O4*MB);
  distr_t amp_B_JJ = (ZB_JJ/(2.0*MB_JJ))*EXP_D(-t_O4*MB_JJ);
  distr_t_list amp_Ksem = (ZK/(2.0*MK))*EXPT_D(-1.0*MK,Corr.Nt)*EXP_D(MK*t_O4);


  double a= 0.0795*fm_to_inv_Gev;

  
 

  Corr.Perform_Nt_t_average=0;
  Corr_JJ.Perform_Nt_t_average=0;
  Corr.Nt= C_O1.size();
  Corr_JJ.Nt=C_O1.size();
  C_O1= complex_distr_t_list(Corr.corr_t( C_O1.Get_vvector(0),""), Corr.corr_t( C_O1.Get_vvector(1), ""));
  C_O2= complex_distr_t_list(Corr.corr_t( C_O2.Get_vvector(0),""), Corr.corr_t( C_O2.Get_vvector(1), ""));

  complex_distr_t_list C_O1_VV_FIERZ_JJ= complex_distr_t_list(Corr_JJ.corr_t( C_O1_VV_FIERZ.Get_vvector(0),""), Corr_JJ.corr_t( C_O1_VV_FIERZ.Get_vvector(1), ""));
  complex_distr_t_list C_O1_AA_FIERZ_JJ= complex_distr_t_list(Corr_JJ.corr_t( C_O1_AA_FIERZ.Get_vvector(0),""), Corr_JJ.corr_t( C_O1_AA_FIERZ.Get_vvector(1), ""));
  complex_distr_t_list C_O2_VV_FIERZ_JJ= complex_distr_t_list(Corr_JJ.corr_t( C_O2_VV_FIERZ.Get_vvector(0),""), Corr_JJ.corr_t( C_O2_VV_FIERZ.Get_vvector(1), ""));
  complex_distr_t_list C_O2_AA_FIERZ_JJ= complex_distr_t_list(Corr_JJ.corr_t( C_O2_AA_FIERZ.Get_vvector(0),""), Corr_JJ.corr_t( C_O2_AA_FIERZ.Get_vvector(1), ""));

  
  C_O1_VV_FIERZ= complex_distr_t_list(Corr.corr_t( C_O1_VV_FIERZ.Get_vvector(0),""), Corr.corr_t( C_O1_VV_FIERZ.Get_vvector(1), ""));
  C_O1_AA_FIERZ= complex_distr_t_list(Corr.corr_t( C_O1_AA_FIERZ.Get_vvector(0),""), Corr.corr_t( C_O1_AA_FIERZ.Get_vvector(1), ""));
  C_O2_VV_FIERZ= complex_distr_t_list(Corr.corr_t( C_O2_VV_FIERZ.Get_vvector(0),""), Corr.corr_t( C_O2_VV_FIERZ.Get_vvector(1), ""));
  C_O2_AA_FIERZ= complex_distr_t_list(Corr.corr_t( C_O2_AA_FIERZ.Get_vvector(0),""), Corr.corr_t( C_O2_AA_FIERZ.Get_vvector(1), ""));

 
  
  C_O1_odd= complex_distr_t_list(Corr.corr_t( C_O1_odd.Get_vvector(0),""), Corr.corr_t( C_O1_odd.Get_vvector(1), ""));
  C_O2_odd= complex_distr_t_list(Corr.corr_t( C_O2_odd.Get_vvector(0),""), Corr.corr_t( C_O2_odd.Get_vvector(1), ""));
  C_O1_even= complex_distr_t_list(Corr.corr_t( C_O1_even.Get_vvector(0),""), Corr.corr_t( C_O1_even.Get_vvector(1), ""));
  C_O2_even= complex_distr_t_list(Corr.corr_t( C_O2_even.Get_vvector(0),""), Corr.corr_t( C_O2_even.Get_vvector(1), ""));

  C_O1_VV= complex_distr_t_list(Corr.corr_t( C_O1_VV.Get_vvector(0),""), Corr.corr_t( C_O1_VV.Get_vvector(1), ""));
  C_O1_AV= complex_distr_t_list(Corr.corr_t( C_O1_AV.Get_vvector(0),""), Corr.corr_t( C_O1_AV.Get_vvector(1), ""));
  C_O1_VA= complex_distr_t_list(Corr.corr_t( C_O1_VA.Get_vvector(0),""), Corr.corr_t( C_O1_VA.Get_vvector(1), ""));
  C_O1_AA= complex_distr_t_list(Corr.corr_t( C_O1_AA.Get_vvector(0),""), Corr.corr_t( C_O1_AA.Get_vvector(1), ""));

  C_O1_not_ave_VV= complex_distr_t_list(Corr.corr_t( C_O1_not_ave_VV.Get_vvector(0),""), Corr.corr_t( C_O1_not_ave_VV.Get_vvector(1), ""));
   C_O1_not_ave_AA= complex_distr_t_list(Corr.corr_t( C_O1_not_ave_AA.Get_vvector(0),""), Corr.corr_t( C_O1_not_ave_AA.Get_vvector(1), ""));

  C_O2_VV= complex_distr_t_list(Corr.corr_t( C_O2_VV.Get_vvector(0),""), Corr.corr_t( C_O2_VV.Get_vvector(1), ""));
  C_O2_AV= complex_distr_t_list(Corr.corr_t( C_O2_AV.Get_vvector(0),""), Corr.corr_t( C_O2_AV.Get_vvector(1), ""));
  C_O2_VA= complex_distr_t_list(Corr.corr_t( C_O2_VA.Get_vvector(0),""), Corr.corr_t( C_O2_VA.Get_vvector(1), ""));
  C_O2_AA= complex_distr_t_list(Corr.corr_t( C_O2_AA.Get_vvector(0),""), Corr.corr_t( C_O2_AA.Get_vvector(1), ""));

  C_O2_not_ave_VV= complex_distr_t_list(Corr.corr_t( C_O2_not_ave_VV.Get_vvector(0),""), Corr.corr_t( C_O2_not_ave_VV.Get_vvector(1), ""));
  C_O2_not_ave_AA= complex_distr_t_list(Corr.corr_t( C_O2_not_ave_AA.Get_vvector(0),""), Corr.corr_t( C_O2_not_ave_AA.Get_vvector(1), ""));

  /*
  C_O1_vacuum_VV= complex_distr_t_list(Corr.corr_t( C_O1_vacuum_VV.Get_vvector(0),""), Corr.corr_t( C_O1_vacuum_VV.Get_vvector(1), ""));
  C_O1_vacuum_AV= complex_distr_t_list(Corr.corr_t( C_O1_vacuum_AV.Get_vvector(0),""), Corr.corr_t( C_O1_vacuum_AV.Get_vvector(1), ""));
  C_O1_vacuum_VA= complex_distr_t_list(Corr.corr_t( C_O1_vacuum_VA.Get_vvector(0),""), Corr.corr_t( C_O1_vacuum_VA.Get_vvector(1), ""));
  C_O1_vacuum_AA= complex_distr_t_list(Corr.corr_t( C_O1_vacuum_AA.Get_vvector(0),""), Corr.corr_t( C_O1_vacuum_AA.Get_vvector(1), ""));
  */

  
  Corr.Nt= C_semlep_S0.size();
  
  C_semlep_S0= V4*complex_distr_t_list(Corr.corr_t( C_semlep_S0.Get_vvector(0),""), Corr.corr_t( C_semlep_S0.Get_vvector(1), ""))/(amp_B*amp_Ksem);
  C_semlep_V0= V4*complex_distr_t_list(Corr.corr_t( C_semlep_V0.Get_vvector(0),""), Corr.corr_t( C_semlep_V0.Get_vvector(1), ""))/(amp_B*amp_Ksem);
  C_semlep_V3= V4*complex_distr_t_list(Corr.corr_t( C_semlep_V3.Get_vvector(0),""), Corr.corr_t( C_semlep_V3.Get_vvector(1), ""))/(amp_B*amp_Ksem);
  C_semlep_A3= V4*complex_distr_t_list(Corr.corr_t( C_semlep_A3.Get_vvector(0),""), Corr.corr_t( C_semlep_A3.Get_vvector(1), ""))/(amp_B*amp_Ksem);
  C_semlep_A3_ave= V4*complex_distr_t_list(Corr.corr_t( C_semlep_A3_ave.Get_vvector(0),""), Corr.corr_t( C_semlep_A3_ave.Get_vvector(1), ""))/(amp_B*amp_Ksem);
    
  C_etacc= complex_distr_t_list(Corr.corr_t( C_etacc.Get_vvector(0),""), Corr.corr_t( C_etacc.Get_vvector(1), ""));

  C_cc_VKV3 = complex_distr_t_list( Corr.corr_t( C_cc_VKV3.Get_vvector(0), ""), Corr.corr_t( C_cc_VKV3.Get_vvector(1), ""));
  C_cc_VKV0 = complex_distr_t_list( Corr.corr_t( C_cc_VKV0.Get_vvector(0), ""), Corr.corr_t( C_cc_VKV0.Get_vvector(1), ""));
  C_cc_VKA3 = complex_distr_t_list( Corr.corr_t( C_cc_VKA3.Get_vvector(0), ""), Corr.corr_t( C_cc_VKA3.Get_vvector(1), ""));
  C_cc_VKA0 = complex_distr_t_list( Corr.corr_t( C_cc_VKA0.Get_vvector(0), ""), Corr.corr_t( C_cc_VKA0.Get_vvector(1), ""));



  Print_To_File({}, { C_semlep_A3.RE_ave(), C_semlep_A3.RE_err(), C_semlep_A3.IM_ave(), C_semlep_A3.IM_err() }, "../data/PINGU/C_sem_A3","","");
  Print_To_File({}, { C_semlep_A3_ave.RE_ave(), C_semlep_A3_ave.RE_err(), C_semlep_A3_ave.IM_ave(), C_semlep_A3_ave.IM_err() }, "../data/PINGU/C_sem_A3_ave","","");
  Print_To_File({}, { C_semlep_V3.RE_ave(), C_semlep_V3.RE_err(), C_semlep_V3.IM_ave(), C_semlep_V3.IM_err() }, "../data/PINGU/C_sem_V3","","");

  
   
  //fit time

  vector<complex_distr_t_list> C_O2_ord;
  vector<complex_distr_t_list> C_O1_ord;
  vector<complex_distr_t_list> C_O1_VV_FIERZ_ord;
  vector<complex_distr_t_list> C_O1_AA_FIERZ_ord;
  vector<complex_distr_t_list> C_O2_VV_FIERZ_ord;
  vector<complex_distr_t_list> C_O2_AA_FIERZ_ord;
  vector<complex_distr_t_list> C_O1_VV_FIERZ_ord_JJ;
  vector<complex_distr_t_list> C_O1_AA_FIERZ_ord_JJ;
  vector<complex_distr_t_list> C_O2_VV_FIERZ_ord_JJ;
  vector<complex_distr_t_list> C_O2_AA_FIERZ_ord_JJ;
  vector<complex_distr_t_list> C_O2_odd_ord;
  vector<complex_distr_t_list> C_O1_odd_ord;
  vector<complex_distr_t_list> C_O2_even_ord;
  vector<complex_distr_t_list> C_O1_even_ord;

  vector<complex_distr_t_list> C_O1_VV_ord, C_O1_AA_ord, C_O1_AV_ord, C_O1_VA_ord;
  vector<complex_distr_t_list> C_O2_VV_ord, C_O2_AA_ord, C_O2_AV_ord, C_O2_VA_ord;

  vector<complex_distr_t_list> C_O1_not_ave_VV_ord, C_O1_not_ave_AA_ord;
  vector<complex_distr_t_list> C_O2_not_ave_VV_ord, C_O2_not_ave_AA_ord;

   
  for(int t=0;t<T;t++) {
    C_O2_ord.emplace_back(UseJack); C_O1_ord.emplace_back(UseJack);
    C_O1_AA_FIERZ_ord.emplace_back(UseJack); C_O1_VV_FIERZ_ord.emplace_back(UseJack);
    C_O2_AA_FIERZ_ord.emplace_back(UseJack); C_O2_VV_FIERZ_ord.emplace_back(UseJack);
    C_O1_AA_FIERZ_ord_JJ.emplace_back(UseJack); C_O1_VV_FIERZ_ord_JJ.emplace_back(UseJack);
    C_O2_AA_FIERZ_ord_JJ.emplace_back(UseJack); C_O2_VV_FIERZ_ord_JJ.emplace_back(UseJack);
    C_O2_odd_ord.emplace_back(UseJack); C_O1_odd_ord.emplace_back(UseJack);
    C_O2_even_ord.emplace_back(UseJack); C_O1_even_ord.emplace_back(UseJack);

    C_O1_VV_ord.emplace_back(UseJack); C_O1_AA_ord.emplace_back(UseJack); C_O1_AV_ord.emplace_back(UseJack); C_O1_VA_ord.emplace_back(UseJack);
    C_O2_VV_ord.emplace_back(UseJack); C_O2_AA_ord.emplace_back(UseJack); C_O2_AV_ord.emplace_back(UseJack); C_O2_VA_ord.emplace_back(UseJack);

    C_O1_not_ave_VV_ord.emplace_back(UseJack); C_O1_not_ave_AA_ord.emplace_back(UseJack); 
    C_O2_not_ave_VV_ord.emplace_back(UseJack); C_O2_not_ave_AA_ord.emplace_back(UseJack);

        
  }

  for(int ts=0;ts<T;ts++) {

   

    for(int tem=0; tem<T;tem++) {

      int tt = tem + ts*T;

      distr_t amp_K = (ZK/(2.0*MK))*EXP_D(-1.0*MK*(ts-tem));
     
      
      distr_t amp_K_JJ = (ZK_JJ/(2.0*MK_JJ))*EXP_D(-1.0*MK_JJ*(ts-tem));

      C_O1_ord[ts].distr_list.push_back( V4*V4*C_O1.distr_list[tt]/(amp_B*amp_K) );
      C_O2_ord[ts].distr_list.push_back( V4*V4*C_O2.distr_list[tt]/(amp_B*amp_K) );

      C_O1_odd_ord[ts].distr_list.push_back( V4*V4*C_O1_odd.distr_list[tt]/(amp_B*amp_K) );
      C_O2_odd_ord[ts].distr_list.push_back( V4*V4*C_O2_odd.distr_list[tt]/(amp_B*amp_K) );


      C_O1_even_ord[ts].distr_list.push_back( V4*V4*C_O1_even.distr_list[tt]/(amp_B*amp_K) );
      C_O2_even_ord[ts].distr_list.push_back( V4*V4*C_O2_even.distr_list[tt]/(amp_B*amp_K) );

      C_O1_VV_FIERZ_ord[ts].distr_list.push_back( V4*V4*C_O1_VV_FIERZ.distr_list[tt]/(amp_B*amp_K) );
      C_O1_AA_FIERZ_ord[ts].distr_list.push_back( V4*V4*C_O1_AA_FIERZ.distr_list[tt]/(amp_B*amp_K) );
      
      C_O2_VV_FIERZ_ord[ts].distr_list.push_back( V4*V4*C_O2_VV_FIERZ.distr_list[tt]/(amp_B*amp_K) );
      C_O2_AA_FIERZ_ord[ts].distr_list.push_back( V4*V4*C_O2_AA_FIERZ.distr_list[tt]/(amp_B*amp_K) );

      C_O1_VV_FIERZ_ord_JJ[ts].distr_list.push_back( V4*V4*C_O1_VV_FIERZ_JJ.distr_list[tt]/(amp_B_JJ*amp_K_JJ) );
      C_O1_AA_FIERZ_ord_JJ[ts].distr_list.push_back( V4*V4*C_O1_AA_FIERZ_JJ.distr_list[tt]/(amp_B_JJ*amp_K_JJ) );
      
      C_O2_VV_FIERZ_ord_JJ[ts].distr_list.push_back( V4*V4*C_O2_VV_FIERZ_JJ.distr_list[tt]/(amp_B_JJ*amp_K_JJ) );
      C_O2_AA_FIERZ_ord_JJ[ts].distr_list.push_back( V4*V4*C_O2_AA_FIERZ_JJ.distr_list[tt]/(amp_B_JJ*amp_K_JJ) );

      C_O1_VV_ord[ts].distr_list.push_back( V4*V4*C_O1_VV.distr_list[tt]/(amp_B*amp_K));
      C_O1_AV_ord[ts].distr_list.push_back( V4*V4*C_O1_AV.distr_list[tt]/(amp_B*amp_K));
      C_O1_AA_ord[ts].distr_list.push_back( V4*V4*C_O1_AA.distr_list[tt]/(amp_B*amp_K));
      C_O1_VA_ord[ts].distr_list.push_back( V4*V4*C_O1_VA.distr_list[tt]/(amp_B*amp_K));

      C_O1_not_ave_VV_ord[ts].distr_list.push_back( V4*V4*C_O1_not_ave_VV.distr_list[tt]/(amp_B*amp_K));
      C_O1_not_ave_AA_ord[ts].distr_list.push_back( V4*V4*C_O1_not_ave_AA.distr_list[tt]/(amp_B*amp_K));
    

      C_O2_VV_ord[ts].distr_list.push_back( V4*V4*C_O2_VV.distr_list[tt]/(amp_B*amp_K));
      C_O2_AV_ord[ts].distr_list.push_back( V4*V4*C_O2_AV.distr_list[tt]/(amp_B*amp_K));
      C_O2_AA_ord[ts].distr_list.push_back( V4*V4*C_O2_AA.distr_list[tt]/(amp_B*amp_K));
      C_O2_VA_ord[ts].distr_list.push_back( V4*V4*C_O2_VA.distr_list[tt]/(amp_B*amp_K));

      C_O2_not_ave_VV_ord[ts].distr_list.push_back( V4*V4*C_O2_not_ave_VV.distr_list[tt]/(amp_B*amp_K));
      C_O2_not_ave_AA_ord[ts].distr_list.push_back( V4*V4*C_O2_not_ave_AA.distr_list[tt]/(amp_B*amp_K));


    }
  }

 
 
  cout<<"MB: "<<MB.ave()/a<<" "<<MB.err()/a<<endl;
  cout<<"MK: "<<MK.ave()/a<<" "<<MK.err()/a<<endl;

  //cout<<"ampB: "<<amp_B.ave()<<" "<<amp_B.err()<<endl;


  C_etacc = C_etacc*V4;

  C_cc_VKV3 = C_cc_VKV3*V4;
  C_cc_VKV0 = C_cc_VKV0*V4;
  C_cc_VKA3 = C_cc_VKA3*V4;
  C_cc_VKA0 = C_cc_VKA0*V4;

  Corr.Perform_Nt_t_average=0;

  distr_t_list C_etac_NEW(UseJack);
  distr_t_list C_etac_NEW_IM(UseJack);

  for(int t=0;t<T;t++) C_etac_NEW.distr_list.push_back( C_etacc.distr_list[(t+int_t_O4)%T].RE );
  for(int t=0;t<T;t++) C_etac_NEW_IM.distr_list.push_back( C_etacc.distr_list[(t+int_t_O4)%T].IM );

  distr_t_list C_cc_VKV3_NEW(UseJack), C_cc_VKV3_NEW_IM(UseJack);
  distr_t_list C_cc_VKV0_NEW(UseJack), C_cc_VKV0_NEW_IM(UseJack);
  distr_t_list C_cc_VKA3_NEW(UseJack), C_cc_VKA3_NEW_IM(UseJack);
  distr_t_list C_cc_VKA0_NEW(UseJack), C_cc_VKA0_NEW_IM(UseJack);

 
  for(int t=0;t<T;t++) {

    C_cc_VKV3_NEW.distr_list.push_back( C_cc_VKV3.distr_list[(t+int_t_O4)%T].RE);  C_cc_VKV3_NEW_IM.distr_list.push_back( C_cc_VKV3.distr_list[(t+int_t_O4)%T].IM);
    C_cc_VKV0_NEW.distr_list.push_back( C_cc_VKV0.distr_list[(t+int_t_O4)%T].RE);  C_cc_VKV0_NEW_IM.distr_list.push_back( C_cc_VKV0.distr_list[(t+int_t_O4)%T].IM);
    C_cc_VKA3_NEW.distr_list.push_back( C_cc_VKA3.distr_list[(t+int_t_O4)%T].RE);  C_cc_VKA3_NEW_IM.distr_list.push_back( C_cc_VKA3.distr_list[(t+int_t_O4)%T].IM);
    C_cc_VKA0_NEW.distr_list.push_back( C_cc_VKA0.distr_list[(t+int_t_O4)%T].RE);  C_cc_VKA0_NEW_IM.distr_list.push_back( C_cc_VKA0.distr_list[(t+int_t_O4)%T].IM);

 
  }
  
  
  

  //symmetryze

  for(int t=0; t <= T/2; t++) {

    distr_t x = 0.5*( C_etac_NEW[t] + C_etac_NEW[(T-t)%T] );
    C_etac_NEW.distr_list[t] = x;
    C_etac_NEW.distr_list[(T-t)%T] = x;
  }

  

  Corr.Tmin=20;
  Corr.Tmax=40;

  
  distr_t_list f_etac = Corr.decay_constant_t( pow(0.463134,2)*C_etac_NEW , "../data/PINGU/f_etac");

  distr_t_list metac = Corr.effective_mass_t( C_etac_NEW, "../data/PINGU/m_etac");
  distr_t_list metac_IM = Corr.effective_mass_t( C_etac_NEW_IM, "../data/PINGU/m_etac_IM");

  distr_t_list m_V3V3 = Corr.effective_mass_t( C_cc_VKV3_NEW, "../data/PINGU/m_V3V3"); distr_t_list m_V3V3_IM = Corr.effective_mass_t( C_cc_VKV3_NEW_IM, "../data/PINGU/m_V3V3_IM");
  distr_t_list m_V3V0 = Corr.effective_mass_t( C_cc_VKV0_NEW, "../data/PINGU/m_V3V0"); distr_t_list m_V3V0_IM = Corr.effective_mass_t( C_cc_VKV0_NEW_IM, "../data/PINGU/m_V3V0_IM");
  distr_t_list m_V3A3 = Corr.effective_mass_t( C_cc_VKA3_NEW, "../data/PINGU/m_V3A3"); distr_t_list m_V3A3_IM = Corr.effective_mass_t( C_cc_VKA3_NEW_IM, "../data/PINGU/m_V3A3_IM");
  distr_t_list m_V3A0 = Corr.effective_mass_t( C_cc_VKA0_NEW, "../data/PINGU/m_V3A0"); distr_t_list m_V3A0_IM = Corr.effective_mass_t( C_cc_VKA0_NEW_IM, "../data/PINGU/m_V3A0_IM");

  Print_To_File({} , { C_cc_VKV3_NEW.ave(), C_cc_VKV3_NEW.err(), C_cc_VKV3_NEW_IM.ave(), C_cc_VKV3_NEW_IM.err() }, "../data/PINGU/C_cc_V3V3","","");
  Print_To_File({} , { C_cc_VKV0_NEW.ave(), C_cc_VKV0_NEW.err(), C_cc_VKV0_NEW_IM.ave(), C_cc_VKV3_NEW_IM.err() }, "../data/PINGU/C_cc_V3V0","","");
  Print_To_File({} , { C_cc_VKA3_NEW.ave(), C_cc_VKA3_NEW.err(), C_cc_VKA3_NEW_IM.ave(), C_cc_VKA3_NEW_IM.err() }, "../data/PINGU/C_cc_V3A3","","");
  Print_To_File({} , { C_cc_VKA0_NEW.ave(), C_cc_VKA0_NEW.err(), C_cc_VKA0_NEW_IM.ave(), C_cc_VKA0_NEW_IM.err() }, "../data/PINGU/C_cc_V3A0","","");
  


  //construct semileptonic form factors
  double aml=0.00072;
  double amh=0.463134;
  double ams=0.0182782;
  double ak= 2.0*M_PI/64;
  distr_t MK0 = SQRT_D( MK*MK - ak*ak);

  distr_t_list f0_WI = C_semlep_S0.RE()*(amh - ams)/(MB*MB -MK0*MK0);
  distr_t_list f_par = -1.0*C_semlep_V0.RE();
  distr_t_list f_perp = (1.0/ak)*C_semlep_V3.IM();
  distr_t_list f0 = (MB-MK)*f_par/(MB*MB - MK0*MK0) + (ak*ak)*f_perp/(MB*MB - MK0*MK0);
  distr_t_list f_plus = (0.5/MB)*f_par + (0.5/MB)*(MB-MK)*f_perp;

  Print_To_File({}, { f0_WI.ave(), f0_WI.err()}, "../data/PINGU/f0_WI", "", "");
  Print_To_File({}, { f_par.ave(), f_par.err()}, "../data/PINGU/f_par", "", "");
  Print_To_File({}, { f_perp.ave(), f_perp.err()}, "../data/PINGU/f_perp", "", "");
  Print_To_File({}, { f0.ave(), f0.err()}, "../data/PINGU/f0", "", "");
  Print_To_File({}, { f_plus.ave(), f_plus.err()}, "../data/PINGU/f_plus", "", "");

  Corr.Perform_Nt_t_average=1;


  //print
  for(int ts=0;ts<T;ts++) {

    complex_distr_t_list C_Oplus_ord= C_O1_ord[ts] + C_O2_ord[ts];
    complex_distr_t_list C_Ominus_ord = C_O1_ord[ts] - C_O2_ord[ts];

    complex_distr_t_list C_O1_even_FIERZ= C_O1_VV_FIERZ_ord[ts] + C_O1_AA_FIERZ_ord[ts];
    complex_distr_t_list C_O2_even_FIERZ= C_O2_VV_FIERZ_ord[ts] + C_O2_AA_FIERZ_ord[ts];
    
    
    Print_To_File({}, {C_O1_ord[ts].RE_ave(), C_O1_ord[ts].RE_err(), C_O1_ord[ts].IM_ave(), C_O1_ord[ts].IM_err() },    "../data/PINGU/C_O1_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O2_ord[ts].RE_ave(), C_O2_ord[ts].RE_err(), C_O2_ord[ts].IM_ave(), C_O2_ord[ts].IM_err() },    "../data/PINGU/C_O2_ts_"+to_string(ts) , "", "");
    
    Print_To_File({}, {C_O1_odd_ord[ts].RE_ave(), C_O1_odd_ord[ts].RE_err(), C_O1_odd_ord[ts].IM_ave(), C_O1_odd_ord[ts].IM_err() },    "../data/PINGU/C_O1_odd_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O2_odd_ord[ts].RE_ave(), C_O2_odd_ord[ts].RE_err(), C_O2_odd_ord[ts].IM_ave(), C_O2_odd_ord[ts].IM_err() },    "../data/PINGU/C_O2_odd_ts_"+to_string(ts) , "", "");

    Print_To_File({}, {C_O1_even_ord[ts].RE_ave(), C_O1_even_ord[ts].RE_err(), C_O1_even_ord[ts].IM_ave(), C_O1_even_ord[ts].IM_err() },    "../data/PINGU/C_O1_even_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O2_even_ord[ts].RE_ave(), C_O2_even_ord[ts].RE_err(), C_O2_even_ord[ts].IM_ave(), C_O2_even_ord[ts].IM_err() },    "../data/PINGU/C_O2_even_ts_"+to_string(ts) , "", "");

    Print_To_File({}, {C_O1_VV_FIERZ_ord[ts].RE_ave(), C_O1_VV_FIERZ_ord[ts].RE_err(), C_O1_VV_FIERZ_ord[ts].IM_ave(), C_O1_VV_FIERZ_ord[ts].IM_err() },    "../data/PINGU/C_O1_VV_FIERZ_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O1_AA_FIERZ_ord[ts].RE_ave(), C_O1_AA_FIERZ_ord[ts].RE_err(), C_O1_AA_FIERZ_ord[ts].IM_ave(), C_O1_AA_FIERZ_ord[ts].IM_err() },    "../data/PINGU/C_O1_AA_FIERZ_ts_"+to_string(ts) , "", "");

    Print_To_File({}, {C_O1_even_FIERZ.RE_ave(), C_O1_even_FIERZ.RE_err(), C_O1_even_FIERZ.IM_ave(), C_O1_even_FIERZ.IM_err() },    "../data/PINGU/C_O1_even_FIERZ_ts_"+to_string(ts) , "", "");

    Print_To_File({}, {C_O2_even_FIERZ.RE_ave(), C_O2_even_FIERZ.RE_err(), C_O2_even_FIERZ.IM_ave(), C_O2_even_FIERZ.IM_err() },    "../data/PINGU/C_O2_even_FIERZ_ts_"+to_string(ts) , "", "");

    Print_To_File({}, {C_O2_VV_FIERZ_ord[ts].RE_ave(), C_O2_VV_FIERZ_ord[ts].RE_err(), C_O2_VV_FIERZ_ord[ts].IM_ave(), C_O2_VV_FIERZ_ord[ts].IM_err() },    "../data/PINGU/C_O2_VV_FIERZ_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O2_AA_FIERZ_ord[ts].RE_ave(), C_O2_AA_FIERZ_ord[ts].RE_err(), C_O2_AA_FIERZ_ord[ts].IM_ave(), C_O2_AA_FIERZ_ord[ts].IM_err() },    "../data/PINGU/C_O2_AA_FIERZ_ts_"+to_string(ts) , "", "");

    Print_To_File({}, {C_O1_VV_ord[ts].RE_ave(), C_O1_VV_ord[ts].RE_err(), C_O1_VV_ord[ts].IM_ave(), C_O1_VV_ord[ts].IM_err() },    "../data/PINGU/C_O1_VV_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O1_AA_ord[ts].RE_ave(), C_O1_AA_ord[ts].RE_err(), C_O1_AA_ord[ts].IM_ave(), C_O1_AA_ord[ts].IM_err() },    "../data/PINGU/C_O1_AA_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O1_AV_ord[ts].RE_ave(), C_O1_AV_ord[ts].RE_err(), C_O1_AV_ord[ts].IM_ave(), C_O1_AV_ord[ts].IM_err() },    "../data/PINGU/C_O1_AV_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O1_VA_ord[ts].RE_ave(), C_O1_VA_ord[ts].RE_err(), C_O1_VA_ord[ts].IM_ave(), C_O1_VA_ord[ts].IM_err() },    "../data/PINGU/C_O1_VA_ts_"+to_string(ts) , "", "");
    
    Print_To_File({}, {C_O1_not_ave_VV_ord[ts].RE_ave(), C_O1_not_ave_VV_ord[ts].RE_err(), C_O1_not_ave_VV_ord[ts].IM_ave(), C_O1_not_ave_VV_ord[ts].IM_err() },    "../data/PINGU/C_O1_not_ave_VV_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O1_not_ave_AA_ord[ts].RE_ave(), C_O1_not_ave_AA_ord[ts].RE_err(), C_O1_not_ave_AA_ord[ts].IM_ave(), C_O1_not_ave_AA_ord[ts].IM_err() },    "../data/PINGU/C_O1_not_ave_AA_ts_"+to_string(ts) , "", "");

    Print_To_File({}, {C_O2_VV_ord[ts].RE_ave(), C_O2_VV_ord[ts].RE_err(), C_O2_VV_ord[ts].IM_ave(), C_O2_VV_ord[ts].IM_err() },    "../data/PINGU/C_O2_VV_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O2_AA_ord[ts].RE_ave(), C_O2_AA_ord[ts].RE_err(), C_O2_AA_ord[ts].IM_ave(), C_O2_AA_ord[ts].IM_err() },    "../data/PINGU/C_O2_AA_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O2_AV_ord[ts].RE_ave(), C_O2_AV_ord[ts].RE_err(), C_O2_AV_ord[ts].IM_ave(), C_O2_AV_ord[ts].IM_err() },    "../data/PINGU/C_O2_AV_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O2_VA_ord[ts].RE_ave(), C_O2_VA_ord[ts].RE_err(), C_O2_VA_ord[ts].IM_ave(), C_O2_VA_ord[ts].IM_err() },    "../data/PINGU/C_O2_VA_ts_"+to_string(ts) , "", "");
    
    Print_To_File({}, {C_O2_not_ave_VV_ord[ts].RE_ave(), C_O2_not_ave_VV_ord[ts].RE_err(), C_O2_not_ave_VV_ord[ts].IM_ave(), C_O2_not_ave_VV_ord[ts].IM_err() },    "../data/PINGU/C_O2_not_ave_VV_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_O2_not_ave_AA_ord[ts].RE_ave(), C_O2_not_ave_AA_ord[ts].RE_err(), C_O2_not_ave_AA_ord[ts].IM_ave(), C_O2_not_ave_AA_ord[ts].IM_err() },    "../data/PINGU/C_O2_not_ave_AA_ts_"+to_string(ts) , "", "");

       
    Print_To_File({}, {C_Oplus_ord.RE_ave(), C_Oplus_ord.RE_err(), C_Oplus_ord.IM_ave(), C_Oplus_ord.IM_err() },    "../data/PINGU/C_Oplus_ts_"+to_string(ts) , "", "");
    Print_To_File({}, {C_Ominus_ord.RE_ave(), C_Ominus_ord.RE_err(), C_Ominus_ord.IM_ave(), C_Ominus_ord.IM_err() },    "../data/PINGU/C_Ominus_ts_"+to_string(ts) , "", "");

    distr_t_list C_O1_meff(UseJack), C_O2_meff(UseJack), C_O1_even_meff(UseJack), C_O1_odd_meff(UseJack), C_O2_even_meff(UseJack), C_O2_odd_meff(UseJack);
    distr_t_list C_O1_meff_IM(UseJack), C_O2_meff_IM(UseJack), C_O1_even_meff_IM(UseJack), C_O1_odd_meff_IM(UseJack), C_O2_even_meff_IM(UseJack), C_O2_odd_meff_IM(UseJack);

    for(int t=0;t<T;t++) {
  
	C_O1_meff.distr_list.push_back( LOG_D( C_O1_ord[ts].distr_list[(t+int_t_O4)%T].RE/C_O1_ord[ts].distr_list[ (t+int_t_O4+1)%T ].RE ));
	C_O2_meff.distr_list.push_back( LOG_D( C_O2_ord[ts].distr_list[(t+int_t_O4)%T].RE/C_O2_ord[ts].distr_list[ (t+int_t_O4+1)%T ].RE ));
	C_O1_even_meff.distr_list.push_back( LOG_D( C_O1_even_ord[ts].distr_list[(t+int_t_O4)%T].RE/C_O1_even_ord[ts].distr_list[ (t+int_t_O4+1)%T ].RE ));
	C_O1_odd_meff.distr_list.push_back( LOG_D( C_O1_odd_ord[ts].distr_list[(t+int_t_O4)%T].RE/C_O1_odd_ord[ts].distr_list[ (t+int_t_O4+1)%T ].RE ));
	C_O2_even_meff.distr_list.push_back( LOG_D( C_O2_even_ord[ts].distr_list[(t+int_t_O4)%T].RE/C_O2_even_ord[ts].distr_list[ (t+int_t_O4+1)%T ].RE ));
	C_O2_odd_meff.distr_list.push_back( LOG_D( C_O2_odd_ord[ts].distr_list[(t+int_t_O4)%T].RE/C_O2_odd_ord[ts].distr_list[ (t+int_t_O4+1)%T ].RE ));

	C_O1_meff_IM.distr_list.push_back( LOG_D( C_O1_ord[ts].distr_list[(t+int_t_O4)%T].IM/C_O1_ord[ts].distr_list[ (t+int_t_O4+1)%T ].IM ));
	C_O2_meff_IM.distr_list.push_back( LOG_D( C_O2_ord[ts].distr_list[(t+int_t_O4)%T].IM/C_O2_ord[ts].distr_list[ (t+int_t_O4+1)%T ].IM ));
	C_O1_even_meff_IM.distr_list.push_back( LOG_D( C_O1_even_ord[ts].distr_list[(t+int_t_O4)%T].IM/C_O1_even_ord[ts].distr_list[ (t+int_t_O4+1)%T ].IM ));
	C_O1_odd_meff_IM.distr_list.push_back( LOG_D( C_O1_odd_ord[ts].distr_list[(t+int_t_O4)%T].IM/C_O1_odd_ord[ts].distr_list[ (t+int_t_O4+1)%T ].IM ));
	C_O2_even_meff_IM.distr_list.push_back( LOG_D( C_O2_even_ord[ts].distr_list[(t+int_t_O4)%T].IM/C_O2_even_ord[ts].distr_list[ (t+int_t_O4+1)%T ].IM ));
	C_O2_odd_meff_IM.distr_list.push_back( LOG_D( C_O2_odd_ord[ts].distr_list[(t+int_t_O4)%T].IM/C_O2_odd_ord[ts].distr_list[ (t+int_t_O4+1)%T ].IM ));
      
    
    }

    Print_To_File({}, {C_O1_meff.ave(), C_O1_meff.err(), C_O2_meff.ave(), C_O2_meff.err() },  "../data/PINGU/meff_ts_"+to_string(ts) , "", "");

    Print_To_File({}, {C_O1_even_meff.ave(), C_O1_even_meff.err(), C_O2_even_meff.ave(), C_O2_even_meff.err() },  "../data/PINGU/meff_even_ts_"+to_string(ts) , "", "");

    Print_To_File({}, {C_O1_odd_meff.ave(), C_O1_odd_meff.err(), C_O2_odd_meff.ave(), C_O2_odd_meff.err() },  "../data/PINGU/meff_odd_ts_"+to_string(ts) , "", "");

    Print_To_File({}, {C_O1_meff_IM.ave(), C_O1_meff_IM.err(), C_O2_meff_IM.ave(), C_O2_meff_IM.err() },  "../data/PINGU/meff_IM_ts_"+to_string(ts) , "", "");
     
    Print_To_File({}, {C_O1_even_meff_IM.ave(), C_O1_even_meff_IM.err(), C_O2_even_meff_IM.ave(), C_O2_even_meff_IM.err() },  "../data/PINGU/meff_IM_even_ts_"+to_string(ts) , "", "");

    Print_To_File({}, {C_O1_odd_meff_IM.ave(), C_O1_odd_meff_IM.err(), C_O2_odd_meff_IM.ave(), C_O2_odd_meff_IM.err() },  "../data/PINGU/meff_IM_odd_ts_"+to_string(ts) , "", "");

   
    
  }


  Print_To_File({}, {C_etacc.RE_ave(), C_etacc.RE_err(), C_etacc.IM_ave(), C_etacc.IM_err() }, "../data/PINGU/C_etacc", "", "");

  Print_To_File({}, {C_cc_VKV3.RE_ave(), C_cc_VKV3.RE_err(), C_cc_VKV3.IM_ave(), C_cc_VKV3.IM_err() }, "../data/PINGU/C_VKV3", "", "");
  Print_To_File({}, {C_cc_VKV0.RE_ave(), C_cc_VKV0.RE_err(), C_cc_VKV0.IM_ave(), C_cc_VKV0.IM_err() }, "../data/PINGU/C_VKV0", "", "");
  Print_To_File({}, {C_cc_VKA3.RE_ave(), C_cc_VKA3.RE_err(), C_cc_VKA3.IM_ave(), C_cc_VKA3.IM_err() }, "../data/PINGU/C_VKA3", "", "");
  Print_To_File({}, {C_cc_VKA0.RE_ave(), C_cc_VKA0.RE_err(), C_cc_VKA0.IM_ave(), C_cc_VKA0.IM_err() }, "../data/PINGU/C_VKA0", "", "");
  

  
  cout<<"Bye:"<<endl;


  //store the correlator

  boost::filesystem::create_directory("../data/PINGU/jackknife");

  boost::filesystem::create_directory("../data/PINGU/jackknife/"+HS);
  for(int ts=0;ts<T;ts++)
    for(int t=0;t<Corr.Nt;t++) {
      
     
      ofstream print_O1("../data/PINGU/jackknife/"+HS+"/C_O1_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat");
      print_O1.precision(10);  
      for(int ijack=0;ijack<Nconfs;ijack++) print_O1 << C_O1_AA_FIERZ_ord[ts].distr_list[t].RE.distr[ijack]<<" "<<C_O1_AA_FIERZ_ord[ts].distr_list[t].IM.distr[ijack]<<" "<<C_O1_VV_FIERZ_ord[ts].distr_list[t].RE.distr[ijack]<<" "<<C_O1_VV_FIERZ_ord[ts].distr_list[t].IM.distr[ijack]<<endl;
      print_O1.close();

    
      ofstream print_O2("../data/PINGU/jackknife/"+HS+"/C_O2_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat");
      print_O2.precision(10);
      for(int ijack=0;ijack<Nconfs;ijack++) print_O2 << C_O2_AA_FIERZ_ord[ts].distr_list[t].RE.distr[ijack]<<" "<<C_O2_AA_FIERZ_ord[ts].distr_list[t].IM.distr[ijack]<<" "<<C_O2_VV_FIERZ_ord[ts].distr_list[t].RE.distr[ijack]<<" "<<C_O2_VV_FIERZ_ord[ts].distr_list[t].IM.distr[ijack]<<endl;
      print_O2.close();


      print_O1.open("../data/PINGU/jackknife/"+HS+"/C_JJ_O1_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat");
      print_O1.precision(10);  
      for(int ijack=0;ijack<Njacks_JJ;ijack++) print_O1 << C_O1_AA_FIERZ_ord_JJ[ts].distr_list[t].RE.distr[ijack]<<" "<<C_O1_AA_FIERZ_ord_JJ[ts].distr_list[t].IM.distr[ijack]<<" "<<C_O1_VV_FIERZ_ord_JJ[ts].distr_list[t].RE.distr[ijack]<<" "<<C_O1_VV_FIERZ_ord_JJ[ts].distr_list[t].IM.distr[ijack]<<endl;
      print_O1.close();

    
      print_O2.open("../data/PINGU/jackknife/"+HS+"/C_JJ_O2_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat");
      print_O2.precision(10);
      for(int ijack=0;ijack<Njacks_JJ;ijack++) print_O2 << C_O2_AA_FIERZ_ord_JJ[ts].distr_list[t].RE.distr[ijack]<<" "<<C_O2_AA_FIERZ_ord_JJ[ts].distr_list[t].IM.distr[ijack]<<" "<<C_O2_VV_FIERZ_ord_JJ[ts].distr_list[t].RE.distr[ijack]<<" "<<C_O2_VV_FIERZ_ord_JJ[ts].distr_list[t].IM.distr[ijack]<<endl;
      print_O2.close();


  }
    
    
      

  return;
}




void penguin_spectral_reco() {


  //define model for charming penguin

  auto model_PENGUIN = [](double E, double eps, string mode, int reim) {

    //Vfloat Masses({3.0969, 3.68610,  3.7737, 4.039}); //GeV
    //Vfloat Gammas({ 0.0926e-3, 2.94e-4, 2.72e-2, 8e-2}); //GeV
    //Vfloat Bs({0.05961, 8e-3, 9.6e-6, 10.7e-6 });

    Vfloat Bs({0.05961, 8e-3, 9.6e-6, 10.7e-6, 6.9e-6 });
    Vfloat Gammas({ 0.0926e-3, 2.94e-4, 2.72e-2, 8e-2, 7e-2}); //GeV
    Vfloat Masses({3.05, 3.68610,  3.7737, 4.039, 4.191}); //GeV

    double mh= (mode=="HSS")?mhS:mhL;
    double mp= (mode=="HSS")?meta:mk;
    
    double F_plus= 1.44;   //extracted from data at mh~3 GeV  k =0.244 MeV and a ~ 0.0795 fm
    double F_zero = 0.87;  //extracted from data at mh~3 GeV  k =0.244 MeV and a ~ 0.0795 fm

    double Zv= 0.70637654;

    double sum=0.0;

    double Ep = sqrt( pow(mp,2) + pow(mom,2));
    Vfloat pbmu({ mh, 0.0});
    Vfloat pkmu({ sqrt( pow(mp,2) + pow(mom,2)), -mom});
    Vfloat qmu({ pbmu[0]-pkmu[0], pbmu[1]-pkmu[1]});
    double q2= qmu[0]*qmu[0] - qmu[1]*qmu[1];
    double ad=  0.07948*fm_to_inv_Gev;  //lattice spacing in GeV^-1

    //cout<<"mode: "<<mode<<" mh: "<<mh<<" mp: "<<mp<<endl;

    for(int im=0; im< (signed)Masses.size();im++) {

      double mv= Masses[im];
      double Ev= sqrt( mv*mv + mom*mom);
      double G= Gammas[im];
      double fV= (3.0/2.0)*sqrt((3.0/(4.0*M_PI*alpha*alpha))*Gammas[im]*Bs[im]*mv);

      //cout<<"fV["<<Masses[im]<<"]: "<<fV<<endl;
   
      
      Vfloat kmu({sqrt( pow(mv,2)+ pow(mom,2)), mom});
      double kq= kmu[0]*qmu[0] - kmu[1]*qmu[1];
      double kpb= kmu[0]*pbmu[0] -kmu[1]*pbmu[1];
      double kpk= kmu[0]*pkmu[0] -kmu[1]*pkmu[1];
    
      
      double A1 = (kq*(mh*mh -mp*mp)/(mv*mv*q2)) - (kpb+kpk)/(mv*mv) - 1 -(mh*mh - mp*mp)/q2;
      double A2 = (mh*mh- mp*mp)/(q2) - (kq*(mh*mh -mp*mp)/(q2*mv*mv));

  
      double x= (Ev+Ep-E)*ad;
      double x0 = (Ev+Ep)*ad;
      double xrev=(Ev+Ep+E)*ad;
      double xrev2=(Ev+Ep+2*E)*ad;

      auto func_re = [&ad, &eps, &G](double x) {
	double eps_p = eps + G/2.0;
	//double m= (exp(x) + exp(-3*x))/(1-exp(-2*x));
	//return 2*ad*cos(eps*ad)/(m - cos(2*eps*ad)/sinh(x));
	return (x/ad)/((x/ad)*(x/ad) + eps_p*eps_p);
      };

      auto func_im = [&ad, &eps, &G](double x) {
	double eps_p = eps + G/2.0;
	//double m= (exp(x) + exp(-3*x))/(1+exp(-2*x));
	//return 2*ad*sin(eps*ad)/(m - cos(2*eps*ad)/cosh(x));

	return eps_p/( (x/ad)*(x/ad) + eps_p*eps_p);
      };
       
          
      double im_fact = func_im(x) + 3.0*func_im(xrev) -3.0*func_im(x0) -func_im(xrev2); //   3.0*im_fact_REV - 3.0*im_fact0 -im_fact_REV2 ;
      double re_fact = func_re(x) + 3.0*func_re(xrev) -3.0*func_re(x0) -func_re(xrev2); // re_fact += 3.0*re_fact_REV - 3.0*re_fact0 -re_fact_REV2;


     
      
    
      if(Use_conserved_current) sum +=  -1.0*(pow(ad,3)/pow(Zv,2))*(2.0)*(pow(fV*mv,2)/(2*Ev))*mom*(A1*F_plus+A2*F_zero)*((reim==1)?re_fact:im_fact);
      
      else sum += -1.0*(2.0*pow(ad,3)/pow(Zv,3))*(pow(fV*mv,2)/(2*Ev))*mom*(A1*F_plus+A2*F_zero)*((reim==1)?re_fact:im_fact);
      
    }

   

    //perturbative part

    double th=3.8;

    int which_kernel=1;

     auto pert = [&](double x) { // x should be energy in physical units

	Vfloat Qs({x,mom});

	double Qq= Qs[0]*qmu[0] - Qs[1]*qmu[1];
	double Qpb= Qs[0]*pbmu[0] -Qs[1]*pbmu[1];
	double Qpk= Qs[0]*pkmu[0] -Qs[1]*pkmu[1];

	double s=x*x -mom*mom;

	double th2_r = th*th - mom*mom;

	if(s < th2_r) return 0.0;

	double PI_scal= (3.0/(12.0*M_PI*M_PI))*sqrt( 1 - th2_r/s)*( 1 + 0.5*th2_r/s) ;
	

	double F1= s*( pkmu[1] - (mh*mh - mp*mp)*qmu[1]/(q2)) -mom*( Qpb + Qpk - (mh*mh-mp*mp)*Qq/q2 );
	double F2= s*(  (mh*mh -mp*mp)*qmu[1]/q2) - mom*( (mh*mh-mp*mp)*Qq/q2 );

	auto kk = [&](double y) {

	  //double RE_A= (exp(y) + exp(-3*y))/(1-exp(-2*y));
	  //double RE_B = 2*ad*cos(eps*ad)/(RE_A - cos(2*eps*ad)/sinh(y));

	  //double IM_A= (exp(y) + exp(-3*y))/(1+exp(-2*y));
	  //double IM_B= 2*ad*sin(eps*ad)/(IM_A - cos(2*eps*ad)/cosh(y));

	  double RE_B= (y/ad)/((y/ad)*(y/ad) + (eps*eps));
	  double IM_B = eps/( (y/ad)*(y/ad) + (eps*eps));

	  return (reim==1)?RE_B:IM_B;
	};

	

	double ker= kk( (x+Ep-E)*ad )  + 3.0*kk( (x+Ep+E)*ad) - 3.0*kk( (x+Ep)*ad) -kk( (x+Ep+2*E)*ad);

	if(which_kernel==0) ker= 1.0;
	if(which_kernel==-1) ker = 3.0*kk( (x+Ep+E)*ad) - 3.0*kk( (x+Ep)*ad) -kk( (x+Ep+2*E)*ad);

	return -(F1*F_plus + F2*F_zero)*PI_scal*ker*(pow(ad,3)/pow(Zv,2))*2.0;

      };


     if( fabs(eps) < 1e-10) {

       if(reim==0) {
	 which_kernel=0;
	 cout<<"IM eps=0: "<<M_PI*pert(E-Ep)<<endl;
	 return sum + M_PI*pert(E-Ep);
       }
       else {
	 which_kernel=0;
	 gsl_function_pp<decltype(pert)> Fp(pert);
	 gsl_integration_workspace * w = gsl_integration_workspace_alloc (30000);
	 gsl_function *G = static_cast<gsl_function*>(&Fp);
	 double pert_res, err_pert_res;
	 //perturbative tail starts at 4 GeV in \bar{c}c
	 gsl_integration_qawc(G, th, 100/ad, E-Ep, 0.0, 5e-5, 30000, w, &pert_res, &err_pert_res);
	 gsl_integration_workspace_free (w);
	 if(err_pert_res/fabs(pert_res) > 1e-4) cout<<"Warning precision achieved  > 1e-3, prec: "<<err_pert_res/fabs(pert_res)<<endl;
	 which_kernel=-1;
	 gsl_function_pp<decltype(pert)> Fp_m1(pert);
	 gsl_integration_workspace * w_m1 = gsl_integration_workspace_alloc (30000);
	 gsl_function *G_m1 = static_cast<gsl_function*>(&Fp_m1);
	 double pert_res_m1, err_pert_res_m1;
	 gsl_integration_qags(G_m1, th, 100/ad, 0.0, 5e-5, 30000, w_m1, &pert_res_m1, &err_pert_res_m1);
	 gsl_integration_workspace_free (w_m1);
	 if(err_pert_res_m1/fabs(pert_res_m1) > 1e-4) cout<<"Warning precision achieved  > 1e-3, prec: "<<err_pert_res_m1/fabs(pert_res_m1)<<endl;

	 return sum + (pert_res + pert_res_m1);

       }

       which_kernel=1;

     }
    
     gsl_function_pp<decltype(pert)> Fp(pert);
     gsl_integration_workspace * w = gsl_integration_workspace_alloc (30000);
     gsl_function *G = static_cast<gsl_function*>(&Fp);
     double pert_res, err_pert_res;
     //perturbative tail starts at 4 GeV in \bar{c}c
     gsl_integration_qags(G, th, 100/ad, 0.0, 5e-5, 30000, w, &pert_res, &err_pert_res);
     gsl_integration_workspace_free (w);
     if(err_pert_res/fabs(pert_res) > 1e-4) cout<<"Warning precision achieved  > 1e-3, prec: "<<err_pert_res/fabs(pert_res)<<endl; 
     return sum + pert_res;
  };
  
  bool L= true;
  //if(!Use_conserved_current) L= true;

  
   auto mult_func_IM = [&L](double erg, double s, double E0) {

     double z= s/fabs((fabs(erg) - E0));
     if(erg <0 ) erg = 0.0;

    assert(L==true);

    if( (erg < 1e-10 )  ) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return L?0.005:0.1;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.005:0.01;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.005:0.01;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.005:0.005;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.005:0.0005;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.005:0.0005;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.005:0.0005;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.005:0.0005;
      else crash("z not found");
    }
    else if( (erg < 3.55 +eps(5) )  && ( erg > 3.55 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return L?10.0:0.1;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:0.01;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.0:0.01;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.5:0.005;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.05:0.0005;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.05:0.0005;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.05:0.0005;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.05:0.0005;
      else crash("z not found");
    }
    else if( (erg < 3.6 +eps(5) )  && ( erg > 3.6 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?2.5:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?10.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.4:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.4:-100;
      else crash("z not found");
    }
    else if( (erg < 3.64 +eps(5) )  && ( erg > 3.64 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?10.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?4.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?2.0:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?1.0:-100;
      else crash("z not found");
    }
    else if( (erg < 3.66 +eps(5) )  && ( erg > 3.66 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?10.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?10.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?10.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?10.0:-100;
      else crash("z not found");
    }
    else if( (erg < 3.75 +eps(5) )  && ( erg > 3.75 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?70.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?5.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?5.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?7.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?7.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?7.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?5.0:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?1.0:-100;
      else crash("z not found");
    }
    else if( (erg < 3.8 +eps(5) )  && ( erg > 3.8 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?3.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?5.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?5.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?6.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?6.0:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?6.0:-100;
      else crash("z not found");
    }
    else if( (erg < 3.9 +eps(5) )  && ( erg > 3.9 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?5.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?10.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?10.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?10.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?6.0:-100;
      else crash("z not found");
    }
    else if( (erg < 4.0 +eps(5) )  && ( erg > 4.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?7.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?7.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?7.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?7.0:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?7.0:-100;
      else crash("z not found");
    }
    else if( (erg < 4.2 +eps(5) )  && ( erg > 4.2 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?15.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?25.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?20.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?20.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?20.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?7.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.1:-100;
      else crash("z not found");
    }
    else if( (erg < 4.4 +eps(5) )  && ( erg > 4.4 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?100.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?100.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?100.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?7.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?3.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.15:-100;
      else crash("z not found");
    }
    else if( (erg < 4.5 +eps(5) )  && ( erg > 4.5 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?20.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?5.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?2.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.07:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:-100;
      else crash("z not found");
    }
    else if( (erg < 4.75 +eps(5) )  && ( erg > 4.75 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?33.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.7:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.3:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.03:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:-100;
      else crash("z not found");
    }
    else if( (erg < 5.0 +eps(5) )  && ( erg > 5.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?23.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?3.6:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.08:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.3:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.09:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.02:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.007:-100;
      else crash("z not found");
    }
    else if( (erg < 5.28 +eps(5) )  && ( erg > 5.28 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.2:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.6:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?2.4:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.08:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.07:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.02:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:0.0001;
      else crash("z not found");
    }
    else crash("Plateaux not found");
    
    return 1e3;
  };



   auto mult_func_RE = [&L](double erg, double s, double E0) {

     
    double z= s/fabs((fabs(erg) - E0));
    if(erg <0 ) erg = 0.0;

    assert(L==true);

    if( erg  < 1e-10 ) {

      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.00001:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.00001:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.00001:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.00001:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.00001:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.00001:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.00001:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.00001:-100;
      else crash("z not found");

    }
    else if( (erg < 3.55 +eps(5) )  && ( erg > 3.55 -eps(5))) {  //here
      
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.0001:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.0001:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0001:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.0001:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0001:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.0001:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0001:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0001:-100;
      else crash("z not found");
    }
    else if( (erg < 3.6 +eps(5) )  && ( erg > 3.6 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.0003:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.0005:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0005:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.0005:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0002:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.0004:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0003:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0001:-100;
      else crash("z not found");
    }
    else if( (erg < 3.64 +eps(5) )  && ( erg > 3.64 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.01:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.006:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.002:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.007:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.005:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
      else crash("z not found");
    }
    else if( (erg < 3.66 +eps(5) )  && ( erg > 3.66 -eps(5))) { 
      
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.02:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.06:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.03:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.02:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.05:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.05:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.04:-100;
      else crash("z not found");
    }
    else if( (erg < 3.75 +eps(5) )  && ( erg > 3.75 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.2:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.5:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.5:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.7:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.9:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.6:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.1:-100;
      else crash("z not found");
    }
    else if( (erg < 3.80 +eps(5) )  && ( erg > 3.80 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.01:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?1.5:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?1.1:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.02:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.01:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.04:-100;
      else crash("z not found");
    }
    else if( (erg < 3.9 +eps(5) )  && ( erg > 3.9 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.4:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.09:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.7:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.5:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.3:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.4:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.2:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.2:-100;
      else crash("z not found");
    }
    else if( (erg < 4.0 +eps(5) )  && ( erg > 4.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.2:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.5:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?1.5:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.2:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.07:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.03:-100;
      else crash("z not found");
    }
    else if( (erg < 4.2 +eps(5) )  && ( erg > 4.2 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.12:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.4:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.05:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:-100;
      else crash("z not found");
    }
    else if( (erg < 4.4 +eps(5) )  && ( erg > 4.4 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.8:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.2:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.2:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.1:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.02:-100;
      else crash("z not found");
    }
    else if( (erg < 4.5 +eps(5) )  && ( erg > 4.5 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?3.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.5:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.3:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.07:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.06:-100;
      else crash("z not found");
    }
    else if( (erg < 4.75 +eps(5) )  && ( erg > 4.75 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?4.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.02:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.15:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.3:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.2:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.05:-100;
      else crash("z not found");
    }
    else if( (erg < 5.0 +eps(5) )  && ( erg > 5.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?4.0:-100; 
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.24:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.125:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.8:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?1.5:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?1.0:-100;
      else crash("z not found");
    }
    else if( (erg < 5.28 +eps(5) )  && ( erg > 5.28 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.01:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.8:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.2:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.4:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.5:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?5.0:-100;
      else crash("z not found");
    }
    else crash("Plateaux not found");
    
    return 1e3;
  };



   //################################   O2   B_s -> eta_ss' ell+ ell- ###################################

   auto mult_func_IM_O2 = [&L](double erg, double s, double E0) {

     double z= s/fabs((fabs(erg) - E0));
     if(erg <0 ) erg = 0.0;
     
     assert(L==true);
     
     if( (erg < 1e-10 )  ) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return L?0.0005:0.1;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.0005:0.01;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0005:0.01;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.0005:0.005;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0005:0.0005;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.0005:0.0005;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0005:0.0005;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0005:0.0005;
       else crash("z not found");
     }
     else if( (erg < 3.55 +eps(5) )  && ( erg > 3.55 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return L?10.0:0.1;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:0.01;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.001:0.01;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.001:0.005;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.001:0.0005;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:0.0005;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.001:0.0005;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:0.0005;
       else crash("z not found");
     }
     else if( (erg < 3.6 +eps(5) )  && ( erg > 3.6 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.01:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
       else crash("z not found");
     }
     else if( (erg < 3.64 +eps(5) )  && ( erg > 3.64 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.006:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
       else crash("z not found");
     }
     else if( (erg < 3.66 +eps(5) )  && ( erg > 3.66 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?50.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?50.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.01:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.01:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.01:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.03:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.01:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:-100;
       else crash("z not found");
     }
     else if( (erg < 3.75 +eps(5) )  && ( erg > 3.75 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?70.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?70.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.01:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.01:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.01:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.01:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.01:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:-100;
       else crash("z not found");
     }
     else if( (erg < 3.8 +eps(5) )  && ( erg > 3.8 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0005:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.005:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0005:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.01:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0135:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0037:-100;
       else crash("z not found");
     }
     else if( (erg < 3.9 +eps(5) )  && ( erg > 3.9 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?50.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?50.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0046:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
       else crash("z not found");
     }
     else if( (erg < 4.0 +eps(5) )  && ( erg > 4.0 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.05:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.05:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.05:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.05:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.003:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
       else crash("z not found");
     }
     else if( (erg < 4.2 +eps(5) )  && ( erg > 4.2 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?15.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?25.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0001:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.0001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0001:-100;
       else crash("z not found");
     }
     else if( (erg < 4.4 +eps(5) )  && ( erg > 4.4 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?50.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?100.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.01:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.01:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.01:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.005:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.003:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0001:-100;
       else crash("z not found");
     }
     else if( (erg < 4.5 +eps(5) )  && ( erg > 4.5 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?20.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.02:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.003:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0003:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.00001:-100;
      else crash("z not found");
     }
     else if( (erg < 4.75 +eps(5) )  && ( erg > 4.75 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?33.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?7.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.00004:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.00004:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.00004:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.0001:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0001:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0001:-100;
      else crash("z not found");
     }
     else if( (erg < 5.0 +eps(5) )  && ( erg > 5.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?23.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0003:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.0001:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0003:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.0001:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.00001:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.00001:-100;
      else crash("z not found");
     }
     else if( (erg < 5.28 +eps(5) )  && ( erg > 5.28 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.2:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.004:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.0001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0001:0.0001;
       else crash("z not found");
     }
     else crash("Plateaux not found");
     
    return 1e3;
   };
   
   
   
   auto mult_func_RE_O2 = [&L](double erg, double s, double E0) {
     
     
     double z= s/fabs((fabs(erg) - E0));
     if(erg <0 ) erg = 0.0;
     
     assert(L==true);
     
     if( erg  < 1e-10 ) {
       
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.000001:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.000001:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.000001:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.000001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.000001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.000001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.000001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.000001:-100;
       else crash("z not found");
       
     }
     else if( (erg < 3.55 +eps(5) )  && ( erg > 3.55 -eps(5))) {  //here
       
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.0001:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.0001:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.00001:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.00001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.00001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.00001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.00001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.00001:-100;
       else crash("z not found");
     }
     else if( (erg < 3.6 +eps(5) )  && ( erg > 3.6 -eps(5))) { 
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.0003:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.0005:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.00001:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.00001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.00001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.00001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.00001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.00001:-100;
       else crash("z not found");
     }
     else if( (erg < 3.64 +eps(5) )  && ( erg > 3.64 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.001:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.001:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.00001:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.00001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.00001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.00001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.00001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.00001:-100;
       else crash("z not found");
     }
     else if( (erg < 3.66 +eps(5) )  && ( erg > 3.66 -eps(5))) { 
       
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.02:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.06:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.00033:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.00033:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.00033:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.00033:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.00033:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.00033:-100;
       else crash("z not found");
     }
     else if( (erg < 3.75 +eps(5) )  && ( erg > 3.75 -eps(5))) { 
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?2.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?2.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.01:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.01:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.01:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.01:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.006:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
       else crash("z not found");
     }
     else if( (erg < 3.8 +eps(5) )  && ( erg > 3.8 -eps(5))) { 
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.001:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.0007:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0004:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.00027:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0001:-100;
       else crash("z not found");
     }
     else if( (erg < 3.9 +eps(5) )  && ( erg > 3.9 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.4:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.9:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.004:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.0007:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.00047:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.0003:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0003:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0001:-100;
      else crash("z not found");
     }
     else if( (erg < 4.0 +eps(5) )  && ( erg > 4.0 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?2.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?2.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0008:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.0003:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.00004:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.00005:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.00001:-100;
       else crash("z not found");
     }
     else if( (erg < 4.2 +eps(5) )  && ( erg > 4.2 -eps(5))) { 
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?2.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.001:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.001:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
       else crash("z not found");
     }
     else if( (erg < 4.4 +eps(5) )  && ( erg > 4.4 -eps(5))) { 
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.000008:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.000008:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.000008:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.0001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0004:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
       else crash("z not found");
     }
     else if( (erg < 4.5 +eps(5) )  && ( erg > 4.5 -eps(5))) { 
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?3.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0001:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.0001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.00015:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.0003:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
       else crash("z not found");
     }
     else if( (erg < 4.75 +eps(5) )  && ( erg > 4.75 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?4.0:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?2.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0006:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0005:-100;
       else crash("z not found");
     }
     else if( (erg < 5.0 +eps(5) )  && ( erg > 5.0 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?4.0:-100; 
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?2.4:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0001:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
       else crash("z not found");
     }
     else if( (erg < 5.28 +eps(5) )  && ( erg > 5.28 -eps(5))) {
       if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.1:-100;
       else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
       else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0001:-100;
       else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.0003:-100;
       else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:-100;
       else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.001:-100;
       else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
       else crash("z not found");
     }
     else crash("Plateaux not found");
     
     return 1e3;
   };
   
   
   //#####################################################################################################



   //################################   O1   B -> K ell+ ell- ###################################

    
   auto mult_func_IM_O1_HSL = [&L](double erg, double s, double E0) {

     double z= s/fabs((fabs(erg) - E0));
     if(erg <0 ) erg = 0.0;

    assert(L==true);

    if( (erg < 1e-10 )  ) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return L?0.005:0.1;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.005:0.01;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.005:0.01;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.005:0.005;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.005:0.0005;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.005:0.0005;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.005:0.0005;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.005:0.0005;
      else crash("z not found");
    }
    else if( (erg < 3.55 +eps(5) )  && ( erg > 3.55 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return L?10.0:0.1;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:0.01;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.0:0.01;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.5:0.005;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.05:0.0005;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.05:0.0005;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.05:0.0005;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.05:0.0005;
      else crash("z not found");
    }
    else if( (erg < 3.6 +eps(5) )  && ( erg > 3.6 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?2.5:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?10.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.4:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.4:-100;
      else crash("z not found");
    }
    else if( (erg < 3.64 +eps(5) )  && ( erg > 3.64 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?1.0:-100;
      else crash("z not found");
    }
    else if( (erg < 3.66 +eps(5) )  && ( erg > 3.66 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?3.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?3.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?1.0:-100;
      else crash("z not found");
    }
    else if( (erg < 3.75 +eps(5) )  && ( erg > 3.75 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?70.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?70.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?70.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?7.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?7.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?3.0:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?1.0:-100;
      else crash("z not found");
    }
    else if( (erg < 3.8 +eps(5) )  && ( erg > 3.8 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.6:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.6:-100;
      else crash("z not found");
    }
    else if( (erg < 3.9 +eps(5) )  && ( erg > 3.9 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?3.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.6:-100;
      else crash("z not found");
    }
    else if( (erg < 4.0 +eps(5) )  && ( erg > 4.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.8:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.7:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.2:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.7:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.35:-100;
      else crash("z not found");
    }
    else if( (erg < 4.2 +eps(5) )  && ( erg > 4.2 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?15.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?25.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?20.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?10.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?10.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.7:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.1:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:-100;
      else crash("z not found");
    }
    else if( (erg < 4.4 +eps(5) )  && ( erg > 4.4 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?100.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?100.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.7:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.03:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.05:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.015:-100;
      else crash("z not found");
    }
    else if( (erg < 4.5 +eps(5) )  && ( erg > 4.5 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?20.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?5.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.7:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?4.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.07:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:-100;
      else crash("z not found");
    }
    else if( (erg < 4.75 +eps(5) )  && ( erg > 4.75 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?33.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?7.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?5.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.03:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:-100;
      else crash("z not found");
    }
    else if( (erg < 5.0 +eps(5) )  && ( erg > 5.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?23.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?3.6:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.08:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.3:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.09:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.02:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.007:-100;
      else crash("z not found");
    }
    else if( (erg < 5.28 +eps(5) )  && ( erg > 5.28 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?20.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?6.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?2.4:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.08:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.07:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.02:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:0.0001;
      else crash("z not found");
    }
    else crash("Plateaux not found");
    
    return 1e3;
  };



   auto mult_func_RE_O1_HSL = [&L](double erg, double s, double E0) {

     
    double z= s/fabs((fabs(erg) - E0));
    if(erg <0 ) erg = 0.0;

    assert(L==true);

    if( erg  < 1e-10 ) {

      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.00001:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.00001:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.00001:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.00001:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.00001:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.00001:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.00001:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.00001:-100;
      else crash("z not found");

    }
    else if( (erg < 3.55 +eps(5) )  && ( erg > 3.55 -eps(5))) {  //here
      
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.0001:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.0001:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0001:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.0001:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0001:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.0001:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0001:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0001:-100;
      else crash("z not found");
    }
    else if( (erg < 3.6 +eps(5) )  && ( erg > 3.6 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.0003:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.0005:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0005:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.0005:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0002:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.0004:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0003:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0001:-100;
      else crash("z not found");
    }
    else if( (erg < 3.64 +eps(5) )  && ( erg > 3.64 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.01:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.006:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?2.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.7:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.3:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:-100;
      else crash("z not found");
    }
    else if( (erg < 3.66 +eps(5) )  && ( erg > 3.66 -eps(5))) { 
      
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.02:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.06:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.3:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.2:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.05:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.05:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.04:-100;
      else crash("z not found");
    }
    else if( (erg < 3.75 +eps(5) )  && ( erg > 3.75 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.8:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.09:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.06:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.1:-100;
      else crash("z not found");
    }
    else if( (erg < 3.8 +eps(5) )  && ( erg > 3.8 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.15:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.11:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.08:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.1:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.04:-100;
      else crash("z not found");
    }
    else if( (erg < 3.9 +eps(5) )  && ( erg > 3.9 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.4:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.9:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.7:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.5:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.3:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.4:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.2:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.05:-100;
      else crash("z not found");
    }
    else if( (erg < 4.0 +eps(5) )  && ( erg > 4.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.5:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?1.5:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.6:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.05:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:-100;
      else crash("z not found");
    }
    else if( (erg < 4.2 +eps(5) )  && ( erg > 4.2 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.2:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.005:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.05:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:-100;
      else crash("z not found");
    }
    else if( (erg < 4.4 +eps(5) )  && ( erg > 4.4 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.8:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.1:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.02:-100;
      else crash("z not found");
    }
    else if( (erg < 4.5 +eps(5) )  && ( erg > 4.5 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?3.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.2:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.3:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.07:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.006:-100;
      else crash("z not found");
    }
    else if( (erg < 4.75 +eps(5) )  && ( erg > 4.75 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?4.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.5:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.03:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.002:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.05:-100;
      else crash("z not found");
    }
    else if( (erg < 5.0 +eps(5) )  && ( erg > 5.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?4.0:-100; 
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?2.4:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.25:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.008:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.015:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?1.0:-100;
      else crash("z not found");
    }
    else if( (erg < 5.28 +eps(5) )  && ( erg > 5.28 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?6.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?2.4:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.03:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.005:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.04:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.5:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?5.0:-100;
      else crash("z not found");
    }
    else crash("Plateaux not found");
    
    return 1e3;
   };
   
   //#####################################################################################################



   //################################   O2   B -> K ell+ ell- ###################################

    
   auto mult_func_IM_O2_HSL = [&L](double erg, double s, double E0) {

     double z= s/fabs((fabs(erg) - E0));
     if(erg <0 ) erg = 0.0;

    assert(L==true);

    if( (erg < 1e-10 )  ) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return L?0.005:0.1;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.005:0.01;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.005:0.01;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.005:0.005;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.005:0.0005;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.005:0.0005;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.005:0.0005;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.005:0.0005;
      else crash("z not found");
    }
    else if( (erg < 3.55 +eps(5) )  && ( erg > 3.55 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return L?10.0:0.1;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:0.01;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.0:0.01;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.05:0.005;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.005:0.0005;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.005:0.0005;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.005:0.0005;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.005:0.0005;
      else crash("z not found");
    }
    else if( (erg < 3.6 +eps(5) )  && ( erg > 3.6 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.25:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.04:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.04:-100;
      else crash("z not found");
    }
    else if( (erg < 3.64 +eps(5) )  && ( erg > 3.64 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.1:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.1:-100;
      else crash("z not found");
    }
    else if( (erg < 3.66 +eps(5) )  && ( erg > 3.66 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?3.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.3:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.1:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.1:-100;
      else crash("z not found");
    }
    else if( (erg < 3.75 +eps(5) )  && ( erg > 3.75 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?7.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?7.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?7.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?.7:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.7:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.7:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.3:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.1:-100;
      else crash("z not found");
    }
    else if( (erg < 3.8 +eps(5) )  && ( erg > 3.8 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.7:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.06:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.06:-100;
      else crash("z not found");
    }
    else if( (erg < 3.9 +eps(5) )  && ( erg > 3.9 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.3:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.1:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.06:-100;
      else crash("z not found");
    }
    else if( (erg < 4.0 +eps(5) )  && ( erg > 4.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.08:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.07:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.02:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.07:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.035:-100;
      else crash("z not found");
    }
    else if( (erg < 4.2 +eps(5) )  && ( erg > 4.2 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?15.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?25.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?20.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?1.0:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.07:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.01:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
      else crash("z not found");
    }
    else if( (erg < 4.4 +eps(5) )  && ( erg > 4.4 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?50.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?100.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?100.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.07:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.003:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.005:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0015:-100;
      else crash("z not found");
    }
    else if( (erg < 4.5 +eps(5) )  && ( erg > 4.5 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?20.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?5.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.07:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.4:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.007:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
      else crash("z not found");
    }
    else if( (erg < 4.75 +eps(5) )  && ( erg > 4.75 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?33.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?7.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?5.0:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.003:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
      else crash("z not found");
    }
    else if( (erg < 5.0 +eps(5) )  && ( erg > 5.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?23.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?3.6:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.008:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.03:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.009:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.002:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0007:-100;
      else crash("z not found");
    }
    else if( (erg < 5.28 +eps(5) )  && ( erg > 5.28 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?20.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?6.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?2.4:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.08:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.007:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.002:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:0.0001;
      else crash("z not found");
    }
    else crash("Plateaux not found");
    
    return 1e3;
  };



   auto mult_func_RE_O2_HSL = [&L](double erg, double s, double E0) {

     
    double z= s/fabs((fabs(erg) - E0));
    if(erg <0 ) erg = 0.0;

    assert(L==true);

    if( erg  < 1e-10 ) {

      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.00001:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.00001:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.00001:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.00001:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.00001:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.00001:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.00001:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.00001:-100;
      else crash("z not found");

    }
    else if( (erg < 3.55 +eps(5) )  && ( erg > 3.55 -eps(5))) {  //here
      
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.0001:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.001:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.001:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.001:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.001:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.001:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
      else crash("z not found");
    }
    else if( (erg < 3.6 +eps(5) )  && ( erg > 3.6 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.0003:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.005:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.005:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.005:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.002:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.004:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.003:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
      else crash("z not found");
    }
    else if( (erg < 3.64 +eps(5) )  && ( erg > 3.64 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.01:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.006:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.2:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.1:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.07:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.03:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
      else crash("z not found");
    }
    else if( (erg < 3.66 +eps(5) )  && ( erg > 3.66 -eps(5))) { 
      
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?0.02:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.06:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.03:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.02:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.005:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.005:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.004:-100;
      else crash("z not found");
    }
    else if( (erg < 3.75 +eps(5) )  && ( erg > 3.75 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.2:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.8:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.03:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.009:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.006:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.01:-100;
      else crash("z not found");
    }
    else if( (erg < 3.8 +eps(5) )  && ( erg > 3.8 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.1:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.015:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.011:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.008:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.01:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.004:-100;
      else crash("z not found");
    }
    else if( (erg < 3.9 +eps(5) )  && ( erg > 3.9 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.4:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?0.9:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.007:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.05:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.03:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.04:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.02:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.005:-100;
      else crash("z not found");
    }
    else if( (erg < 4.0 +eps(5) )  && ( erg > 4.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.015:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.15:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.06:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.005:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
      else crash("z not found");
    }
    else if( (erg < 4.2 +eps(5) )  && ( erg > 4.2 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.2:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.01:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0005:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.005:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.001:-100;
      else crash("z not found");
    }
    else if( (erg < 4.4 +eps(5) )  && ( erg > 4.4 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.008:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.01:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.002:-100;
      else crash("z not found");
    }
    else if( (erg < 4.5 +eps(5) )  && ( erg > 4.5 -eps(5))) { 
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?3.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?1.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.01:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.02:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.03:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.01:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.007:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.0006:-100;
      else crash("z not found");
    }
    else if( (erg < 4.75 +eps(5) )  && ( erg > 4.75 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?4.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?2.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.015:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.10:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.003:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0002:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.005:-100;
      else crash("z not found");
    }
    else if( (erg < 5.0 +eps(5) )  && ( erg > 5.0 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?4.0:-100; 
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?2.4:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0125:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.0008:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.001:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.001:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.0015:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.1:-100;
      else crash("z not found");
    }
    else if( (erg < 5.28 +eps(5) )  && ( erg > 5.28 -eps(5))) {
      if( (z < 0.3 + eps(5)) && (z > 0.3 - eps(5))) return  L?10.0:-100;
      else if( (z < 0.5 + eps(5)) && (z > 0.5 - eps(5))) return  L?6.0:-100;
      else if( (z < 0.6 + eps(5)) && (z > 0.6 - eps(5))) return  L?0.0024:-100;
      else if( (z < 0.8 + eps(5)) && (z > 0.8 - eps(5))) return  L?0.003:-100;
      else if( (z < 1.0 + eps(5)) && (z > 1.0 - eps(5))) return  L?0.0005:-100;
      else if( (z < 1.25 + eps(5)) && (z > 1.25 - eps(5))) return  L?0.004:-100;
      else if( (z < 1.5 + eps(5)) && (z > 1.5 - eps(5))) return  L?0.05:-100;
      else if( (z < 2.0 + eps(5)) && (z > 2.0 - eps(5))) return  L?0.50:-100;
      else crash("z not found");
    }
    else crash("Plateaux not found");
    
    return 1e3;
   };
   
   //#####################################################################################################






   
   

  

  omp_set_num_threads(4);

  //load correlators

  int T=128;

  string pingu_dir="../data/PINGU";

  if(!Use_conserved_current) pingu_dir="../data/PINGU_local";

  string out_dir= "PINGU";
  if(!Use_conserved_current) out_dir = "PINGU_local";


 
  
  
  //load correlators for B -> K ll

  vector<distr_t_list> C1_A, C1_V, C2_A, C2_V;
  vector<distr_t_list> CS1_A, CS1_V, CS2_A, CS2_V;
  
  vector<distr_t_list> C1_A_cov, C1_V_cov, C2_A_cov, C2_V_cov;
  vector<distr_t_list> CS1_A_cov, CS1_V_cov, CS2_A_cov, CS2_V_cov;

  for(int ts=0;ts<T;ts++) { C1_A.emplace_back(UseJack), C1_V.emplace_back(UseJack), C2_A.emplace_back(UseJack), C2_V.emplace_back(UseJack) ;}
  for(int ts=0;ts<T;ts++) { CS1_A.emplace_back(UseJack), CS1_V.emplace_back(UseJack), CS2_A.emplace_back(UseJack), CS2_V.emplace_back(UseJack) ;}

  for(int ts=0;ts<T;ts++) { C1_A_cov.emplace_back(UseJack), C1_V_cov.emplace_back(UseJack), C2_A_cov.emplace_back(UseJack), C2_V_cov.emplace_back(UseJack) ;}
  for(int ts=0;ts<T;ts++) { CS1_A_cov.emplace_back(UseJack), CS1_V_cov.emplace_back(UseJack), CS2_A_cov.emplace_back(UseJack), CS2_V_cov.emplace_back(UseJack) ;}

  for(int ts=0;ts<T;ts++) {
    for(int t=0;t<T; t++) { 
      C1_A[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSL/jackknife/HSL/C_JJ_O1_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 1, 4) );
      C1_V[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSL/jackknife/HSL/C_JJ_O1_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 3, 4) );
      C2_A[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSL/jackknife/HSL/C_JJ_O2_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 1, 4) );
      C2_V[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSL/jackknife/HSL/C_JJ_O2_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 3, 4) );
      
      CS1_A[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSS/jackknife/HSS/C_JJ_O1_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 1, 4) );
      CS1_V[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSS/jackknife/HSS/C_JJ_O1_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 3, 4) );
      CS2_A[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSS/jackknife/HSS/C_JJ_O2_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 1, 4) );
      CS2_V[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSS/jackknife/HSS/C_JJ_O2_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 3, 4) );

      C1_A_cov[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSL/jackknife/HSL/C_O1_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 1, 4) );
      C1_V_cov[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSL/jackknife/HSL/C_O1_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 3, 4) );
      C2_A_cov[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSL/jackknife/HSL/C_O2_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 1, 4) );
      C2_V_cov[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSL/jackknife/HSL/C_O2_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 3, 4) );
      
      CS1_A_cov[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSS/jackknife/HSS/C_O1_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 1, 4) );
      CS1_V_cov[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSS/jackknife/HSS/C_O1_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 3, 4) );
      CS2_A_cov[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSS/jackknife/HSS/C_O2_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 1, 4) );
      CS2_V_cov[ts].distr_list.emplace_back( UseJack, Read_From_File( pingu_dir+"/HSS/jackknife/HSS/C_O2_ts_"+to_string(ts)+"_t_"+to_string(t)+".dat", 3, 4) );

      
    }
  }

  cout<<"Correlator loaded"<<endl;



  //get final correlator

  distr_t_list C1_A_FIN(UseJack), C1_V_FIN(UseJack), C2_A_FIN(UseJack), C2_V_FIN(UseJack);

  distr_t_list CS1_A_FIN(UseJack), CS1_V_FIN(UseJack), CS2_A_FIN(UseJack), CS2_V_FIN(UseJack);


  distr_t_list C1_A_FIN_cov(UseJack), C1_V_FIN_cov(UseJack), C2_A_FIN_cov(UseJack), C2_V_FIN_cov(UseJack);

  distr_t_list CS1_A_FIN_cov(UseJack), CS1_V_FIN_cov(UseJack), CS2_A_FIN_cov(UseJack), CS2_V_FIN_cov(UseJack);

  int tsep=12;

  for(int t=0;t<=T/2;t++) {

    C1_A_FIN.distr_list.push_back(     C1_A[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );
    C1_V_FIN.distr_list.push_back(     C1_V[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );
    C2_A_FIN.distr_list.push_back(     C2_A[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );
    C2_V_FIN.distr_list.push_back(     C2_V[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );

    CS1_A_FIN.distr_list.push_back(     CS1_A[60].distr_list[t] );
    CS1_V_FIN.distr_list.push_back(     CS1_V[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );
    CS2_A_FIN.distr_list.push_back(     CS2_A[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );
    CS2_V_FIN.distr_list.push_back(     CS2_V[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );

    C1_A_FIN_cov.distr_list.push_back(     C1_A_cov[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );
    C1_V_FIN_cov.distr_list.push_back(     C1_V_cov[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );
    C2_A_FIN_cov.distr_list.push_back(     C2_A_cov[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );
    C2_V_FIN_cov.distr_list.push_back(     C2_V_cov[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );

    CS1_A_FIN_cov.distr_list.push_back(     CS1_A_cov[60].distr_list[t] );
    CS1_V_FIN_cov.distr_list.push_back(     CS1_V_cov[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );
    CS2_A_FIN_cov.distr_list.push_back(     CS2_A_cov[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );
    CS2_V_FIN_cov.distr_list.push_back(     CS2_V_cov[(t<int_t_O4)?(int_t_O4+tsep):(t+tsep)].distr_list[t] );
    
  }

  
  //symmetrize by hand
  for(int t=1;t<T/2;t++) {

    C1_A_FIN.distr_list.push_back( C1_A_FIN.distr_list[T/2 - t] );
    C1_V_FIN.distr_list.push_back( C1_V_FIN.distr_list[T/2 - t] );
    C2_A_FIN.distr_list.push_back( C2_A_FIN.distr_list[T/2 - t] );
    C2_V_FIN.distr_list.push_back( C2_V_FIN.distr_list[T/2 - t] );

    CS1_A_FIN.distr_list.push_back( CS1_A_FIN.distr_list[T/2 - t] );
    CS1_V_FIN.distr_list.push_back( CS1_V_FIN.distr_list[T/2 - t] );
    CS2_A_FIN.distr_list.push_back( CS2_A_FIN.distr_list[T/2 - t] );
    CS2_V_FIN.distr_list.push_back( CS2_V_FIN.distr_list[T/2 - t] );

    C1_A_FIN_cov.distr_list.push_back( C1_A_FIN_cov.distr_list[T/2 - t] );
    C1_V_FIN_cov.distr_list.push_back( C1_V_FIN_cov.distr_list[T/2 - t] );
    C2_A_FIN_cov.distr_list.push_back( C2_A_FIN_cov.distr_list[T/2 - t] );
    C2_V_FIN_cov.distr_list.push_back( C2_V_FIN_cov.distr_list[T/2 - t] );

    CS1_A_FIN_cov.distr_list.push_back( CS1_A_FIN_cov.distr_list[T/2 - t] );
    CS1_V_FIN_cov.distr_list.push_back( CS1_V_FIN_cov.distr_list[T/2 - t] );
    CS2_A_FIN_cov.distr_list.push_back( CS2_A_FIN_cov.distr_list[T/2 - t] );
    CS2_V_FIN_cov.distr_list.push_back( CS2_V_FIN_cov.distr_list[T/2 - t] );


  }

 

  //print correlators

  Print_To_File({}, {C1_A_FIN.ave(), C1_A_FIN.err(), C1_V_FIN.ave(), C1_V_FIN.err(), (C1_A_FIN+C1_V_FIN).ave(), (C1_A_FIN+C1_V_FIN).err() }, pingu_dir+"/HSL/C1_final", "", "");
  Print_To_File({}, {C2_A_FIN.ave(), C2_A_FIN.err(), C2_V_FIN.ave(), C2_V_FIN.err(), (C2_A_FIN+C2_V_FIN).ave(), (C2_A_FIN+C2_V_FIN).err() }, pingu_dir+"/HSL/C2_final", "", "");

  Print_To_File({}, {CS1_A_FIN.ave(), CS1_A_FIN.err(), CS1_V_FIN.ave(), CS1_V_FIN.err(),  (CS1_A_FIN+CS1_V_FIN).ave(), (CS1_A_FIN+CS1_V_FIN).err() }, pingu_dir+"/HSS/C1_final", "", "");
  Print_To_File({}, {CS2_A_FIN.ave(), CS2_A_FIN.err(), CS2_V_FIN.ave(), CS2_V_FIN.err(),  (CS2_A_FIN+CS2_V_FIN).ave(), (CS2_A_FIN+CS2_V_FIN).err() }, pingu_dir+"/HSS/C2_final", "", "");


  //define first time-ordering contribution


  
  //do HLT analysis

  //build HLT correlator 

  distr_t_list C1_HLT(UseJack), C2_HLT(UseJack), CS1_HLT(UseJack), CS2_HLT(UseJack);
  distr_t_list C1_HLT_cov(UseJack), C2_HLT_cov(UseJack), CS1_HLT_cov(UseJack), CS2_HLT_cov(UseJack);


  for(int t=0;t<=(T/2)-int_t_O4;t++) {
    C1_HLT.distr_list.push_back( C1_A_FIN.distr_list[t+int_t_O4 ] +  C1_V_FIN.distr_list[t+int_t_O4 ] );
    C2_HLT.distr_list.push_back( C2_A_FIN.distr_list[t+int_t_O4 ] +  C2_V_FIN.distr_list[t+int_t_O4 ] );
    
    CS1_HLT.distr_list.push_back( CS1_A_FIN.distr_list[t+int_t_O4 ] +  CS1_V_FIN.distr_list[t+int_t_O4 ] );
    CS2_HLT.distr_list.push_back( CS2_A_FIN.distr_list[t+int_t_O4 ] +  CS2_V_FIN.distr_list[t+int_t_O4 ] );
    
    C1_HLT_cov.distr_list.push_back( C1_A_FIN_cov.distr_list[t+int_t_O4 ]+ C1_V_FIN_cov.distr_list[t+int_t_O4 ] );
    C2_HLT_cov.distr_list.push_back( C2_A_FIN_cov.distr_list[t+int_t_O4 ]+ C2_V_FIN_cov.distr_list[t+int_t_O4 ] );
    
    CS1_HLT_cov.distr_list.push_back( CS1_A_FIN_cov.distr_list[t+int_t_O4 ]+  CS1_V_FIN_cov.distr_list[t+int_t_O4 ] );
    CS2_HLT_cov.distr_list.push_back( CS2_A_FIN_cov.distr_list[t+int_t_O4 ]+  CS2_V_FIN_cov.distr_list[t+int_t_O4 ] );
  }

  int T_HLT=C1_HLT.size();


  //print effective masses from correlators in the second TO

  distr_t_list C1_eff_mass(UseJack), CS1_eff_mass(UseJack);
  distr_t_list C2_eff_mass(UseJack), CS2_eff_mass(UseJack);
 
  
  for(int t=0; t < T_HLT-1; t++) {

    C1_eff_mass.distr_list.push_back( LOG_D( C1_HLT[t]/C1_HLT[t+1])) ;
    CS1_eff_mass.distr_list.push_back( LOG_D( CS1_HLT[t]/CS1_HLT[t+1]));

    C2_eff_mass.distr_list.push_back( LOG_D( C2_HLT[t]/C2_HLT[t+1])) ;
    CS2_eff_mass.distr_list.push_back( LOG_D( CS2_HLT[t]/CS2_HLT[t+1]));

  }
  
  
  CorrAnalysis Corr(UseJack, Njacks_JJ, 1000);
  Corr.Perform_Nt_t_average=0;
  Corr.Tmin=14; Corr.Tmax=30;
  distr_t mass= Corr.Fit_distr(CS1_eff_mass);
  distr_t res = Corr.Fit_distr( -1.0*CS1_HLT/EXPT_D(-1.0*mass, CS1_HLT.size()));
  double a_GeV =  0.07948*fm_to_inv_Gev;
  cout<<"mass(HSS): "<<(mass/a_GeV).ave()<<" "<<(mass/a_GeV).err()<<endl;



  

  Print_To_File({}, {C1_eff_mass.ave(), C1_eff_mass.err()}, "../data/"+out_dir+"/HSL/meff_C1_final", "", "");
  Print_To_File({}, {CS1_eff_mass.ave(), CS1_eff_mass.err()}, "../data/"+out_dir+"/HSS/meff_C1_final", "", "");

  Print_To_File({}, {C2_eff_mass.ave(), C2_eff_mass.err()}, "../data/"+out_dir+"/HSL/meff_C2_final", "", "");
  Print_To_File({}, {CS2_eff_mass.ave(), CS2_eff_mass.err()}, "../data/"+out_dir+"/HSS/meff_C2_final", "", "");

 
  //compute covariance matrix
  Vfloat Cov_1, Cov_2;
  Vfloat Cov_S1, Cov_S2;

  Vfloat Corr_1,Corr_2;
  Vfloat Corr_S1, Corr_S2;

  Vfloat TT, RR; 

  for(int tt=0; tt <T_HLT; tt++)
    for(int rr=0;rr <T_HLT;rr++) {
      TT.push_back(tt);
      RR.push_back(rr);

      Cov_1.push_back( C1_HLT_cov.distr_list[tt]%C1_HLT_cov.distr_list[rr]) ;
      Cov_2.push_back( C2_HLT_cov.distr_list[tt]%C2_HLT_cov.distr_list[rr]) ;
      
      Cov_S1.push_back( CS1_HLT_cov.distr_list[tt]%CS1_HLT_cov.distr_list[rr]) ;
      Cov_S2.push_back( CS2_HLT_cov.distr_list[tt]%CS2_HLT_cov.distr_list[rr]) ;
      
      Corr_1.push_back( Cov_1[ rr + T_HLT*tt]/(C1_HLT_cov.err(tt)*C1_HLT_cov.err(rr)));
      Corr_2.push_back( Cov_2[ rr + T_HLT*tt]/(C2_HLT_cov.err(tt)*C2_HLT_cov.err(rr)));
      
      Corr_S1.push_back( Cov_S1[ rr + T_HLT*tt]/(CS1_HLT_cov.err(tt)*CS1_HLT_cov.err(rr)));
      Corr_S2.push_back( Cov_S2[ rr + T_HLT*tt]/(CS2_HLT_cov.err(tt)*CS2_HLT_cov.err(rr)));
      
      Cov_1[rr+ T_HLT*tt] = Corr_1[rr +T_HLT*tt]*C1_HLT.err(tt)*C1_HLT.err(rr);
      Cov_2[rr+ T_HLT*tt] = Corr_2[rr +T_HLT*tt]*C2_HLT.err(tt)*C2_HLT.err(rr);
      
      Cov_S1[rr+ T_HLT*tt] = Corr_S1[rr +T_HLT*tt]*CS1_HLT.err(tt)*CS1_HLT.err(rr);
      Cov_S2[rr+ T_HLT*tt] = Corr_S2[rr +T_HLT*tt]*CS2_HLT.err(tt)*CS2_HLT.err(rr);
    }


  cout<<"covariance matrix computed"<<endl;

  //printing covariance matrix
  boost::filesystem::create_directory("../data/"+out_dir+"/spectral_reconstruction");
  boost::filesystem::create_directory("../data/"+out_dir+"/spectral_reconstruction/HSL");
  boost::filesystem::create_directory("../data/"+out_dir+"/spectral_reconstruction/HSS");
  
  boost::filesystem::create_directory("../data/"+out_dir+"/spectral_reconstruction/HSL/Erg");
  boost::filesystem::create_directory("../data/"+out_dir+"/spectral_reconstruction/HSS/Erg");
  
  boost::filesystem::create_directory("../data/"+out_dir+"/spectral_reconstruction/HSL/eps");
  boost::filesystem::create_directory("../data/"+out_dir+"/spectral_reconstruction/HSS/eps");
  
  boost::filesystem::create_directory("../data/"+out_dir+"/spectral_reconstruction/HSL/covariance");
  boost::filesystem::create_directory("../data/"+out_dir+"/spectral_reconstruction/HSS/covariance");
  

  Print_To_File({},{TT, RR, Cov_1, Corr_1 }, "../data/"+out_dir+"/spectral_reconstruction/HSL/covariance/Corr_1", "", "");
  Print_To_File({},{TT, RR, Cov_2, Corr_2 }, "../data/"+out_dir+"/spectral_reconstruction/HSL/covariance/Corr_2", "", "");
  Print_To_File({},{TT, RR, Cov_S1, Corr_S1 }, "../data/"+out_dir+"/spectral_reconstruction/HSS/covariance/Corr_1", "", "");
  Print_To_File({},{TT, RR, Cov_S2, Corr_S2 }, "../data/"+out_dir+"/spectral_reconstruction/HSS/covariance/Corr_2", "", "");


  //perform HLT analysis

  Vfloat Ergs({3.55,3.6,3.64,3.66,3.75,3.8,3.9,4.0, 4.2, 4.40, 4.50, 4.75, 5.0, 5.28});
  Vfloat Ergs_HSL({3.35,3.42,3.57,3.66,3.75,3.8,3.9,4.0, 4.2, 4.40, 4.50, 4.75, 5.0, 5.28    });
  Vfloat Ergs_NZ({3.55,3.6,3.64,3.66,3.75,3.8,3.9,4.0, 4.2, 4.40, 4.50, 4.75, 5.0, 5.28});
  Vfloat Zetas({0.3,0.5,0.6, 0.8, 1.0, 1.25, 1.5,2});


  int id_475= Ergs.size()-1;
  
  //Print model penguin
  Vfloat Zetas_mod({0.0});
  for(auto & z:Zetas) Zetas_mod.push_back(z);


  
  Vfloat Emod;
  int Nsteps=1500;


  VVfloat MOD_HSS_RE(Zetas_mod.size()), MOD_HSS_IM(Zetas_mod.size()), MOD_HSL_RE(Zetas_mod.size()), MOD_HSL_IM(Zetas_mod.size());


  Vfloat alphas_mod;
  for(int i=0;i<200;i++) alphas_mod.push_back( i*0.1) ;

  Vfloat MOD_HSS_528_RE, MOD_HSS_528_IM;

  for(int a=0;a<(signed)alphas_mod.size();a++) {

    double alps=alphas_mod[a];

    double s= fabs(alps*(5.28-3.7+0.150));
    MOD_HSS_528_RE.push_back( model_PENGUIN(5.28,s,"HSS",1));
    MOD_HSS_528_IM.push_back( model_PENGUIN(5.28,s,"HSS",0));

  }

  Print_To_File({}, {alphas_mod, MOD_HSS_528_RE, MOD_HSS_528_IM }, "../data/PINGU/spectral_reconstruction/fit/mod_528","","");
  
  
  for(int istep=0;istep<Nsteps;istep++) Emod.push_back( 3.4 + (5.28-3.4)*(1.0*istep)/(Nsteps-1.0));
  for(int is=0;is<Nsteps;is++) {
    //loop over Zetas

    double Erg= Emod[is];

    cout<<"Computing model n/N= "<<(is+1)*1.0/(1.0*Nsteps)<<endl;

    for(int iz=0;iz<(signed)Zetas_mod.size();iz++) {
      
    
      double s_HSS = fabs(Zetas_mod[iz]*(fabs(Erg-3.7)+0.150));
      double s_HSL = fabs(Zetas_mod[iz]*(fabs(Erg-3.5)+0.150));
      //cout<<"NOW HSS"<<endl;
      MOD_HSS_RE[iz].push_back( model_PENGUIN(Erg,s_HSS,"HSS",1 ));
      MOD_HSS_IM[iz].push_back( model_PENGUIN(Erg,s_HSS,"HSS",0 ));
      //cout<<"NOW HSL"<<endl;
      MOD_HSL_RE[iz].push_back( model_PENGUIN(Erg,s_HSL,"HSL",1 ));
      MOD_HSL_IM[iz].push_back( model_PENGUIN(Erg,s_HSL,"HSL",0 ));

      
      

      
    }
  }
  
  

  
  auto Ag_func = [](double erg, double s, double E0) {

    if(erg < E0) return 1e-4;

    if( (erg-E0)/s < 1) return 1e-3;
    else if( (erg-E0)/s < 2) return 1e-2;
    else if( (erg-E0)/s < 3)  return 3e-2;
    else if( (erg-E0)/s < 4) return 6e-2;
    else return 1e-1;

  };

  
  //define lambda functions

   auto K_RE= [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {
    
     PrecFloat x= (E-m);
	       
     PrecFloat cosh_ov_sinh_half= (exp(x) + exp(-3*x))/(1-exp(-2*x)); 
	       
     return  2*cos(s)/(cosh_ov_sinh_half - cos(2*s)/sinh(x));
	       
   };
   
   auto K_IM = [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {
     
     PrecFloat x= (E-m);
     
     PrecFloat cosh_ov_cosh_half= (exp(x) + exp(-3*x))/(1+exp(-2*x)); 
     
     return 2*sin(s)/(cosh_ov_cosh_half - cos(2*s)/cosh(x));
     
   };


   auto K_RE_0= [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {
    
     PrecFloat x= (E+m);
     PrecFloat x0= E;
     PrecFloat x2= (E+2*m);
     
	       
     PrecFloat cosh_ov_sinh_half= (exp(x) + exp(-3*x))/(1-exp(-2*x));
     PrecFloat cosh_ov_sinh_half_0= (exp(x0) + exp(-3*x0))/(1-exp(-2*x0));
     PrecFloat cosh_ov_sinh_half_2= (exp(x2) + exp(-3*x2))/(1-exp(-2*x2));
	       
     return  2*3*cos(s)/(cosh_ov_sinh_half - cos(2*s)/sinh(x)) +  -2.0*3*cos(s)/(cosh_ov_sinh_half_0 - cos(2*s)/sinh(x0))  -  2*cos(s)/(cosh_ov_sinh_half_2 - cos(2*s)/sinh(x2)) ;
	       
   };
   
   auto K_IM_0 = [](const PrecFloat &E, const PrecFloat &m, const PrecFloat &s, const PrecFloat &E0, int ijack) -> PrecFloat {
     
     PrecFloat x= (E+m);
     PrecFloat x0= E;
     PrecFloat x2= (E+2*m);
     
     PrecFloat cosh_ov_cosh_half= (exp(x) + exp(-3*x))/(1+exp(-2*x));
     PrecFloat cosh_ov_cosh_half_0= (exp(x0) + exp(-3*x0))/(1+exp(-2*x0));
     PrecFloat cosh_ov_cosh_half_2= (exp(x2) + exp(-3*x2))/(1+exp(-2*x2));
     
     return 2*3*sin(s)/(cosh_ov_cosh_half - cos(2*s)/cosh(x)) -  2.0*3*sin(s)/(cosh_ov_cosh_half_0 - cos(2*s)/cosh(x0))   - 2.0*sin(s)/(cosh_ov_cosh_half_2 - cos(2*s)/cosh(x2))  ;
     
   };

   
  


   vector<complex_distr_t_list> H1, H2;
   vector<complex_distr_t_list> HS1, HS2;

   for(auto & s: Zetas) {
     H1.emplace_back( UseJack, Ergs.size());
     H2.emplace_back( UseJack, Ergs.size());
    
     HS1.emplace_back( UseJack, Ergs.size());
     HS2.emplace_back( UseJack, Ergs.size());
   }

   if(! Skip_HLT) {

#pragma omp parallel for schedule(dynamic)
     for(int is=0;is < (signed)Zetas.size(); is++ ) {

    
    
     

        
     for(int eg=0; eg < (signed)Ergs.size(); eg++) {

       
       double Emin=3.5; //GeV
       double smin= 0.150*a_GeV;
       double aE0 = Emin*a_GeV;
       double erg= Ergs[eg]*a_GeV;
       double erg_HSL= Ergs_HSL[eg]*a_GeV;
       double sigma= Zetas[is]*fabs(( fabs(erg_HSL - aE0)+smin));
       double sigma_GeV= sigma/a_GeV;
       
       

       double syst, l;
       
       //######### HLT PARS ##########
       double mult_IM = mult_func_IM_O1_HSL(Ergs[eg],Zetas[is]*fabs((Ergs[eg]-Emin)),Emin);
       double mult_RE = mult_func_RE_O1_HSL(Ergs[eg],Zetas[is]*fabs((Ergs[eg]-Emin)),Emin);
       if(Use_conserved_current) { mult_RE /=100; mult_IM /=100;}

       double mult_IM_0 = mult_func_IM(-Ergs[eg], Zetas[is]*fabs((Ergs[eg]-Emin)),Emin);
       double mult_RE_0 = mult_func_RE(-Ergs[eg], Zetas[is]*fabs((Ergs[eg]-Emin)),Emin);
       
       double Ag_target= Ag_func(erg,sigma, aE0);
       int tmax= 34;
       double Emax=0.0;
       int Is_Emax_Finite=(Emax > 1e-10)?1:0;
       double alpha=1.99;
       int prec=128;
       //############################
       

       //HSL
       H1[is].distr_list[eg].RE =  Get_Laplace_transfo(  erg_HSL,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C1_HLT, syst, mult_RE, l, "TANT", "HSL", "1_RE", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin", Cov_1, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);
       H1[is].distr_list[eg].RE = H1[is].distr_list[eg].RE.ave() + ( H1[is].distr_list[eg].RE - H1[is].distr_list[eg].RE.ave())*sqrt( 1.0 + pow( syst/H1[is].distr_list[eg].RE.err(),2));
       H1[is].distr_list[eg].RE = H1[is].distr_list[eg].RE +  Get_Laplace_transfo(  erg_HSL,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE_0, C1_HLT, syst, mult_RE_0, l, "TANT", "HSL", "1_RE", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin_zero", Cov_1, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);

       H1[is].distr_list[eg].IM =  Get_Laplace_transfo(  erg_HSL,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C1_HLT, syst, mult_IM, l, "TANT", "HSL", "1_IM", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin", Cov_1, fake_func,0, fake_func_d ,  Is_Emax_Finite,  Emax, alpha, 1);
       H1[is].distr_list[eg].IM = H1[is].distr_list[eg].IM.ave() + ( H1[is].distr_list[eg].IM - H1[is].distr_list[eg].IM.ave())*sqrt( 1.0 + pow( syst/H1[is].distr_list[eg].IM.err(),2));
       H1[is].distr_list[eg].IM = H1[is].distr_list[eg].IM +  Get_Laplace_transfo(  erg_HSL,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM_0, C1_HLT, syst, mult_IM_0, l, "TANT", "HSL", "1_IM", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin_zero", Cov_1, fake_func,0, fake_func_d , Is_Emax_Finite, Emax, alpha, 1);


       //update multiplicities for O2
       mult_IM = mult_func_IM_O2_HSL(Ergs[eg],Zetas[is]*fabs((Ergs[eg]-Emin)),Emin);
       mult_RE = mult_func_RE_O2_HSL(Ergs[eg],Zetas[is]*fabs((Ergs[eg]-Emin)),Emin);

       H2[is].distr_list[eg].RE =  Get_Laplace_transfo(  erg_HSL,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, C2_HLT, syst, mult_RE, l, "TANT", "HSL", "2_RE", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin", Cov_2, fake_func,0, fake_func_d ,  Is_Emax_Finite,  Emax, alpha, 1);
       H2[is].distr_list[eg].RE = H2[is].distr_list[eg].RE.ave() + ( H2[is].distr_list[eg].RE - H2[is].distr_list[eg].RE.ave())*sqrt( 1.0 + pow( syst/H2[is].distr_list[eg].RE.err(),2));
       H2[is].distr_list[eg].RE = H2[is].distr_list[eg].RE + Get_Laplace_transfo(  erg_HSL,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE_0, C2_HLT, syst, mult_RE_0, l, "TANT", "HSL", "2_RE", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin_zero", Cov_2, fake_func,0, fake_func_d ,  Is_Emax_Finite,  Emax, alpha, 1);
       
       H2[is].distr_list[eg].IM =  Get_Laplace_transfo(  erg_HSL,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, C2_HLT, syst, mult_IM, l, "TANT", "HSL", "2_IM", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin", Cov_2, fake_func,0, fake_func_d ,  Is_Emax_Finite,  Emax, alpha, 1);
       H2[is].distr_list[eg].IM = H2[is].distr_list[eg].IM.ave() + ( H2[is].distr_list[eg].IM - H2[is].distr_list[eg].IM.ave())*sqrt( 1.0 + pow( syst/H2[is].distr_list[eg].IM.err(),2));
       H2[is].distr_list[eg].IM = H2[is].distr_list[eg].IM +  Get_Laplace_transfo(  erg_HSL,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM_0, C2_HLT, syst, mult_IM_0, l, "TANT", "HSL", "2_IM", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin_zero", Cov_2, fake_func,0, fake_func_d ,  Is_Emax_Finite,  Emax, alpha, 1);


       //######### HLT PARS FOR HSS ###################
       Emin=3.70; //GeV
       aE0 = Emin*a_GeV;
       sigma= Zetas[is]*fabs(( fabs(erg - aE0)+smin));
       sigma_GeV= sigma/a_GeV;
       mult_RE = mult_func_RE(Ergs[eg],Zetas[is]*fabs(Ergs[eg]-Emin),Emin);
       mult_IM = mult_func_IM(Ergs[eg],Zetas[is]*fabs(Ergs[eg]-Emin),Emin);
       mult_RE_0 = mult_func_RE(-Ergs[eg],Zetas[is]*fabs(Ergs[eg]-Emin),Emin);
       mult_IM_0 = mult_func_IM(-Ergs[eg],Zetas[is]*fabs(Ergs[eg]-Emin),Emin);
       if(Use_conserved_current) { mult_RE /=100; mult_IM /=100; mult_RE_0 /=100; mult_IM_0 /= 100;}
       Ag_target= Ag_func(erg,sigma, aE0);
       //##############################################
       
       
       //HSS
       HS1[is].distr_list[eg].RE =  Get_Laplace_transfo(  erg,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, CS1_HLT, syst, mult_RE, l, "TANT", "HSS", "1_RE", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin", Cov_S1, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, alpha, 1);
       HS1[is].distr_list[eg].RE = HS1[is].distr_list[eg].RE.ave() + ( HS1[is].distr_list[eg].RE - HS1[is].distr_list[eg].RE.ave())*sqrt( 1.0 + pow( syst/HS1[is].distr_list[eg].RE.err(),2));
       HS1[is].distr_list[eg].RE = HS1[is].distr_list[eg].RE + Get_Laplace_transfo(  erg,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE_0, CS1_HLT, syst, mult_RE_0, l, "TANT", "HSS", "1_RE", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin_zero", Cov_S1, fake_func,0, fake_func_d ,  Is_Emax_Finite, Emax, alpha, 1); 

             
       HS1[is].distr_list[eg].IM =  Get_Laplace_transfo(  erg,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, CS1_HLT, syst, mult_IM, l, "TANT", "HSS", "1_IM", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin", Cov_S1, fake_func,0, fake_func_d ,  Is_Emax_Finite,  Emax, alpha, 1);
       HS1[is].distr_list[eg].IM = HS1[is].distr_list[eg].IM.ave() + ( HS1[is].distr_list[eg].IM - HS1[is].distr_list[eg].IM.ave())*sqrt( 1.0 + pow( syst/HS1[is].distr_list[eg].IM.err(),2));
       HS1[is].distr_list[eg].IM = HS1[is].distr_list[eg].IM +  Get_Laplace_transfo(  erg,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM_0, CS1_HLT, syst, mult_IM_0, l, "TANT", "HSS", "1_IM", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin_zero", Cov_S1, fake_func,0, fake_func_d ,  Is_Emax_Finite,  Emax, alpha, 1);



       //update multiplicites for O2 Bs-> eta_ss' + ll
       mult_RE = mult_func_RE_O2(Ergs[eg],Zetas[is]*fabs(Ergs[eg]-Emin),Emin);
       mult_IM = mult_func_IM_O2(Ergs[eg],Zetas[is]*fabs(Ergs[eg]-Emin),Emin);
       mult_RE_0 = mult_func_RE_O2(-Ergs[eg],Zetas[is]*fabs(Ergs[eg]-Emin),Emin);
       mult_IM_0 = mult_func_IM_O2(-Ergs[eg],Zetas[is]*fabs(Ergs[eg]-Emin),Emin);
       Ag_target= Ag_func(erg,sigma, aE0);
       
       
             
       HS2[is].distr_list[eg].RE =  Get_Laplace_transfo(  erg,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE, CS2_HLT, syst, mult_RE, l, "TANT", "HSS", "2_RE", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin", Cov_S2, fake_func,0, fake_func_d ,  Is_Emax_Finite,  Emax, alpha, 1);
       HS2[is].distr_list[eg].RE = HS2[is].distr_list[eg].RE.ave() + ( HS2[is].distr_list[eg].RE - HS2[is].distr_list[eg].RE.ave())*sqrt( 1.0 + pow( syst/HS2[is].distr_list[eg].RE.err(),2));
       HS2[is].distr_list[eg].RE = HS2[is].distr_list[eg].RE +  Get_Laplace_transfo(  erg,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_RE_0, CS2_HLT, syst, mult_RE_0, l, "TANT", "HSS", "2_RE", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin_zero", Cov_S2, fake_func,0, fake_func_d ,  Is_Emax_Finite,  Emax, alpha, 1);
       
       HS2[is].distr_list[eg].IM =  Get_Laplace_transfo(  erg,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM, CS2_HLT, syst, mult_IM, l, "TANT", "HSS", "2_IM", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin", Cov_S2, fake_func,0, fake_func_d ,  Is_Emax_Finite,  Emax, alpha, 1);
       HS2[is].distr_list[eg].IM = HS2[is].distr_list[eg].IM.ave() + ( HS2[is].distr_list[eg].IM - HS2[is].distr_list[eg].IM.ave())*sqrt( 1.0 + pow( syst/HS2[is].distr_list[eg].IM.err(),2));
       HS2[is].distr_list[eg].IM = HS2[is].distr_list[eg].IM +  Get_Laplace_transfo(  erg,  sigma, aE0,  T_HLT, tmax , prec, "Erg_"+to_string_with_precision(Ergs[eg],3)+"_s_"+to_string_with_precision(sigma_GeV,3),K_IM_0, CS2_HLT, syst, mult_IM_0, l, "TANT", "HSS", "2_IM", Ag_target,0, -1.0*Get_id_distr(Njacks_JJ,UseJack) , 0.0 , "penguin_zero", Cov_S2, fake_func,0, fake_func_d ,  Is_Emax_Finite,  Emax, alpha, 1);


              
     }


     

     

     cout<<"Zeta: "<<Zetas[is]<<" computed!"<<endl<<flush;
   }

     

   //print the results

   //print as a function of E for fixed sigma

     /*
     vector<complex_distr_t_list> FIN_H1, FIN_H2, FIN_HS1, FIN_HS2;
     
     for(auto & s: Zetas) {
       FIN_H1.emplace_back( UseJack, Ergs_NZ.size());
       FIN_H2.emplace_back( UseJack, Ergs_NZ.size());
       
       FIN_HS1.emplace_back( UseJack, Ergs_NZ.size());
       FIN_HS2.emplace_back( UseJack, Ergs_NZ.size());
     }
     
     
     for(int s=0;s<(signed)Zetas.size();s++) {
       for(int e=0;e<Ergs_NZ.size();e++) {
	 
	 FIN_H1[s].distr_list[e] = H1[s].distr_list[e+1] - H1[s].distr_list[0];
	 FIN_H2[s].distr_list[e] = H2[s].distr_list[e+1] - H2[s].distr_list[0];
	 
	 FIN_HS1[s].distr_list[e] = HS1[s].distr_list[e+1] - HS1[s].distr_list[0];
	 FIN_HS2[s].distr_list[e] = HS2[s].distr_list[e+1] - HS2[s].distr_list[0];
	 
       }
     }
     
     
     */

   for(int is=0; is < (signed)Zetas.size(); is++) {

     /*
     Print_To_File({}, {Ergs_NZ, FIN_H1[is].RE_ave(), FIN_H1[is].RE_err(), FIN_H1[is].IM_ave(), FIN_H1[is].IM_err() }, "../data/"+out_dir+"/spectral_reconstruction/HSL/Erg/H1_sub_z_"+to_string_with_precision(Zetas[is],2), "", "");
     Print_To_File({}, {Ergs_NZ, FIN_H2[is].RE_ave(), FIN_H2[is].RE_err(), FIN_H2[is].IM_ave(), FIN_H2[is].IM_err() }, "../data/"+out_dir+"/spectral_reconstruction/HSL/Erg/H2_sub_z_"+to_string_with_precision(Zetas[is],2), "", "");
     
     Print_To_File({}, {Ergs_NZ, FIN_HS1[is].RE_ave(), FIN_HS1[is].RE_err(), FIN_HS1[is].IM_ave(), FIN_HS1[is].IM_err() }, "../data/"+out_dir+"/spectral_reconstruction/HSS/Erg/H1_sub_z_"+to_string_with_precision(Zetas[is],2), "", "");
     Print_To_File({}, {Ergs_NZ, FIN_HS2[is].RE_ave(), FIN_HS2[is].RE_err(), FIN_HS2[is].IM_ave(), FIN_HS2[is].IM_err() }, "../data/"+out_dir+"/spectral_reconstruction/HSS/Erg/H2_sub_z_"+to_string_with_precision(Zetas[is],2), "", "");
     */

     Print_To_File({}, {Ergs_HSL, H1[is].RE_ave(), H1[is].RE_err(), H1[is].IM_ave(), H1[is].IM_err() }, "../data/"+out_dir+"/spectral_reconstruction/HSL/Erg/H1_z_"+to_string_with_precision(Zetas[is],2), "", "");
     Print_To_File({}, {Ergs_HSL, H2[is].RE_ave(), H2[is].RE_err(), H2[is].IM_ave(), H2[is].IM_err() }, "../data/"+out_dir+"/spectral_reconstruction/HSL/Erg/H2_z_"+to_string_with_precision(Zetas[is],2), "", "");
     
     Print_To_File({}, {Ergs, HS1[is].RE_ave(), HS1[is].RE_err(), HS1[is].IM_ave(), HS1[is].IM_err() }, "../data/"+out_dir+"/spectral_reconstruction/HSS/Erg/H1_z_"+to_string_with_precision(Zetas[is],2), "", "");
     Print_To_File({}, {Ergs, HS2[is].RE_ave(), HS2[is].RE_err(), HS2[is].IM_ave(), HS2[is].IM_err() }, "../data/"+out_dir+"/spectral_reconstruction/HSS/Erg/H2_z_"+to_string_with_precision(Zetas[is],2), "", "");
     
     
   }
   


   //print as a function of sigma at fixed energy

   vector<complex_distr_t_list> H1_ERG,  H2_ERG;
   vector<complex_distr_t_list> HS1_ERG, HS2_ERG;


   for(int ie=0; ie< (signed)Ergs.size(); ie++) {
     H1_ERG.emplace_back( UseJack, Zetas.size());
     H2_ERG.emplace_back( UseJack, Zetas.size());
     
     HS1_ERG.emplace_back( UseJack, Zetas.size());
     HS2_ERG.emplace_back( UseJack, Zetas.size());
     
     for(int is=0; is< (signed)Zetas.size(); is++) {

       H1_ERG[ie].distr_list[is] = H1[is].distr_list[ie];
       H2_ERG[ie].distr_list[is] = H2[is].distr_list[ie];
              
       HS1_ERG[ie].distr_list[is] = HS1[is].distr_list[ie];
       HS2_ERG[ie].distr_list[is] = HS2[is].distr_list[ie];
                
     }     
   }

   for(int ie=0; ie < (signed)Ergs.size(); ie++) {

     Print_To_File({}, {Zetas, H1_ERG[ie].RE_ave(), H1_ERG[ie].RE_err(), H1_ERG[ie].IM_ave(), H1_ERG[ie].IM_err() }, "../data/"+out_dir+"/spectral_reconstruction/HSL/eps/H1_erg_"+to_string_with_precision(Ergs_HSL[ie],2), "", "");
     Print_To_File({}, {Zetas, H2_ERG[ie].RE_ave(), H2_ERG[ie].RE_err(), H2_ERG[ie].IM_ave(), H2_ERG[ie].IM_err() }, "../data/"+out_dir+"/spectral_reconstruction/HSL/eps/H2_erg_"+to_string_with_precision(Ergs_HSL[ie],2), "", "");
     
     Print_To_File({}, {Zetas, HS1_ERG[ie].RE_ave(), HS1_ERG[ie].RE_err(), HS1_ERG[ie].IM_ave(), HS1_ERG[ie].IM_err() }, "../data/"+out_dir+"/spectral_reconstruction/HSS/eps/H1_erg_"+to_string_with_precision(Ergs[ie],2), "", "");
     Print_To_File({}, {Zetas, HS2_ERG[ie].RE_ave(), HS2_ERG[ie].RE_err(), HS2_ERG[ie].IM_ave(), HS2_ERG[ie].IM_err() }, "../data/"+out_dir+"/spectral_reconstruction/HSS/eps/H2_erg_"+to_string_with_precision(Ergs[ie],2), "", "");
     

   }



   //print eps->0 extrapolation at 4.75 GeV


   
     
   class ipar {
   public:
     ipar() : val_RE(0.0), err_RE(0.0), val_IM(0.0), err_IM(0.0), eps(0.0) {}
     double val_RE,err_RE,val_IM, err_IM, eps;
   };
   
     class fpar {
     public:
       fpar() {}
       fpar(const Vfloat &par) {
	 if((signed)par.size() != 3) crash("In class fpar  class constructor Vfloat par has size != 2");
	 R=par[0];
	 A=par[1];
	 A2=par[2];
       }
       double R,A,A2;
     };


     int Nmeas= Zetas.size() -1;
     
     //fit on mean values to get ch2
     bootstrap_fit<fpar,ipar> bf(Njacks_JJ);
     bf.set_warmup_lev(1); //sets warmup
     bf.Set_number_of_measurements(Nmeas);
     bf.Set_verbosity(1);
     bf.Add_par("R", 1.0, 0.1);
     bf.Add_par("A", 1.0, 0.1);
     bf.Add_par("A2", 1.0, 0.1);
     //for mean values
     bootstrap_fit<fpar,ipar> bf_ch2(1);
     bf_ch2.set_warmup_lev(1); //sets warmup
     bf_ch2.Set_number_of_measurements(Nmeas);
     bf_ch2.Set_verbosity(1);
     bf_ch2.Add_par("R", 1.0, 0.1);
     bf_ch2.Add_par("A", 1.0, 0.1);
     bf_ch2.Add_par("A2", 1.0, 0.1);
     bf.Fix_par("A2", 0.0);
     bf_ch2.Fix_par("A2",0.0);
     
     //##############################//

          
     int only_four_finest=true;
     int RE_IM=1;
     
     //ansatz
     bf.ansatz=  [&only_four_finest](const fpar &p, const ipar &ip) {
       if(only_four_finest && ip.eps > 1.1) { return 0.0 ; }
       return p.R + p.A*pow(ip.eps,1) + p.A2*pow(ip.eps,2);
       
     };


     
     bf.measurement=  [&RE_IM, &only_four_finest ](const fpar &p, const ipar &ip) {
       if(only_four_finest && ip.eps > 1.1) { return 0.0 ;}
       if(RE_IM==0) return ip.val_IM;
       return ip.val_RE;       
     };
     bf.error=  [&RE_IM, &only_four_finest ](const fpar &p, const ipar &ip) {
      
       if(RE_IM==0) return ip.err_IM;
       return ip.err_RE;       
     };

     bf_ch2.ansatz= bf.ansatz;
     bf_ch2.measurement = bf.measurement;
     bf_ch2.error = bf.error;

         
  
     //fill the data
     vector<vector<ipar>> data(Njacks_JJ);
     vector<vector<ipar>> data_ch2(1);
     //allocate space for output result

     cout<<"Nmeas: "<<Nmeas<<endl;
     
     for(auto &data_iboot: data) data_iboot.resize(Nmeas);
     for(auto &data_iboot: data_ch2) data_iboot.resize(Nmeas);
     for(int ijack=0;ijack<Njacks_JJ;ijack++) {
       for(int imeas=0;imeas<Nmeas;imeas++) {

	 data[ijack][imeas].val_RE= HS1_ERG[id_475].distr_list[imeas+1].RE.distr[ijack];
	 data[ijack][imeas].err_RE = HS1_ERG[id_475].RE_err(imeas+1);
	 data[ijack][imeas].val_IM= HS1_ERG[id_475].distr_list[imeas+1].IM.distr[ijack];
	 data[ijack][imeas].err_IM = HS1_ERG[id_475].IM_err(imeas+1);
	 data[ijack][imeas].eps = Zetas[imeas+1];
	 //mean values
	 if(ijack==0) {
	   data_ch2[ijack][imeas].val_RE = HS1_ERG[id_475].RE_ave(imeas+1);
	   data_ch2[ijack][imeas].val_IM = HS1_ERG[id_475].IM_ave(imeas+1);
	   data_ch2[ijack][imeas].err_RE = HS1_ERG[id_475].RE_err(imeas+1);
	   data_ch2[ijack][imeas].err_IM = HS1_ERG[id_475].IM_err(imeas+1);
	   data_ch2[ijack][imeas].eps=Zetas[imeas+1];
	 }
       }
     }
     
         
     //append
     bf.Append_to_input_par(data);
     bf_ch2.Append_to_input_par(data_ch2);
     //fit
     boot_fit_data<fpar> Bt_fit_RE;
     boot_fit_data<fpar> Bt_fit_RE_ch2;
     only_four_finest=true;
     Bt_fit_RE= bf.Perform_bootstrap_fit();
     Bt_fit_RE_ch2= bf_ch2.Perform_bootstrap_fit();
     bf.Release_par("A2");
     bf_ch2.Release_par("A2");
     only_four_finest=false;
     boot_fit_data<fpar> Bt_fit_RE_eps2;
     boot_fit_data<fpar> Bt_fit_RE_eps2_ch2;
     Bt_fit_RE_eps2= bf.Perform_bootstrap_fit();
     Bt_fit_RE_eps2_ch2= bf_ch2.Perform_bootstrap_fit();
     
     bf.Fix_par("A2",0.0);
     bf_ch2.Fix_par("A2",0.0);
     only_four_finest=true;
     RE_IM=0;
     boot_fit_data<fpar> Bt_fit_IM;
     boot_fit_data<fpar> Bt_fit_IM_ch2;
     Bt_fit_IM= bf.Perform_bootstrap_fit();
     Bt_fit_IM_ch2= bf_ch2.Perform_bootstrap_fit();
     bf.Release_par("A2");
     bf_ch2.Release_par("A2");
     only_four_finest=false;
     boot_fit_data<fpar> Bt_fit_IM_eps2;
     boot_fit_data<fpar> Bt_fit_IM_eps2_ch2;
     Bt_fit_IM_eps2= bf.Perform_bootstrap_fit();
     Bt_fit_IM_eps2_ch2= bf.Perform_bootstrap_fit();

     //print ch2/dof

     cout<<"ch2(RE,linear): "<<Bt_fit_RE_ch2.get_ch2_ave()/(Nmeas - 2)<<endl;
     cout<<"ch2(RE,quad): "<<Bt_fit_RE_eps2_ch2.get_ch2_ave()/(Nmeas-3)<<endl;
     cout<<"ch2(IM,linear): "<<Bt_fit_IM_ch2.get_ch2_ave()/(Nmeas - 2)<<endl;
     cout<<"ch2(IM,quad): "<<Bt_fit_IM_eps2_ch2.get_ch2_ave()/(Nmeas-3)<<endl;


     //retrieve parameters

     distr_t H_RE(UseJack), H_IM(UseJack), H_RE_Q(UseJack), H_IM_Q(UseJack);
     distr_t A1_RE(UseJack), A1_IM(UseJack), A1_RE_Q(UseJack), A1_IM_Q(UseJack);
     distr_t A2_RE(UseJack), A2_IM(UseJack), A2_RE_Q(UseJack), A2_IM_Q(UseJack);


     for(int ijack=0;ijack<Njacks_JJ;ijack++) {

       H_RE.distr.push_back( Bt_fit_RE.par[ijack].R);
       H_IM.distr.push_back( Bt_fit_IM.par[ijack].R);
       H_RE_Q.distr.push_back( Bt_fit_RE_eps2.par[ijack].R);
       H_IM_Q.distr.push_back( Bt_fit_IM_eps2.par[ijack].R);

       A1_RE.distr.push_back( Bt_fit_RE.par[ijack].A);
       A1_IM.distr.push_back( Bt_fit_IM.par[ijack].A);
       A1_RE_Q.distr.push_back( Bt_fit_RE_eps2.par[ijack].A);
       A1_IM_Q.distr.push_back( Bt_fit_IM_eps2.par[ijack].A);

       
       A2_RE.distr.push_back( Bt_fit_RE.par[ijack].A2);
       A2_IM.distr.push_back( Bt_fit_IM.par[ijack].A2);
       A2_RE_Q.distr.push_back( Bt_fit_RE_eps2.par[ijack].A2);
       A2_IM_Q.distr.push_back( Bt_fit_IM_eps2.par[ijack].A2);
  
     }


     //fit func

     Vfloat alpha_to_print;
     distr_t_list func_RE(UseJack), func_IM(UseJack), func_RE_Q(UseJack), func_IM_Q(UseJack);
     int Nsteps=200;
     for(int istep=0;istep<Nsteps;istep++) { alpha_to_print.push_back( istep*3.0/(Nsteps-1.0) );}
     for(int a=0;a<(signed)alpha_to_print.size(); a++) {
       double alpha=alpha_to_print[a];
       func_RE.distr_list.push_back( H_RE + A1_RE*alpha + A2_RE*pow(alpha,2)) ;
       func_IM.distr_list.push_back( H_IM + A1_IM*alpha + A2_IM*pow(alpha,2)) ;
       func_RE_Q.distr_list.push_back( H_RE_Q + A1_RE_Q*alpha + A2_RE_Q*pow(alpha,2));
       func_IM_Q.distr_list.push_back( H_IM_Q + A1_IM_Q*alpha + A2_IM_Q*pow(alpha,2));
     }


     //print
     boost::filesystem::create_directory("../data/PINGU/spectral_reconstruction/fit");
     
     Print_To_File({}, { alpha_to_print, func_RE.ave(), func_RE.err(), func_RE_Q.ave(), func_RE_Q.err() }, "../data/PINGU/spectral_reconstruction/fit/H1_RE_5.28","","");
     Print_To_File({}, { alpha_to_print, func_IM.ave(), func_IM.err(), func_IM_Q.ave(), func_IM_Q.err() }, "../data/PINGU/spectral_reconstruction/fit/H1_IM_5.28","","");
     

     
   
   

   }



   //print model_to_file

     
   

   for(int is=0; is < (signed)Zetas_mod.size(); is++) {

     Print_To_File({}, {Emod, MOD_HSL_RE[is], MOD_HSL_IM[is] }, "../data/"+out_dir+"/spectral_reconstruction/HSL/Erg/H1_mod_z_"+to_string_with_precision(Zetas_mod[is],2), "", "");
     Print_To_File({}, {Emod, MOD_HSS_RE[is], MOD_HSS_IM[is] }, "../data/"+out_dir+"/spectral_reconstruction/HSS/Erg/H1_mod_z_"+to_string_with_precision(Zetas_mod[is],2), "", "");
     
   }
   

   
   cout<<"done!"<<endl;




   return ;
}



void TEST_PI_Q2() {


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





    data_t C_1, C_2, C_3;


    C_1.Read("../PI_Q2", "mes_contr_PT2_CONS_1", "S0V1", Sort_light_confs);
    C_2.Read("../PI_Q2", "mes_contr_PT2_CONS_2", "S0V2", Sort_light_confs);
    C_3.Read("../PI_Q2", "mes_contr_PT2_CONS_3", "S0V3", Sort_light_confs);

    int Njacks=50;

    //set up statistical analysis
    CorrAnalysis Corr(UseJack,Njacks,100);
    Corr.Nt=128;
    int T=128;
    Corr.Reflection_sign=1;
    Corr.Perform_Nt_t_average=1;

    boost::filesystem::create_directory("../data/PI_Q2_TEST");
    distr_t_list V= (1.0/3.0)*Corr.corr_t( summ_master( C_1.col(0)[0], C_2.col(0)[0], C_3.col(0)[0] ), "../data/PI_Q2_TEST/V");


    double a_B= 0.0795/0.197327;
    int Nqs=100;
    Vfloat Qs;
    for(int i=0;i<Nqs;i++) { Qs.push_back( 0.100 + i*0.02); }

    distr_t_list PI_NAIVE(UseJack,Nqs,Njacks), PI_SUB(UseJack,Nqs,Njacks);

    for(int iq=0;iq<Nqs;iq++) {
      distr_t p_sum_naive(UseJack,Njacks);
      distr_t p_sum_sub(UseJack, Njacks);
      for(int t=0;t<44;t++) {
	p_sum_naive = p_sum_naive + -(( cos(Qs[iq]*a_B*t))/(a_B*a_B*Qs[iq]*Qs[iq]))*(V.distr_list[t]+ (t==0?(0.0*Get_id_jack_distr(Njacks)):V.distr_list[T-t]));
	p_sum_sub = p_sum_sub + (1/(a_B*a_B*Qs[iq]*Qs[iq]))*( 1- cos(Qs[iq]*a_B*t) )*(V.distr_list[t]+ (t==0?(0.0*Get_id_jack_distr(Njacks)):V.distr_list[T-t]));
      }

      PI_NAIVE.distr_list[iq] = p_sum_naive;
      PI_SUB.distr_list[iq] = p_sum_sub;

    }


    //print

    Print_To_File({}, { Qs, PI_NAIVE.ave(), PI_NAIVE.err(), PI_SUB.ave(), PI_SUB.err() } , "../data/PI_Q2_TEST/PI_Q2", "", "");
    
    
    

    











  exit(-1);
  return;
}

  


  
  
   
