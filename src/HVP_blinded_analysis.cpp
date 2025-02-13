#include "../include/HVP_blinded_analysis.h"
#include "Corr_analysis.h"
#include "numerics.h"


const double DTT = 0.5;
const double alpha = 1.0 / 137.035999;
const bool UseJack = true;
const int Njacks = 100;
const bool Skip_GS_analysis = true;
const double l1ph= -0.4; //-0.4
const double l2ph= 4.3; //4.3
const double l3ph= 3.2; //3.2
const double l4ph= 4.4; //4.4
const double s0= 2.0-M_PI/2.0;
const double s1 = M_PI/4.0 - 0.5;
const double s2 = 0.5 - M_PI/8.0;
const double s3 = 3.0*M_PI/16.0 - 0.5;

using namespace std;


void mcorr_HVP() {

  double fm_to_inv_Gev= 1.0/0.197327;
  double qu= 2.0/3.0;
  double qd= -1.0/3.0;
  
   //generate lattice spacings and RCs
   //############################################################################################
   //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
   GaussianMersenne GM(36551294);
   LatticeInfo a_info;
   distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack), a_Z(UseJack), a_E(UseJack);
   distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
   distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
   double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err, a_Z_ave, a_Z_err, a_E_ave, a_E_err;
   double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err;
   double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err;
   a_info.LatInfo_new_ens("cA211a.53.24");
   a_A_ave= a_info.a_from_afp_FLAG;
   a_A_err= a_info.a_from_afp_FLAG_err;
   ZA_A_ave = a_info.Za_WI_strange;
   ZA_A_err = a_info.Za_WI_strange_err;
   ZV_A_ave = a_info.Zv_WI_strange;
   ZV_A_err = a_info.Zv_WI_strange_err;
   a_info.LatInfo_new_ens("cB211b.072.64");
   a_B_ave= a_info.a_from_afp_FLAG;
   a_B_err= a_info.a_from_afp_FLAG_err;
   ZA_B_ave = a_info.Za_WI_FLAG;
   ZA_B_err = a_info.Za_WI_FLAG_err;
   ZV_B_ave = a_info.Zv_WI_FLAG;
   ZV_B_err = a_info.Zv_WI_FLAG_err;
   a_info.LatInfo_new_ens("cC211a.06.80");
   a_C_ave= a_info.a_from_afp_FLAG;
   a_C_err= a_info.a_from_afp_FLAG_err;
   ZA_C_ave = a_info.Za_WI_FLAG;
   ZA_C_err = a_info.Za_WI_FLAG_err;
   ZV_C_ave = a_info.Zv_WI_FLAG;
   ZV_C_err = a_info.Zv_WI_FLAG_err;
   a_info.LatInfo_new_ens("cD211a.054.96");
   a_D_ave= a_info.a_from_afp_FLAG;
   a_D_err= a_info.a_from_afp_FLAG_err;
   ZA_D_ave = a_info.Za_WI_FLAG;
   ZA_D_err = a_info.Za_WI_FLAG_err;
   ZV_D_ave = a_info.Zv_WI_FLAG;
   ZV_D_err = a_info.Zv_WI_FLAG_err;
   a_info.LatInfo_new_ens("cZ211a.077.64");
   a_Z_ave= a_info.a_from_afp_FLAG;
   a_Z_err= a_info.a_from_afp_FLAG_err;
   a_info.LatInfo_new_ens("cE211a.044.112");
   a_E_ave= a_info.a_from_afp_FLAG;
   a_E_err= a_info.a_from_afp_FLAG_err;
   
   for(int ijack=0;ijack<Njacks;ijack++) {


     a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
     a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
     a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
     a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
     a_Z.distr.push_back( fm_to_inv_Gev*( a_Z_ave + GM()*a_Z_err*(1.0/sqrt(Njacks-1.0))));
     a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err*(1.0/sqrt(Njacks-1.0))));
     ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err*(1.0/sqrt(Njacks -1.0)));
     ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err*(1.0/sqrt(Njacks -1.0)));
     ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err*(1.0/sqrt(Njacks -1.0)));
     ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err*(1.0/sqrt(Njacks -1.0)));
     ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err*(1.0/sqrt(Njacks -1.0)));
     ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err*(1.0/sqrt(Njacks -1.0)));
     ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err*(1.0/sqrt(Njacks -1.0)));
     ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err*(1.0/sqrt(Njacks -1.0)));
  
   }


   /* vector<string> Corr_tags({
       "VkVk_OS_mu6.6750e-04_eexIR.bin",
       "VkVk_OS_mu6.6750e-04_stoch.bin",
       "VkVk_OS_mu7.2000e-04_eexIR.bin",
       "VkVk_OS_mu7.2000e-04_stoch.bin",
       "VkVk_tm_mu6.6750e-04_eexIR.bin",
       "VkVk_tm_mu6.6750e-04_stoch.bin",
       "VkVk_tm_mu7.2000e-04_eexIR.bin",
       "VkVk_tm_mu7.2000e-04_stoch.bin",
     });
   */

   
     vector<string> Corr_tags({
       "VkVk_OS_mu6.6750e-04.bin",
       "VkVk_OS_mu7.2000e-04.bin",
       "VkVk_tm_mu6.6750e-04.bin",
       "VkVk_tm_mu7.2000e-04.bin",
       });


   VVVfloat Corrs(Corr_tags.size());

   cout<<"Number of correlators: "<<Corr_tags.size()<<endl;

   int Nconfs=770;
   int Nhits=128;
   int Nsubs=1;
   int TT=128;


   for(int id=0; id<(signed)Corr_tags.size(); id++) {


     
     string ch= Corr_tags[id];
     
  	
     FILE *stream = fopen(("../mass_corr_LMA/B64/"+ch).c_str(), "rb");

     /*
     double Nconfs, TT, Nhits, Nsubs;
     bin_read(Nconfs, stream);
     D(3);
     cout<<"Nconfs: "<<Nconfs<<endl;

     bin_read(Nhits, stream);
     cout<<"Nhits: "<<Nhits<<endl;
     bin_read(TT, stream);
     cout<<"TT: "<<TT<<endl;
     bin_read(Nsubs, stream);
     cout<<"Nsubs: "<<Nsubs<<endl;
     */

     Corrs[id].resize(TT);
     for(auto & cc: Corrs[id]) cc.resize(Nconfs,0.0);
      
     cout<<"Nconfs: "<<Nconfs<<endl;
     cout<<"TT: "<<TT<<" "<<TT/2+1<<endl;
     cout<<"Nhits: "<<Nhits<<endl;
     cout<<"Nsubs: "<<Nsubs<<endl;
     
     
     
     for(size_t iconf=0;iconf<Nconfs;iconf++) {
       for(size_t t=0;t<TT/2+1;t++) {
	 for(size_t is=0;is<Nsubs;is++) {
	   
	  double c;
	  bin_read(c, stream);
	  Corrs[id][t][iconf] += c/Nsubs;
	 }
	 //symmetrize
	 if( t != 0) {
	  Corrs[id][TT-t][iconf] = Corrs[id][t][iconf];
	 }
       }
    }
     
     
     fclose(stream);
     
   }
   
   
   //analyze correlators

   
   VVfloat C_dm_OS_raw= summ_master( Corrs[0], Multiply_Vvector_by_scalar( Corrs[1], -1.0));
   VVfloat C_dm_TM_raw= summ_master( Corrs[2], Multiply_Vvector_by_scalar( Corrs[3], -1.0));

   //make jackk analysis

   CorrAnalysis Corr(UseJack, Njacks,100);
   Corr.Perform_Nt_t_average=0;
   Corr.Nt = 128;

   distr_t a_distr= a_B;
   distr_t Zv= ZV_B;
   distr_t Za= ZA_A;

   boost::filesystem::create_directory("../data/dm_corr");
   boost::filesystem::create_directory("../data/dm_corr/B64");
   distr_t_list C_dm_OS = Zv*Zv*1e10*( qu*qu + qd*qd)*Corr.corr_t( C_dm_OS_raw, "../data/dm_corr/B64/C_dm_OS");
   distr_t_list C_dm_TM = Za*Za*1e10*( qu*qu + qd*qd)*Corr.corr_t( C_dm_TM_raw, "../data/dm_corr/B64/C_dm_TM");


   //calculate HVP up to 2.8 fm

   auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
   distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);


   distr_t p_sum_TM= 0.0*Get_id_jack_distr(Njacks);
   distr_t p_sum_OS= 0.0*Get_id_jack_distr(Njacks);

   distr_t amu_TM(UseJack), amu_OS(UseJack);

   distr_t_list amu_tcut_TM(UseJack), amu_tcut_OS(UseJack);

   for(int t=0;t<Corr.Nt/2;t++) {

     p_sum_TM = p_sum_TM +  4.0*w(t,1)*pow(alpha,2)*C_dm_TM.distr_list[t]*Ker.distr_list[t] ;
     p_sum_OS = p_sum_OS +  4.0*w(t,1)*pow(alpha,2)*C_dm_OS.distr_list[t]*Ker.distr_list[t] ;

     amu_tcut_TM.distr_list.push_back( p_sum_TM);
     amu_tcut_OS.distr_list.push_back( p_sum_OS);


   }


   Print_To_File({}, { amu_tcut_TM.ave(), amu_tcut_TM.err()}, "../data/dm_corr/B64/damu_TM", "", "");
   Print_To_File({}, { amu_tcut_OS.ave(), amu_tcut_OS.err()}, "../data/dm_corr/B64/damu_OS", "", "");











   return;
   
   




}


void HVP_blinded_analysis() {
  
  double fm_to_inv_Gev= 1.0/0.197327;
  double qu= 2.0/3.0;
  double qd= -1.0/3.0;
  
   //generate lattice spacings and RCs
   //############################################################################################
   //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
   GaussianMersenne GM(36551294);
   LatticeInfo a_info;
   distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack), a_Z(UseJack), a_E(UseJack);
   distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
   distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
   double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err, a_Z_ave, a_Z_err, a_E_ave, a_E_err;
   double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err;
   double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err;
   a_info.LatInfo_new_ens("cA211a.53.24");
   a_A_ave= a_info.a_from_afp_FLAG;
   a_A_err= a_info.a_from_afp_FLAG_err;
   ZA_A_ave = a_info.Za_WI_strange;
   ZA_A_err = a_info.Za_WI_strange_err;
   ZV_A_ave = a_info.Zv_WI_strange;
   ZV_A_err = a_info.Zv_WI_strange_err;
   a_info.LatInfo_new_ens("cB211b.072.64");
   a_B_ave= a_info.a_from_afp_FLAG;
   a_B_err= a_info.a_from_afp_FLAG_err;
   ZA_B_ave = a_info.Za_WI_FLAG;
   ZA_B_err = a_info.Za_WI_FLAG_err;
   ZV_B_ave = a_info.Zv_WI_FLAG;
   ZV_B_err = a_info.Zv_WI_FLAG_err;
   a_info.LatInfo_new_ens("cC211a.06.80");
   a_C_ave= a_info.a_from_afp_FLAG;
   a_C_err= a_info.a_from_afp_FLAG_err;
   ZA_C_ave = a_info.Za_WI_FLAG;
   ZA_C_err = a_info.Za_WI_FLAG_err;
   ZV_C_ave = a_info.Zv_WI_FLAG;
   ZV_C_err = a_info.Zv_WI_FLAG_err;
   a_info.LatInfo_new_ens("cD211a.054.96");
   a_D_ave= a_info.a_from_afp_FLAG;
   a_D_err= a_info.a_from_afp_FLAG_err;
   ZA_D_ave = a_info.Za_WI_FLAG;
   ZA_D_err = a_info.Za_WI_FLAG_err;
   ZV_D_ave = a_info.Zv_WI_FLAG;
   ZV_D_err = a_info.Zv_WI_FLAG_err;
   a_info.LatInfo_new_ens("cZ211a.077.64");
   a_Z_ave= a_info.a_from_afp_FLAG;
   a_Z_err= a_info.a_from_afp_FLAG_err;
   a_info.LatInfo_new_ens("cE211a.044.112");
   a_E_ave= a_info.a_from_afp_FLAG;
   a_E_err= a_info.a_from_afp_FLAG_err;
   
   for(int ijack=0;ijack<Njacks;ijack++) {


     a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
     a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
     a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
     a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
     a_Z.distr.push_back( fm_to_inv_Gev*( a_Z_ave + GM()*a_Z_err*(1.0/sqrt(Njacks-1.0))));
     a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err*(1.0/sqrt(Njacks-1.0))));
     ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err*(1.0/sqrt(Njacks -1.0)));
     ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err*(1.0/sqrt(Njacks -1.0)));
     ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err*(1.0/sqrt(Njacks -1.0)));
     ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err*(1.0/sqrt(Njacks -1.0)));
     ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err*(1.0/sqrt(Njacks -1.0)));
     ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err*(1.0/sqrt(Njacks -1.0)));
     ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err*(1.0/sqrt(Njacks -1.0)));
     ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err*(1.0/sqrt(Njacks -1.0)));
  
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

    auto Sort_light_confs_h5 = [](string A, string B) {


			     //return A<B;
			     
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
							       
		      
			    string rA = A_bis.substr(A_bis.length()-5);
			    string rB = B_bis.substr(B_bis.length()-5);
			    if(rA.substr(0,1) == "r") { 
			      int n1 = stoi(rA.substr(1,1));
			      int n2 = stoi(rB.substr(1,1));
			      if(rA == rB) {
			      if(rA=="r0.h5" || rA=="r2.h5") return conf_num_A > conf_num_B;
			      else if(rA=="r1.h5" || rA=="r3.h5") return conf_num_A < conf_num_B;
			      else crash("stream not recognized");
			      }
			      else return n1<n2;
			    }
			    return A_bis<B_bis;
			  };


  
  VVfloat C_tm, C_OS;
  VVfloat C_tm_C80, C_OS_C80;
  VVfloat C_tm_C80_c1, C_OS_C80_c1;
  VVfloat C_tm_C80_c2, C_OS_C80_c2;
  VVfloat C_tm_C80_c3, C_OS_C80_c3;
  VVfloat C_tm_C80_c4, C_OS_C80_c4;
  VVfloat C_tm_D96, C_OS_D96;

  //###################################################
  //###########    READING  BINARY FILE ################
  
  
  //read binary
  vector<string> Corr_tags({"VkVk_OS.bin", "VkVk_tm.bin"});
  /*vector<string> confs_B64_id;
  ifstream READ_GAUGE_CONFS("../../PEAKY_BLINDER/blinded_confs_RM3/B64/confs_B64.txt");
  while(!READ_GAUGE_CONFS.eof()) {
    string tag;
    READ_GAUGE_CONFS>>tag;
    if(!READ_GAUGE_CONFS.eof()) confs_B64_id.push_back(tag);
  }
  READ_GAUGE_CONFS.close();

  //sort
  vector<string> ordered_confs_B64_id=confs_B64_id;
  sort(ordered_confs_B64_id.begin(), ordered_confs_B64_id.end(), Sort_light_confs);
  cout<<"Ordering confs on B64"<<endl;
  for(auto &c: ordered_confs_B64_id) cout<<c<<endl;
  cout<<"done!"<<endl;
  vector<int> B64_T(confs_B64_id.size());
  for(int i=0;i<(signed)confs_B64_id.size();i++) {
    int found_ens=-1;
    for(int j=0;j<(signed)confs_B64_id.size(); j++) {
      if(confs_B64_id[i] == ordered_confs_B64_id[j]) { B64_T[i] = j; found_ens ++;}
    }
    if(found_ens != 0) crash("cannot match configs after reordering");
  }
  */
  for(int id=0; id<(signed)Corr_tags.size(); id++) {

    
    string ch= Corr_tags[id];
	
    FILE *stream = fopen(("../../PEAKY_BLINDER/blinded_confs_RM3/B64/"+ch).c_str(), "rb");
    size_t Nconfs, TT, Nhits, Nsubs;
    bin_read(Nconfs, stream);
    bin_read(Nhits, stream);
    bin_read(TT, stream);
    bin_read(Nsubs, stream);

    if(id==0) {
      C_OS.resize(TT);
      for(auto & cc: C_OS) cc.resize(Nconfs,0);
    }
    else {
      C_tm.resize(TT);
      for(auto & cc: C_tm) cc.resize(Nconfs,0);
    }
   
     
    
   
     cout<<"Nconfs: "<<Nconfs<<endl;
     cout<<"TT: "<<TT<<" "<<TT/2+1<<endl;
     cout<<"Nhits: "<<Nhits<<endl;
     cout<<"Nsubs: "<<Nsubs<<endl;
    
  
    
    
    for(size_t iconf=0;iconf<Nconfs;iconf++) {
      for(size_t t=0;t<TT/2+1;t++) {
	for(size_t is=0;is<Nsubs;is++) {
	  
	  double c;
	  bin_read(c, stream);
	  if(id==0) {
	    C_OS[t][iconf] += c/Nsubs;
	  }
	  else {
	    C_tm[t][iconf] += c/Nsubs;
	  }
	}
	//symmetrize
	if( t != 0) {
	  if(id==0) {
	    C_OS[TT-t][iconf] = C_OS[t][iconf];
	  }
	  else {
	    C_tm[TT-t][iconf] = C_tm[t][iconf];
	  }
	}
      }
    }
    
   
    
    fclose(stream);
    
  }


  for(int id=0; id<(signed)Corr_tags.size(); id++) {

    
    string ch= Corr_tags[id];
	
    FILE *stream = fopen(("../../PEAKY_BLINDER/blinded_confs_RM3/C80/"+ch).c_str(), "rb");
    size_t Nconfs, TT, Nhits, Nsubs;
    bin_read(Nconfs, stream);
    bin_read(Nhits, stream);
    bin_read(TT, stream);
    bin_read(Nsubs, stream);

    

    if(id==0) {
      C_OS_C80.resize(TT);
      C_OS_C80_c1.resize(TT);
      C_OS_C80_c2.resize(TT);
      C_OS_C80_c3.resize(TT);
      C_OS_C80_c4.resize(TT);
      for(auto & cc: C_OS_C80) cc.resize(Nconfs-218,0);
      for(auto & cc: C_OS_C80_c1) cc.resize(Nconfs-218,0);
      for(auto & cc: C_OS_C80_c2) cc.resize(Nconfs-218,0);
      for(auto & cc: C_OS_C80_c3) cc.resize(Nconfs-218,0);
      for(auto & cc: C_OS_C80_c4) cc.resize(Nconfs-218,0);
    }
    else {
      C_tm_C80.resize(TT);
      C_tm_C80_c1.resize(TT);
      C_tm_C80_c2.resize(TT);
      C_tm_C80_c3.resize(TT);
      C_tm_C80_c4.resize(TT);
      for(auto & cc: C_tm_C80) cc.resize(Nconfs-218,0);
      for(auto & cc: C_tm_C80_c1) cc.resize(Nconfs-218,0);
      for(auto & cc: C_tm_C80_c2) cc.resize(Nconfs-218,0);
      for(auto & cc: C_tm_C80_c3) cc.resize(Nconfs-218,0);
      for(auto & cc: C_tm_C80_c4) cc.resize(Nconfs-218,0);
    }
   
     
    
   
     cout<<"Nconfs: "<<Nconfs<<endl;
     cout<<"TT: "<<TT<<" "<<TT/2+1<<endl;
     cout<<"Nhits: "<<Nhits<<endl;
     cout<<"Nsubs: "<<Nsubs<<endl;
    
  
    
    
    for(size_t iconf=0;iconf<Nconfs;iconf++) {
      for(size_t t=0;t<TT/2+1;t++) {
	for(size_t is=0;is<Nsubs;is++) {
	  
	  double c;
	  bin_read(c, stream);
	  if(iconf >=218 ) {
	  if(id==0) {
	    C_OS_C80[t][iconf-218] += c/Nsubs;
	    if(is < Nsubs/2) C_OS_C80_c1[t][iconf-218] += 2*c/Nsubs;
	    else C_OS_C80_c4[t][iconf-218] += 2*c/Nsubs;
	    
	    if(is < (Nsubs/2) + 1) C_OS_C80_c2[t][iconf-218] += c/( (Nsubs/2.0) +1);
	    if(is < (Nsubs/2) + 2) C_OS_C80_c3[t][iconf-218] += c/( (Nsubs/2.0) +2);
	  }
	  else {
	    C_tm_C80[t][iconf-218] += c/Nsubs;
	    if(is < Nsubs/2) C_tm_C80_c1[t][iconf-218] += 2*c/Nsubs;
	    else C_tm_C80_c4[t][iconf-218] += 2*c/Nsubs;

	    if(is < (Nsubs/2) + 1) C_tm_C80_c2[t][iconf-218] += c/( (Nsubs/2.0) +1);
	    if(is < (Nsubs/2) + 2) C_tm_C80_c3[t][iconf-218] += c/( (Nsubs/2.0) +2);
	  }
	  }
	}
	//symmetrize
	if( t != 0 && iconf >=218) {
	  if(id==0) {
	    C_OS_C80[TT-t][iconf-218] = C_OS_C80[t][iconf-218];
	    C_OS_C80_c1[TT-t][iconf-218] = C_OS_C80_c1[t][iconf-218];
	    C_OS_C80_c2[TT-t][iconf-218] = C_OS_C80_c2[t][iconf-218];
	    C_OS_C80_c3[TT-t][iconf-218] = C_OS_C80_c3[t][iconf-218];
	    C_OS_C80_c4[TT-t][iconf-218] = C_OS_C80_c4[t][iconf-218];
	  }
	  else {
	    C_tm_C80[TT-t][iconf-218] = C_tm_C80[t][iconf-218];
	    C_tm_C80_c1[TT-t][iconf-218] = C_tm_C80_c1[t][iconf-218];
	    C_tm_C80_c2[TT-t][iconf-218] = C_tm_C80_c2[t][iconf-218];
	    C_tm_C80_c3[TT-t][iconf-218] = C_tm_C80_c3[t][iconf-218];
	    C_tm_C80_c4[TT-t][iconf-218] = C_tm_C80_c4[t][iconf-218];
	  }
	}
      }
    }
    
   
    
    fclose(stream);
    
   }

  


  for(int id=0; id<(signed)Corr_tags.size(); id++) {

    
    string ch= Corr_tags[id];
    
    FILE *stream = fopen(("../../PEAKY_BLINDER/blinded_confs_RM3/D96/"+ch).c_str(), "rb");
    size_t Nconfs, TT, Nhits, Nsubs;
    bin_read(Nconfs, stream);
    bin_read(Nhits, stream);
    bin_read(TT, stream);
    bin_read(Nsubs, stream);
    
    if(id==0) {
      C_OS_D96.resize(TT);
      for(auto & cc: C_OS_D96) cc.resize(Nconfs,0);
    }
    else {
      C_tm_D96.resize(TT);
      for(auto & cc: C_tm_D96) cc.resize(Nconfs,0);
    }
   
     
    
   
    cout<<"Nconfs: "<<Nconfs<<endl;
    cout<<"TT: "<<TT<<" "<<TT/2+1<<endl;
    cout<<"Nhits: "<<Nhits<<endl;
    cout<<"Nsubs: "<<Nsubs<<endl;
    
  
    
    
    for(size_t iconf=0;iconf<Nconfs;iconf++) {
      for(size_t t=0;t<TT/2+1;t++) {
	for(size_t is=0;is<Nsubs;is++) {
	  
	  double c;
	  bin_read(c, stream);
	  if(id==0) {
	    C_OS_D96[t][iconf] += c/Nsubs;
	  }
	  else {
	    C_tm_D96[t][iconf] += c/Nsubs;
	  }
	}
	//symmetrize
	if( t != 0) {
	  if(id==0) {
	    C_OS_D96[TT-t][iconf] = C_OS_D96[t][iconf];
	  }
	  else {
	    C_tm_D96[TT-t][iconf] = C_tm_D96[t][iconf];
	  }
	}
      }
    }
    
   
    
    fclose(stream);
    
  }

 

  data_t data_disc_JJ_light, data_disc_JJ, data_disc_JJ_light_stoch, data_disc_JJ_ls; 
  data_t data_disc_OS_JJ_light;
  data_disc_JJ_light.Read("../magnetic_susc_disco/disco", "JJ_light.txt", "VKVK", Sort_light_confs);
  data_disc_JJ.Read("../magnetic_susc_disco/disco", "JJ.txt", "VKVK", Sort_light_confs);
  data_disc_JJ_ls.Read("../magnetic_susc_disco/disco", "JJ_ls.txt", "VKVK", Sort_light_confs);
  data_disc_JJ_light_stoch.Read("../magnetic_susc_disco/light_light_stoch", "disco", "", Sort_light_confs);
  data_disc_OS_JJ_light.Read("../magnetic_susc_disco/disco_P5", "VKVK_diff.txt","VKVK", Sort_light_confs);

   data_t data_disc_JJ_light_impr;
   data_t data_disc_JJ_light_impr_stoch;
   data_disc_JJ_light_impr.Read("../new_isoQCD_loops/disco", "JJ_light.txt", "VKVK", Sort_light_confs);
   data_disc_JJ_light_impr_stoch.Read("../new_isoQCD_loops/disco", "JJ_light_stoch.txt", "VKVK", Sort_light_confs);


   data_t data_disc_JJ_light_LMA21;
   data_t data_disc_JJ_light_stoch_LMA21;
   data_t data_disc_JJ_light_impr_LMA21;
   data_t data_disc_P5_light_LMA21;
   data_t data_disc_P5_light_impr_LMA21;
   data_t data_disc_A0_light_LMA21;
   data_t data_disc_A0_light_impr_LMA21;
   data_t data_disc_bub_P5_light_LMA21;
   data_t data_disc_P5_light;
   data_t data_disc_eta_light;
   data_t data_disc_bub_P5_light;
   data_t data_disc_bub_eta_light;
   data_disc_JJ_light_LMA21.Read("../new_isoQCD_loops/disco_LMA21", "JJ_light.txt", "VKVK", Sort_light_confs);
   data_disc_P5_light_LMA21.Read("../new_isoQCD_loops/disco_LMA21", "P5P5_light.txt", "P5P5", Sort_light_confs);
   data_disc_P5_light_impr_LMA21.Read("../new_isoQCD_loops/disco_LMA21", "P5P5_light_impr.txt", "P5P5", Sort_light_confs);
   data_disc_A0_light_LMA21.Read("../new_isoQCD_loops/disco_LMA21", "A0A0_light.txt", "A0A0", Sort_light_confs);
   data_disc_A0_light_impr_LMA21.Read("../new_isoQCD_loops/disco_LMA21", "A0A0_light_impr.txt", "A0A0", Sort_light_confs);
   data_disc_JJ_light_stoch_LMA21.Read("../new_isoQCD_loops/disco_LMA21", "JJ_light_stoch.txt", "VKVK", Sort_light_confs);
   data_disc_JJ_light_impr_LMA21.Read("../new_isoQCD_loops/disco_LMA21", "JJ_light_impr.txt", "VKVK", Sort_light_confs); 
   data_disc_bub_P5_light_LMA21.Read("../new_isoQCD_loops/disco_LMA21", "bub_P5_light.txt", "P5", Sort_light_confs);
   data_disc_P5_light.Read("../magnetic_susc_disco/disco_P5", "P5P5_diff.txt", "P5P5", Sort_light_confs);
   data_disc_eta_light.Read("../magnetic_susc_disco/disco_P5", "P5P5.txt", "P5P5", Sort_light_confs);
   data_disc_bub_P5_light.Read("../magnetic_susc_disco/disco_P5", "bub_P5_diff.txt", "P5P5", Sort_light_confs);
   data_disc_bub_eta_light.Read("../magnetic_susc_disco/disco_P5", "bub_P5.txt", "P5P5", Sort_light_confs);
  
   data_t data_disc_A0A0_light, data_disc_A0P5_light;

   data_disc_A0A0_light.Read("../magnetic_susc_disco/disco_P5", "A0A0.txt", "", Sort_light_confs);
   data_disc_A0P5_light.Read("../magnetic_susc_disco/disco_P5", "A0P5.txt", "", Sort_light_confs);
    


   data_t data_P5P5_OS, data_A0A0_OS, data_A0P5_OS;
   data_P5P5_OS.Read("../R_ratio_data/light", "mes_contr_2pts_ll_2", "P5P5", Sort_light_confs_h5);
   data_A0A0_OS.Read("../R_ratio_data/light", "mes_contr_2pts_ll_2", "A0A0", Sort_light_confs_h5);
   data_A0P5_OS.Read("../R_ratio_data/light", "mes_contr_2pts_ll_2", "P5A0", Sort_light_confs_h5);

   auto foo= [&data_disc_JJ_light_LMA21](string x) -> int {

     int j=-1;
     for(int i=0;i<(signed)data_disc_JJ_light_LMA21.Tag.size();i++) {
       if(x== data_disc_JJ_light_LMA21.Tag[i]) return i;
     }
     
     crash("you should not be here");
     return j;
   };
   
   //#######################        GOUNARIS SAKURAI   MODEL   ###################################


   //Gounaris Sakurai model

   if(! Skip_GS_analysis) {

   int npts_spline= 400;
   int Luscher_num_zeroes= 18;
   int Nresonances=16;
   //Init LL_functions;
   //find first  zeros of the Lusher functions
   Vfloat Luscher_zeroes;
   Zeta_function_zeroes(Luscher_num_zeroes, Luscher_zeroes);
   
   //############################################INTERPOLATE PHI FUNCTION AND DERIVATIVES#############################
   
   
   VVfloat phi_data, phi_der_data;
   Vfloat sx_int;
   Vfloat sx_der, dx_der;
   Vfloat Dz;
   
   for(int L_zero=0;L_zero<Nresonances+1;L_zero++) {
     cout<<"Computing n(Lusch): "<<L_zero<<endl;
     double sx, dx;
     //interpolating between the Luscher_zero[L_zero-1] and Luscher_zero[L_zero];
     if(L_zero==0) { sx_int.push_back(0.0); sx=0.0;}
     else {sx=Luscher_zeroes[L_zero-1];  sx_int.push_back(sx);}
     dx= Luscher_zeroes[L_zero];
     phi_data.resize(L_zero+1);
     phi_der_data.resize(L_zero+1);
     phi_data[L_zero].push_back(L_zero==0?0.0:-M_PI/2.0);
     //divide interval into thousand points;
     double dz = (dx-sx)/npts_spline;
     Dz.push_back(dz);
     
     
     for(int istep=1;istep<=npts_spline-1;istep++) { double pt= sx+dz*istep;
       phi_data[L_zero].push_back( phi(sqrt(pt)));}
     
     phi_data[L_zero].push_back(M_PI/2.0);
     double sx_der_loc =  phi_der_for_back(sqrt(sx)+1e-14, 1);
     double dx_der_loc =  phi_der_for_back(sqrt(dx)-1e-14, -1);
     sx_der.push_back(sx_der_loc);
     dx_der.push_back(dx_der_loc);
     
     phi_der_data[L_zero].push_back(sx_der_loc);
     for(int istep=1;istep<=npts_spline-1;istep++) { double pt= sx+dz*istep;
       phi_der_data[L_zero].push_back( phi_der(sqrt(pt)));}
     phi_der_data[L_zero].push_back(dx_der_loc);
     
   }

   
   
   
   
   
   LL_functions LL(phi_data,phi_der_data,sx_der, dx_der, sx_int, Dz, Nresonances, Luscher_zeroes);
   
   //###########################################END INTERPOLATION PHI FUNCTION AND DERIVATIVES################################
   cout<<"####Spline for phi(x) and phi'(x) successfully generated!"<<endl;
   
   double vol= 5.45*fm_to_inv_Gev;

   double Mpi= 0.135;
   double c = (0.100*0.100 - 0.135*0.135)/(a_B.ave()*a_B.ave());
   double BarM_B = sqrt( pow(Mpi,2) + 0.5*c*a_B_ave*a_B_ave*fm_to_inv_Gev*fm_to_inv_Gev  );
   double BarM_C = sqrt( pow(Mpi,2) + 0.5*c*a_C_ave*a_C_ave*fm_to_inv_Gev*fm_to_inv_Gev  ); 
   double BarM_D = sqrt( pow(Mpi,2) + 0.5*c*a_D_ave*a_D_ave*fm_to_inv_Gev*fm_to_inv_Gev  );
   double BarM_E = sqrt( pow(Mpi,2) + 0.5*c*a_E_ave*a_E_ave*fm_to_inv_Gev*fm_to_inv_Gev  ); 
   Vfloat En;
   LL.Find_pipi_energy_lev(vol, 0.770, 5.5, Mpi, 0.0, En);
   Vfloat EnB;
   LL.Find_pipi_energy_lev(vol, 0.770, 5.5, BarM_B , 0.0, EnB);
   Vfloat EnC;
   LL.Find_pipi_energy_lev(vol, 0.770, 5.5, BarM_C , 0.0, EnC);
   Vfloat EnD;
   LL.Find_pipi_energy_lev(vol, 0.770, 5.5, BarM_D , 0.0, EnD);
   Vfloat EnE;
   LL.Find_pipi_energy_lev(vol, 0.770, 5.5, BarM_E , 0.0, EnE);
   
  

   auto C = [&](double t) { return (10.0/9.0)*4.0*pow(alpha,2)*(1e10)*LL.V_pipi(t, vol, 0.770, 5.5, Mpi, 0.0,  En)*kernel_K(t,1);}; //t is in GeV^-1
   auto C_B = [&](double t) { return (10.0/9.0)*4.0*pow(alpha,2)*(1e10)*LL.V_pipi(t, vol, 0.770, 5.5, BarM_B, 0.0,  EnB)*kernel_K(t,1);}; //t is in GeV^-1
   auto C_C = [&](double t) { return (10.0/9.0)*4.0*pow(alpha,2)*(1e10)*LL.V_pipi(t, vol, 0.770, 5.5, BarM_C, 0.0,  EnC)*kernel_K(t,1);}; //t is in GeV^-1
   auto C_D = [&](double t) { return (10.0/9.0)*4.0*pow(alpha,2)*(1e10)*LL.V_pipi(t, vol, 0.770, 5.5, BarM_D, 0.0,  EnD)*kernel_K(t,1);}; //t is in GeV^-1
   auto C_E = [&](double t) { return (10.0/9.0)*4.0*pow(alpha,2)*(1e10)*LL.V_pipi(t, vol, 0.770, 5.5, BarM_E, 0.0,  EnE)*kernel_K(t,1);}; //t is in GeV^-1
  
  //get HVP
  gsl_function_pp<decltype(C)> F1(C);
  gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (10000);
  gsl_function *G1 = static_cast<gsl_function*>(&F1);
  double amu, amu_err;
  gsl_integration_qagiu(G1, 0.0, 0.0, 1e-6, 10000, w1, &amu, &amu_err);
  
  gsl_function_pp<decltype(C_B)> F_B(C_B);
  gsl_integration_workspace * w_B = gsl_integration_workspace_alloc (10000);
  gsl_function *G_B = static_cast<gsl_function*>(&F_B);
  double amu_B, amu_B_err;
  gsl_integration_qagiu(G_B, 0.0, 0.0, 1e-6, 10000, w_B, &amu_B, &amu_B_err);

  gsl_function_pp<decltype(C_C)> F_C(C_C);
  gsl_integration_workspace * w_C = gsl_integration_workspace_alloc (10000);
  gsl_function *G_C = static_cast<gsl_function*>(&F_C);
  double amu_C, amu_C_err;
  gsl_integration_qagiu(G_C, 0.0, 0.0, 1e-6, 10000, w_C, &amu_C, &amu_C_err);

  gsl_function_pp<decltype(C_D)> F_D(C_D);
  gsl_integration_workspace * w_D = gsl_integration_workspace_alloc (10000);
  gsl_function *G_D = static_cast<gsl_function*>(&F_D);
  double amu_D, amu_D_err;
  gsl_integration_qagiu(G_D, 0.0, 0.0, 1e-6, 10000, w_D, &amu_D, &amu_D_err);

  gsl_function_pp<decltype(C_E)> F_E(C_E);
  gsl_integration_workspace * w_E = gsl_integration_workspace_alloc (10000);
  gsl_function *G_E = static_cast<gsl_function*>(&F_E);
  double amu_E, amu_E_err;
  gsl_integration_qagiu(G_E, 0.0, 0.0, 1e-6, 10000, w_E, &amu_E, &amu_E_err);
  

  cout<<"B: "<< a_B_ave <<" "<<amu_B - amu<<endl;
  cout<<"C: "<< a_C_ave <<" "<<amu_C - amu<<endl;
  cout<<"D: "<< a_D_ave <<" "<<amu_D - amu<<endl;
  cout<<"E: "<< a_E_ave <<" "<<amu_E - amu<<endl;

   }
 

  // exit(-1);












   //##############################################################################################
 

   /*
   GaussianMersenne RN_TEST(41965);
   int T_SIZE=4000;
   //Read covariance matrix
   Eigen::MatrixXd Cov_TEST(T_SIZE, T_SIZE);

   for(int i=0;i<T_SIZE;i++) {
     for(int j=0;j<T_SIZE;j++) {

       Cov_TEST(i, j) = 0.04*exp( -1.0*abs(i-j)/10.0 );

     }
   }

   Eigen::VectorXd par_TEST(T_SIZE);
   for(int i=0;i<T_SIZE;i++) {
     par_TEST(i) = 1.0;
   }

   cout<<"GENERATED_COV_MATRIX_TEST"<<endl;


   Vfloat STREAM_TEST=Covariate(Cov_TEST, par_TEST, RN_TEST);

   cout<<"GENERATED_STREAM_TEST"<<endl;
   */

   //find B64 ensemble
   int id_B64_disco=-1;
   for(int id=0;id<(signed)data_disc_JJ_light.Tag.size();id++) {
     if(data_disc_JJ_light.Tag[id] == "cB211b.072.64") id_B64_disco=id;
   }
  
   CorrAnalysis Corr_BL(UseJack, Njacks,100);
   Corr_BL.Nt = 128;
   
   boost::filesystem::create_directory("../data/HVP_blinded");
   boost::filesystem::create_directory("../data/HVP_blinded/B64");
   boost::filesystem::create_directory("../data/HVP_blinded/B64/Bounding");
   boost::filesystem::create_directory("../data/HVP_blinded/B64/autocorr_Ct_conn");
   boost::filesystem::create_directory("../data/HVP_blinded/B64/autocorr_Ct_disc");
   boost::filesystem::create_directory("../data/HVP_blinded/B64/autocorr_amu_tc_conn");
   boost::filesystem::create_directory("../data/HVP_blinded/B64/autocorr_amu_tc_disc");
   boost::filesystem::create_directory("../data/HVP_blinded/B64/disco");
   boost::filesystem::create_directory("../data/HVP_blinded/C80");
   boost::filesystem::create_directory("../data/HVP_blinded/C80/Bounding");
   boost::filesystem::create_directory("../data/HVP_blinded/C80/disco");
   boost::filesystem::create_directory("../data/HVP_blinded/B96");
   boost::filesystem::create_directory("../data/HVP_blinded/B96/Bounding");
   boost::filesystem::create_directory("../data/HVP_blinded/B96/disco");
   boost::filesystem::create_directory("../data/HVP_blinded/D96");
   boost::filesystem::create_directory("../data/HVP_blinded/D96/Bounding");
   boost::filesystem::create_directory("../data/HVP_blinded/D96/disco");
  
   distr_t_list C_tm_BLIND= Corr_BL.corr_t( C_tm, "../data/HVP_blinded/B64/C_tm");
   distr_t_list C_OS_BLIND= Corr_BL.corr_t( C_OS, "../data/HVP_blinded/B64/C_OS");

   //evaluate partial integral of HVP on single conf

   VVfloat HVP_psum_tm(64), HVP_psum_OS(64);
   for(auto &H: HVP_psum_tm) H.resize(C_tm[0].size(),0.0);
   for(auto &H: HVP_psum_OS) H.resize(C_OS[0].size(),0.0);
   
   for(int iconf=0;iconf<(signed)C_tm[0].size(); iconf++) {

     double acc_tm=0; double acc_OS=0;
     for(int t=1;t<64;t++) {

       acc_tm += C_tm[t][iconf]*kernel_K(t, fm_to_inv_Gev*a_B_ave  );
       acc_OS += C_OS[t][iconf]*kernel_K(t, fm_to_inv_Gev*a_B_ave  );
       HVP_psum_tm[t][iconf] = acc_tm;
       HVP_psum_OS[t][iconf] = acc_OS;

     }

   }


   //evaluate autocorrelation for tm and OS
   for(int t=1;t<64;t++) {

     Compute_autocorrelation_time(C_tm[t],"../data/HVP_blinded/B64/autocorr_Ct_conn", "B64_tm_t_"+to_string(t));




     Compute_autocorrelation_time(C_OS[t],"../data/HVP_blinded/B64/autocorr_Ct_conn", "B64_OS_t_"+to_string(t));

     //Compute_autocorrelation_time(STREAM_TEST,"../data/HVP_blinded/B64/autocorr_Ct_conn", "TEST");

     Compute_autocorrelation_time( data_disc_JJ_light.col(0)[id_B64_disco][t],"../data/HVP_blinded/B64/autocorr_Ct_disc", "B64_t_"+to_string(t));
     
     Compute_autocorrelation_time(HVP_psum_tm[t],"../data/HVP_blinded/B64/autocorr_amu_tc_conn", "B64_tm_t_"+to_string(t));
     Compute_autocorrelation_time(HVP_psum_OS[t],"../data/HVP_blinded/B64/autocorr_amu_tc_conn", "B64_OS_t_"+to_string(t));

   }
   
   
   Corr_BL.Nt=80*2;
   
   distr_t_list C_tm_BLIND_C80= Corr_BL.corr_t( C_tm_C80, "../data/HVP_blinded/C80/C_tm");
   distr_t_list C_OS_BLIND_C80= Corr_BL.corr_t( C_OS_C80, "../data/HVP_blinded/C80/C_OS");

   
   distr_t_list C_tm_BLIND_C80_c1= Corr_BL.corr_t( C_tm_C80_c1, "../data/HVP_blinded/C80/C_tm_c1");
   distr_t_list C_OS_BLIND_C80_c1 = Corr_BL.corr_t( C_OS_C80_c1, "../data/HVP_blinded/C80/C_OS_c1");

   distr_t_list C_tm_BLIND_C80_c2= Corr_BL.corr_t( C_tm_C80_c2, "../data/HVP_blinded/C80/C_tm_c2");
   distr_t_list C_OS_BLIND_C80_c2 = Corr_BL.corr_t( C_OS_C80_c2, "../data/HVP_blinded/C80/C_OS_c2");

   distr_t_list C_tm_BLIND_C80_c3= Corr_BL.corr_t( C_tm_C80_c3, "../data/HVP_blinded/C80/C_tm_c3");
   distr_t_list C_OS_BLIND_C80_c3 = Corr_BL.corr_t( C_OS_C80_c3, "../data/HVP_blinded/C80/C_OS_c3");

   distr_t_list C_tm_BLIND_C80_c4= Corr_BL.corr_t( C_tm_C80_c4, "../data/HVP_blinded/C80/C_tm_c4");
   distr_t_list C_OS_BLIND_C80_c4 = Corr_BL.corr_t( C_OS_C80_c4, "../data/HVP_blinded/C80/C_OS_c4");



   Corr_BL.Nt=96*2;
   
   distr_t_list C_tm_BLIND_D96= Corr_BL.corr_t( C_tm_D96, "../data/HVP_blinded/D96/C_tm");
   distr_t_list C_OS_BLIND_D96= Corr_BL.corr_t( C_OS_D96, "../data/HVP_blinded/D96/C_OS");
   
   
   distr_t X_B64(UseJack,Njacks), X_C80(UseJack,Njacks), X_D96(UseJack,Njacks);

   for(int t=1;t<=30;t++) {
     
     X_B64 = X_B64 + ZA_B_ave*ZA_B_ave*C_tm_BLIND[t] - ZV_B_ave*ZV_B_ave*C_OS_BLIND[t];
     X_C80 = X_C80 + ZA_C_ave*ZA_C_ave*C_tm_BLIND_C80[t] - ZV_C_ave*ZV_C_ave*C_OS_BLIND_C80[t];
     X_D96 = X_D96 + ZA_D_ave*ZA_D_ave*C_tm_BLIND_D96[t] - ZV_D_ave*ZV_D_ave*C_OS_BLIND_D96[t];
   }
  

  

   
   
  C_tm_BLIND = 1e10*( qu*qu + qd*qd   )*ZA_B*ZA_B*C_tm_BLIND;
  C_OS_BLIND = 1e10*( qu*qu + qd*qd   )*ZV_B*ZV_B*C_OS_BLIND;
  
  C_tm_BLIND_C80 = 1e10*( qu*qu + qd*qd   )*ZA_C*ZA_C*C_tm_BLIND_C80;
  C_OS_BLIND_C80 = 1e10*( qu*qu + qd*qd   )*ZV_C*ZV_C*C_OS_BLIND_C80;
  
  C_tm_BLIND_D96 = 1e10*( qu*qu + qd*qd   )*ZA_D*ZA_D*C_tm_BLIND_D96;
  C_OS_BLIND_D96 = 1e10*( qu*qu + qd*qd   )*ZV_D*ZV_D*C_OS_BLIND_D96;
  

  Print_To_File( {} , {C_tm_BLIND.ave(), C_tm_BLIND.err(), C_OS_BLIND.ave(), C_OS_BLIND.err()}, "../data/HVP_blinded/B64/C_ren" , "" , ""); 
  Print_To_File( {} , {C_tm_BLIND_C80.ave(), C_tm_BLIND_C80.err(), C_OS_BLIND_C80.ave(), C_OS_BLIND_C80.err()}, "../data/HVP_blinded/C80/C_ren" , "" , "");
  Print_To_File( {} , {C_tm_BLIND_D96.ave(), C_tm_BLIND_D96.err(), C_OS_BLIND_D96.ave(), C_OS_BLIND_D96.err()}, "../data/HVP_blinded/D96/C_ren" , "" , "");

  
  
  
  distr_t p2_mot= 2*SQRT_D( a_B*a_B*(0.140*0.140) + pow( 2*M_PI/64,2));
  distr_t p2_mot_C80= 2*SQRT_D( a_C*a_C*(0.136*0.136) + pow( 2*M_PI/80,2));
  distr_t p2_mot_D96= 2*SQRT_D( a_D*a_D*(0.140*0.140) + pow( 2*M_PI/96,2));

  distr_t amu_HVP_tm, amu_HVP_OS;
  int Tcut_opt_tm, Tcut_opt_OS;

  
  Bounding_HVP(amu_HVP_tm, Tcut_opt_tm,  C_tm_BLIND, a_B,"../data/HVP_blinded/B64/Bounding/tm" , p2_mot);
  Bounding_HVP(amu_HVP_OS, Tcut_opt_OS,  C_OS_BLIND, a_B, "../data/HVP_blinded/B64/Bounding/OS" , p2_mot);
  
  cout<<"B64 tm-OS correlation: "<<(amu_HVP_tm%amu_HVP_OS)/(amu_HVP_tm.err()*amu_HVP_OS.err())<<endl;
  
  
  distr_t amu_HVP_tm_C80, amu_HVP_OS_C80;
  int Tcut_opt_tm_C80, Tcut_opt_OS_C80;

  Bounding_HVP(amu_HVP_tm_C80, Tcut_opt_tm_C80,  C_tm_BLIND_C80, a_C,"../data/HVP_blinded/C80/Bounding/tm" , p2_mot_C80);
  Bounding_HVP(amu_HVP_OS_C80, Tcut_opt_OS_C80,  C_OS_BLIND_C80, a_C, "../data/HVP_blinded/C80/Bounding/OS" , p2_mot_C80);

  cout<<"C80 tm-OS correlation: "<<(amu_HVP_tm_C80%amu_HVP_OS_C80)/(amu_HVP_tm_C80.err()*amu_HVP_OS_C80.err())<<endl;


  distr_t amu_HVP_tm_D96, amu_HVP_OS_D96;
  int Tcut_opt_tm_D96, Tcut_opt_OS_D96;

  Bounding_HVP(amu_HVP_tm_D96, Tcut_opt_tm_D96,  C_tm_BLIND_D96, a_D,"../data/HVP_blinded/D96/Bounding/tm" , p2_mot_D96);
  Bounding_HVP(amu_HVP_OS_D96, Tcut_opt_OS_D96,  C_OS_BLIND_D96, a_D, "../data/HVP_blinded/D96/Bounding/OS" , p2_mot_D96);

  cout<<"D96 tm-OS correlation: "<<(amu_HVP_tm_D96%amu_HVP_OS_D96)/(amu_HVP_tm_D96.err()*amu_HVP_OS_D96.err())<<endl;
  
  
  cout<<"B64:"<<endl;
  cout<<"amu(TM): "<<amu_HVP_tm.ave()<<" +- "<<amu_HVP_tm.err()<<endl;
  cout<<"amu(OS): "<<amu_HVP_OS.ave()<<" +- "<<amu_HVP_OS.err()<<endl;
  cout<<"average: "<<(0.5*amu_HVP_tm+0.5*amu_HVP_OS).ave()<<" +- "<<(0.5*amu_HVP_tm + 0.5*amu_HVP_OS).err()<<endl;
  cout<<"X: "<<X_B64.ave()<<" +- "<<X_B64.err()<<endl;
  cout<<"C80:"<<endl;
  cout<<"amu(TM): "<<amu_HVP_tm_C80.ave()<<" +- "<<amu_HVP_tm_C80.err()<<endl;
  cout<<"amu(OS): "<<amu_HVP_OS_C80.ave()<<" +- "<<amu_HVP_OS_C80.err()<<endl;
  cout<<"average: "<<(0.5*amu_HVP_tm_C80+0.5*amu_HVP_OS_C80).ave()<<" +- "<<(0.5*amu_HVP_tm_C80 + 0.5*amu_HVP_OS_C80).err()<<endl;
  cout<<"X: "<<X_C80.ave()<<" +- "<<X_C80.err()<<endl;
  cout<<"D96:"<<endl;
  cout<<"amu(TM): "<<amu_HVP_tm_D96.ave()<<" +- "<<amu_HVP_tm_D96.err()<<endl;
  cout<<"amu(OS): "<<amu_HVP_OS_D96.ave()<<" +- "<<amu_HVP_OS_D96.err()<<endl;
  cout<<"average: "<<(0.5*amu_HVP_tm_D96+0.5*amu_HVP_OS_D96).ave()<<" +- "<<(0.5*amu_HVP_tm_D96 + 0.5*amu_HVP_OS_D96).err()<<endl;
  cout<<"X: "<<X_D96.ave()<<" +- "<<X_D96.err()<<endl;


  //analyzing disconnected contribution

  int Nens=data_disc_JJ_light.Tag.size();

  for(int iens=0;iens<Nens;iens++) {


    int jens=-1;
    for(int b=0;b<(signed)data_P5P5_OS.Tag.size();b++) {
      if(data_P5P5_OS.Tag[b] == data_disc_JJ_light.Tag[iens]) jens=b;
    }
    assert(jens != -1);

    CorrAnalysis Corr(UseJack, Njacks,100);

    Corr.Nt = data_disc_JJ_light.nrows[iens];

    string Ens= data_disc_JJ_light.Tag[iens];
    string Ens_sh;
    distr_t Zv(UseJack), Za(UseJack), a_distr(UseJack);
    
    if(Ens=="cB211b.072.64") {
      Ens_sh="B64";
      Zv=ZV_B;
      Za=ZA_B;
      a_distr=a_B;
    }
    else if(Ens=="cB211b.072.96") {
      Ens_sh="B96";
      Zv=ZV_B;
      Za=ZA_B;
      a_distr=a_B;
    }
    else if(Ens=="cC211a.06.80") {
      Ens_sh="C80";
      Zv=ZV_C;
      Za=ZA_C;
      a_distr=a_C;
    }
    else if(Ens=="cD211a.054.96") {
      Ens_sh="D96";
      Zv=ZV_D;
      Za=ZA_D;
      a_distr=a_D;
    }
    else crash("Ensemble: "+Ens+" not found");

    

    //disco JJ
    distr_t_list Corr_disco_JJ_light_impr(UseJack);

  
    
    if(Ens=="cB211b.072.64")  Corr_disco_JJ_light_impr=  1e10*Zv*Zv*(qu+qd)*(qu+qd)*Corr.corr_t( data_disc_JJ_light_impr.col(0)[0], "../data/HVP_blinded/"+Ens_sh+"/disco/JJ_light_impr.dat");

    

    distr_t_list Corr_disco_JJ_light_impr_stoch(UseJack);
    if(Ens=="cB211b.072.64")  Corr_disco_JJ_light_impr_stoch=  1e10*Zv*Zv*(qu+qd)*(qu+qd)*Corr.corr_t( data_disc_JJ_light_impr_stoch.col(0)[0], "../data/HVP_blinded/"+Ens_sh+"/disco/JJ_light_impr_stoch.dat");
    //disco JJ LMA21

    
    distr_t_list Corr_disco_JJ_light_LMA21(UseJack), Corr_disco_JJ_light_impr_LMA21(UseJack);
    if(Ens=="cB211b.072.64" || Ens=="cD211a.054.96") {
      Corr_disco_JJ_light_LMA21=  1e10*Zv*Zv*(qu+qd)*(qu+qd)*Corr.corr_t( data_disc_JJ_light_LMA21.col(0)[foo(Ens)], "../data/HVP_blinded/"+Ens_sh+"/disco/JJ_light_LMA21.dat");
      Corr_disco_JJ_light_impr_LMA21=  1e10*Zv*Zv*(qu+qd)*(qu+qd)*Corr.corr_t( data_disc_JJ_light_impr_LMA21.col(0)[foo(Ens)], "../data/HVP_blinded/"+Ens_sh+"/disco/JJ_light_impr_LMA21.dat");

    }
    
    auto F_disco= [&](const Vfloat& par) { if((signed)par.size() != 2) crash("Lambda function F_disco expects par[2], but par["+to_string((signed)par.size())+"] provided"); return par[0] -par[1]*par[1];};


    
    
    //disco P5P5 LMA21
    distr_t_list P5P5_disc_LMA21(UseJack);
    distr_t_list P5P5_disc_impr_LMA21(UseJack);
    distr_t_list P5P5_disc_VEVs_LMA21(UseJack);
    distr_t_list P5P5_disc_impr_VEVs_LMA21(UseJack);
    //disco A0A0 LMA21
    distr_t_list A0A0_disc_LMA21(UseJack);
    distr_t_list A0A0_disc_impr_LMA21(UseJack);
    
    distr_t_list P5P5_conn = Corr.corr_t(data_P5P5_OS.col(0)[jens], "../data/HVP_blinded/"+Ens_sh+"/disco/P5P5_light_OS.dat");
    distr_t_list A0A0_conn = -1.0*Corr.corr_t(data_A0A0_OS.col(0)[jens], "../data/HVP_blinded/"+Ens_sh+"/disco/A0A0_light_OS.dat");
    distr_t_list A0P5_conn = Corr.corr_t(data_A0P5_OS.col(0)[jens], "../data/HVP_blinded/"+Ens_sh+"/disco/A0P5_light_OS.dat");
    

    //A0P5 and A0A0 correlators
    Corr.Reflection_sign=1;
    distr_t_list A0A0_disc = Corr.corr_t( data_disc_A0A0_light.col(0)[iens], "../data/HVP_blinded/"+Ens_sh+"/disco/A0A0.dat");
    Corr.Reflection_sign=-1;
    distr_t_list A0P5_disc = Corr.corr_t( data_disc_A0P5_light.col(0)[iens], "../data/HVP_blinded/"+Ens_sh+"/disco/A0P5.dat");
    Corr.Reflection_sign=1;

   
    
    int V = pow(Corr.Nt/2, 3);
   
    if(Ens=="cB211b.072.64" || Ens=="cD211a.054.96") {
      int Nconf_BUB= data_disc_P5_light_LMA21.col(0)[foo(Ens)][0].size();
      Vfloat BUB_P5_LMA21( Nconf_BUB,  0.0);
      
      P5P5_disc_LMA21=  Corr.corr_t( data_disc_P5_light_LMA21.col(0)[foo(Ens)], "../data/HVP_blinded/"+Ens_sh+"/disco/P5P5_light_LMA21.dat");
      P5P5_disc_impr_LMA21=  Corr.corr_t( data_disc_P5_light_impr_LMA21.col(0)[foo(Ens)], "../data/HVP_blinded/"+Ens_sh+"/disco/P5P5_light_impr_LMA21.dat");
      A0A0_disc_LMA21=  Corr.corr_t( data_disc_A0_light_LMA21.col(0)[foo(Ens)], "../data/HVP_blinded/"+Ens_sh+"/disco/A0A0_light_LMA21.dat");
      A0A0_disc_impr_LMA21=  Corr.corr_t( data_disc_A0_light_impr_LMA21.col(0)[foo(Ens)], "../data/HVP_blinded/"+Ens_sh+"/disco/A0A0_light_impr_LMA21.dat");
      
      for(int t=0;t<Corr.Nt;t++) 
	for(int iconf=0;iconf<Nconf_BUB;iconf++) BUB_P5_LMA21[iconf] += data_disc_bub_P5_light_LMA21.col(0)[foo(Ens)][t][iconf]/(sqrt(1.0*V)*Corr.Nt);

      
      Jackknife J(10000,Njacks);
      for(int t=0;t<Corr.Nt;t++) {
	P5P5_disc_VEVs_LMA21.distr_list.push_back(J.DoJack(F_disco, 2, data_disc_P5_light_LMA21.col(0)[foo(Ens)][t], BUB_P5_LMA21));
	P5P5_disc_impr_VEVs_LMA21.distr_list.push_back(J.DoJack(F_disco, 2, data_disc_P5_light_impr_LMA21.col(0)[foo(Ens)][t], BUB_P5_LMA21));
      }
      
    }

    distr_t_list P5P5_disc(UseJack);
    distr_t_list P5P5_disc_VEVs(UseJack);
    distr_t_list eta_disc = Corr.corr_t(data_disc_eta_light.col(0)[iens], "../data/HVP_blinded/"+Ens_sh+"/disco/eta_disc.dat");
    P5P5_disc = Corr.corr_t( data_disc_P5_light.col(0)[iens], "../data/HVP_blinded/"+Ens_sh+"/disco/P5P5_light.dat");
    int Nconf_BUB= data_disc_P5_light.col(0)[iens][0].size();
    Vfloat BUB_P5( Nconf_BUB,  0.0);
    for(int t=0;t<Corr.Nt;t++) 
      for(int iconf=0;iconf<Nconf_BUB;iconf++) BUB_P5[iconf] += data_disc_bub_P5_light.col(0)[iens][t][iconf]/(sqrt(1.0*V)*Corr.Nt);
    
    
    Jackknife J(10000,Njacks);
    for(int t=0;t<Corr.Nt;t++) {
      P5P5_disc_VEVs.distr_list.push_back(J.DoJack(F_disco, 2, data_disc_P5_light.col(0)[iens][t], BUB_P5));
    }

    distr_t_list P5P5_full = P5P5_conn + 2*P5P5_disc_VEVs;

    distr_t_list eta_full = P5P5_conn - 2*eta_disc;

    distr_t_list A0A0_full = Za*Za*(A0A0_conn - 2*A0A0_disc);

    //antysymmetrize
    for(int t=1;t<Corr.Nt/2;t++) A0P5_conn.distr_list[Corr.Nt/2 + t] = -1.0*A0P5_conn[Corr.Nt/2 + t];
    
    distr_t_list A0P5_full = A0P5_conn  +2*A0P5_disc;

    Print_To_File({}, {P5P5_disc_VEVs.ave(), P5P5_disc_VEVs.err(), P5P5_full.ave(), P5P5_full.err()}, "../data/HVP_blinded/"+Ens_sh+"/disco/P5P5_conn_&_disc.dat", "", "");
    Print_To_File({}, {A0A0_disc.ave(), A0A0_disc.err(), A0A0_full.ave(), A0A0_full.err()}, "../data/HVP_blinded/"+Ens_sh+"/disco/A0A0_conn_&_disc.dat", "", "");
    Print_To_File({}, {A0P5_disc.ave(), A0P5_disc.err(), A0P5_full.ave(), A0P5_full.err()}, "../data/HVP_blinded/"+Ens_sh+"/disco/A0P5_conn_&_disc.dat", "", "");
    Print_To_File({}, {eta_full.ave(), eta_full.err()}, "../data/HVP_blinded/"+Ens_sh+"/disco/eta_corr.dat", "", "");

    distr_t_list MPI_OS = Corr.effective_mass_t(P5P5_conn, "../data/HVP_blinded/"+Ens_sh+"/disco/Mpi_OS.dat");
    distr_t_list META=  Corr.effective_mass_t(eta_full, "../data/HVP_blinded/"+Ens_sh+"/disco/Meta.dat");
    
  
    

    if(Ens=="cB211b.072.64"  || Ens=="cD211a.054.96" ) {

      distr_t_list P5_full_LMA21= P5P5_conn + 2*P5P5_disc_VEVs_LMA21;
      distr_t_list P5_full_impr_LMA21=  P5P5_conn + 2*P5P5_disc_impr_VEVs_LMA21;
        
      distr_t_list Neutral_pion_LMA21 = Corr.effective_mass_t( P5_full_LMA21, "../data/HVP_blinded/"+Ens_sh+"/disco/Mpi_0_CD_LMA21.dat");
      distr_t_list Neutral_pion_impr_LMA21 = Corr.effective_mass_t( P5_full_impr_LMA21, "../data/HVP_blinded/"+Ens_sh+"/disco/Mpi_0_CD_impr_LMA21.dat");


      distr_t_list A0_full_LMA21= Za*Za*(A0A0_conn - 2*A0A0_disc_LMA21);
      distr_t_list A0_full_impr_LMA21=  Za*Za*(A0A0_conn - 2*A0A0_disc_impr_LMA21);
        
      distr_t_list Neutral_pion_A0_LMA21 = Corr.effective_mass_t( A0_full_LMA21, "../data/HVP_blinded/"+Ens_sh+"/disco/Mpi_0_CD_A0A0_LMA21.dat");
      distr_t_list Neutral_pion_A0_impr_LMA21 = Corr.effective_mass_t( A0_full_impr_LMA21, "../data/HVP_blinded/"+Ens_sh+"/disco/Mpi_0_CD_A0A0_impr_LMA21.dat");
    
    }

   
    
    double Mpic2 = 0.1408;
    if(Ens=="cC211a.06.80") Mpic2 = 0.1367;
    else if(Ens=="cB211b.072.64") Mpic2=0.1402;
    else if(Ens=="cD211a.054.96") Mpic2= 0.1408;
    else crash("CRASHHHH");

   
    
    
      distr_t_list Neutral_pion =  Corr.effective_mass_t( P5P5_disc_VEVs, "");
      distr_t_list Neutral_pion_CD = Corr.effective_mass_t( P5P5_full, "");
      distr_t_list Neutral_pion_A0A0 = Corr.effective_mass_t(A0A0_disc, "");
      distr_t_list Neutral_pion_A0A0_CD = Corr.effective_mass_t(A0A0_full, "");
      Corr.Reflection_sign=-1;
      distr_t_list Neutral_pion_A0P5 = Corr.effective_mass_t(A0P5_disc, "");
  
      distr_t_list Neutral_pion_A0P5_CD = Corr.effective_mass_t(A0P5_full, "");
  
      Corr.Reflection_sign=1;
      Print_To_File({}, { (Neutral_pion).ave(), (Neutral_pion).err()} , "../data/HVP_blinded/"+Ens_sh+"/disco/Mpi_0.dat.t", "", "");
      Print_To_File({}, { (Neutral_pion_CD).ave(), (Neutral_pion_CD).err()} , "../data/HVP_blinded/"+Ens_sh+"/disco/Mpi_0_CD.dat.t", "", "");
      Print_To_File({}, { (Neutral_pion_A0A0).ave(), (Neutral_pion_A0A0).err()} , "../data/HVP_blinded/"+Ens_sh+"/disco/Mpi_0_A0A0.dat.t", "", "");
      Print_To_File({}, { (Neutral_pion_A0A0_CD).ave(), (Neutral_pion_A0A0_CD).err()} , "../data/HVP_blinded/"+Ens_sh+"/disco/Mpi_0_CD_A0A0.dat.t", "", "");
      Print_To_File({}, { (Neutral_pion_A0P5).ave(), (Neutral_pion_A0P5).err()} , "../data/HVP_blinded/"+Ens_sh+"/disco/Mpi_0_A0P5.dat.t", "", "");
      Print_To_File({}, { (Neutral_pion_A0P5_CD).ave(), (Neutral_pion_A0P5_CD).err()} , "../data/HVP_blinded/"+Ens_sh+"/disco/Mpi_0_CD_A0P5.dat.t", "", "");
      distr_t_list C2 = (Neutral_pion_CD/a_distr)*(Neutral_pion_CD/a_distr) - Mpic2*Mpic2;
      distr_t_list C2_A0 = (Neutral_pion_A0A0_CD/a_distr)*(Neutral_pion_A0A0_CD/a_distr) - Mpic2*Mpic2;
      distr_t_list C2_AP = (Neutral_pion_A0P5_CD/a_distr)*(Neutral_pion_A0P5_CD/a_distr) - Mpic2*Mpic2;
      Print_To_File({}, { (C2/(a_distr*a_distr)).ave(), (C2/(a_distr*a_distr)).err()} , "../data/HVP_blinded/"+Ens_sh+"/disco/C2.dat.t", "", "");
      Print_To_File({}, { (C2_A0/(a_distr*a_distr)).ave(), (C2_A0/(a_distr*a_distr)).err()} , "../data/HVP_blinded/"+Ens_sh+"/disco/C2_A0.dat.t", "", "");
      Print_To_File({}, { (C2_AP/(a_distr*a_distr)).ave(), (C2_AP/(a_distr*a_distr)).err()} , "../data/HVP_blinded/"+Ens_sh+"/disco/C2_AP.dat.t", "", "");
      distr_t_list Corr_disco_JJ_light_stoch_LMA21(UseJack), Corr_disco_JJ_light_SET1_LMA21(UseJack);
      if(Ens=="cB211b.072.64"  || Ens=="cD211a.054.96") {
	Corr_disco_JJ_light_stoch_LMA21=  1e10*Zv*Zv*(qu+qd)*(qu+qd)*Corr.corr_t( data_disc_JJ_light_stoch_LMA21.col(0)[foo(Ens)], "../data/HVP_blinded/"+Ens_sh+"/disco/JJ_light_stoch_LMA21.dat");
	Corr_disco_JJ_light_impr_LMA21=  1e10*Zv*Zv*(qu+qd)*(qu+qd)*Corr.corr_t( data_disc_JJ_light_impr_LMA21.col(0)[foo(Ens)], "../data/HVP_blinded/"+Ens_sh+"/disco/JJ_light_impr_LMA21.dat");
      }
 

     
      distr_t_list Corr_disco_OS_JJ_light = 2*1e10*Zv*Zv*(qu*qu+qd*qd)*Corr.corr_t( data_disc_OS_JJ_light.col(0)[iens], "../data/HVP_blinded/"+Ens_sh+"/disco/JJ_OS_light.dat");
      distr_t_list Corr_disco_JJ_light=  1e10*Zv*Zv*(qu+qd)*(qu+qd)*Corr.corr_t( data_disc_JJ_light.col(0)[iens], "../data/HVP_blinded/"+Ens_sh+"/disco/JJ_light.dat");
      distr_t_list Corr_disco_JJ=  1e10*Zv*Zv*Corr.corr_t( data_disc_JJ.col(0)[iens], "../data/HVP_blinded/"+Ens_sh+"/disco/JJ.dat");
      distr_t_list Corr_disco_JJ_ls=  1e10*Zv*Zv*Corr.corr_t( data_disc_JJ_ls.col(0)[iens], "../data/HVP_blinded/"+Ens_sh+"/disco/JJ_ls.dat");
      int id_stoch=-1;
      for(int jens=0;jens<(signed)data_disc_JJ_light_stoch.Tag.size();jens++) {
	if(data_disc_JJ_light_stoch.Tag[jens]==Ens) id_stoch=jens;
      }
      
      distr_t_list Corr_disco_JJ_light_stoch(UseJack);
      if(id_stoch != -1)  Corr_disco_JJ_light_stoch=  1e10*Zv*Zv*(qu+qd)*(qu+qd)*Corr.corr_t( data_disc_JJ_light_stoch.col(0)[id_stoch], "../data/HVP_blinded/"+Ens_sh+"/disco/JJ_light_stoch.dat");
      
      //compute HVP

      distr_t_list HVP_light_tcut(UseJack);
      distr_t_list HVP_ls_tcut(UseJack);
      distr_t_list Win_light_tcut(UseJack);
      distr_t_list HVP_light_impr_tcut(UseJack);
      distr_t_list Win_light_impr_tcut(UseJack);
      distr_t_list HVP_light_impr_stoch_tcut(UseJack);
      distr_t_list Win_light_impr_stoch_tcut(UseJack);
      
      distr_t_list HVP_light_LMA21_tcut(UseJack);
      distr_t_list Win_light_LMA21_tcut(UseJack);
      distr_t_list HVP_light_stoch_LMA21_tcut(UseJack);
      distr_t_list Win_light_stoch_LMA21_tcut(UseJack);
      distr_t_list HVP_light_impr_LMA21_tcut(UseJack);
      distr_t_list Win_light_impr_LMA21_tcut(UseJack);
      
      distr_t_list HVP_light_stoch_tcut(UseJack);
      distr_t_list Win_light_stoch_tcut(UseJack);
      distr_t_list HVP_tcut(UseJack);
      distr_t_list Win_tcut(UseJack);

      distr_t_list HVP_OS_light_tcut(UseJack);
    
    distr_t p_sum_light= 0.0*Get_id_jack_distr(Njacks);
    distr_t p_sum_win_light= p_sum_light;

    distr_t p_sum_light_stoch= 0.0*Get_id_jack_distr(Njacks);
    distr_t p_sum_win_light_stoch= p_sum_light_stoch;

    distr_t p_sum_light_impr= 0.0*Get_id_jack_distr(Njacks);
    distr_t p_sum_win_light_impr= p_sum_light_impr;

    
    distr_t p_sum_light_impr_stoch= 0.0*Get_id_jack_distr(Njacks);
    distr_t p_sum_win_light_impr_stoch= p_sum_light_impr_stoch;

    distr_t p_sum_light_LMA21= 0.0*Get_id_jack_distr(Njacks);
    distr_t p_sum_win_light_LMA21= p_sum_light_LMA21;

    
    distr_t p_sum_light_stoch_LMA21= 0.0*Get_id_jack_distr(Njacks);
    distr_t p_sum_win_light_stoch_LMA21= p_sum_light_stoch_LMA21;

    distr_t p_sum_light_impr_LMA21= 0.0*Get_id_jack_distr(Njacks);
    distr_t p_sum_win_light_impr_LMA21= p_sum_light_impr_LMA21;

    distr_t p_sum_light_OS = 0.0*Get_id_jack_distr(Njacks);
  

    distr_t p_sum= 0.0*Get_id_jack_distr(Njacks);
    distr_t p_sum_win= p_sum;

    distr_t p_sum_ls= 0.0*Get_id_jack_distr(Njacks);
   

    auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
    distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);
    double t0=  0.4*fm_to_inv_Gev;
    double t1 = 1.0*fm_to_inv_Gev;
    double Delta= 0.15*fm_to_inv_Gev;
    auto th0 = [&t0, &Delta](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t0)/Delta));};
    auto th1 = [&t1, &Delta](double ta) ->double { return 1.0/(1.0 + exp(-2.0*(ta-t1)/Delta));};



    if(iens==0) {
      //rhad
      //load R_ratio at order alpha_s^4(1/tmin)

      Vfloat R_ratio_alpha4;
      Vfloat R_ratio_ergs_alpha4;
      
    	
      R_ratio_alpha4= Read_From_File("../../rhad/rhad/PROVA_RHAD", 3, 7);
      R_ratio_ergs_alpha4 =  Read_From_File("../../rhad/rhad/PROVA_RHAD", 0, 7);
	
      
          
      auto integrand_pt_cont = [ &th0, &R_ratio_alpha4, &R_ratio_ergs_alpha4](double t) {
	 //return (t<1e-10)?0.0:(4.0*pow(alpha,2)*( pow(qu,2)+ pow(qd,2))*free_vector_corr_cont(3, 0.0,t)*(1.0 - th0(t))*kernel_K(t,1.0));
	 double fact= 0.005*(1.0/(12*M_PI*M_PI));
	 double corr=0.0;
	 if(t<1e-10) return 0.0;
	 for(int j=0;j<(signed)R_ratio_ergs_alpha4.size();j++) {
	   double Erg= R_ratio_ergs_alpha4[j];
	   double R_rat= R_ratio_alpha4[j];
	   if( Erg > 0.2) corr+= fact*pow(Erg,2)*exp(-Erg*t)*R_rat;
	 }
	 return 4.0*pow(alpha,2)*1e10*corr*(1.0-th0(t))*kernel_K(t,1.0);
      };
	 
    
       

       
       double agm2_pert_up_to_t0, agm2_pert_up_to_t0_err;
       
       
      //set numerical integration with gsl
       
       gsl_function_pp<decltype(integrand_pt_cont)> Fp(integrand_pt_cont);
       gsl_integration_workspace * w_t0 = gsl_integration_workspace_alloc (10000);
       gsl_function *G = static_cast<gsl_function*>(&Fp);
       gsl_integration_qags(G, 0.0, 2.0*fm_to_inv_Gev, 0.0, 1e-8, 10000, w_t0, &agm2_pert_up_to_t0, &agm2_pert_up_to_t0_err);
       gsl_integration_workspace_free (w_t0);
       if(agm2_pert_up_to_t0_err/agm2_pert_up_to_t0 > 1e-5) crash("In determining continuum integral between 0 and tmin, accuracy achieved is only: "+to_string_with_precision(agm2_pert_up_to_t0_err/agm2_pert_up_to_t0,6));

       cout<<"PERT amu(s, tmin=0.198775fm): "<<agm2_pert_up_to_t0<<endl;
      
       
    }
    
    exit(-1);
     
    for(int t=1;t<Corr.Nt/2;t++) {

      p_sum_light = p_sum_light +  4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light.distr_list[t]*Ker.distr_list[t] ;
      p_sum_win_light = p_sum_win_light  +4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr) );
      p_sum_light_OS = p_sum_light_OS +  4.0*w(t,1)*pow(alpha,2)*Corr_disco_OS_JJ_light.distr_list[t]*Ker.distr_list[t] ;
      
      if(Ens == "cB211b.072.64") {

	p_sum_light_impr = p_sum_light_impr +  4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light_impr.distr_list[t]*Ker.distr_list[t] ;
	p_sum_win_light_impr = p_sum_win_light_impr  +4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light_impr.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr) );
	HVP_light_impr_tcut.distr_list.push_back(p_sum_light_impr);
	Win_light_impr_tcut.distr_list.push_back( p_sum_win_light_impr);

	p_sum_light_impr_stoch = p_sum_light_impr_stoch +  4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light_impr_stoch.distr_list[t]*Ker.distr_list[t] ;
	p_sum_win_light_impr_stoch = p_sum_win_light_impr_stoch  +4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light_impr_stoch.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr) );
	HVP_light_impr_stoch_tcut.distr_list.push_back(p_sum_light_impr_stoch);
	Win_light_impr_stoch_tcut.distr_list.push_back( p_sum_win_light_impr_stoch);


      }
      if(Ens == "cB211b.072.64"  || Ens=="cD211a.054.96") {

	//LMA21

	p_sum_light_LMA21 = p_sum_light_LMA21 +  4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light_LMA21.distr_list[t]*Ker.distr_list[t] ;
	p_sum_win_light_LMA21 = p_sum_win_light_LMA21  +4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light_LMA21.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr) );
	HVP_light_LMA21_tcut.distr_list.push_back(p_sum_light_LMA21);
	Win_light_LMA21_tcut.distr_list.push_back( p_sum_win_light_LMA21);



	p_sum_light_stoch_LMA21 = p_sum_light_stoch_LMA21 +  4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light_stoch_LMA21.distr_list[t]*Ker.distr_list[t] ;
	p_sum_win_light_stoch_LMA21 = p_sum_win_light_stoch_LMA21  +4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light_stoch_LMA21.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr) );
	HVP_light_stoch_LMA21_tcut.distr_list.push_back(p_sum_light_stoch_LMA21);
	Win_light_stoch_LMA21_tcut.distr_list.push_back( p_sum_win_light_stoch_LMA21);

	p_sum_light_impr_LMA21 = p_sum_light_impr_LMA21 +  4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light_impr_LMA21.distr_list[t]*Ker.distr_list[t] ;
	p_sum_win_light_impr_LMA21 = p_sum_win_light_impr_LMA21  +4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light_impr_LMA21.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr) );
	HVP_light_impr_LMA21_tcut.distr_list.push_back(p_sum_light_impr_LMA21);
	Win_light_impr_LMA21_tcut.distr_list.push_back( p_sum_win_light_impr_LMA21);



	
      }
      
      HVP_light_tcut.distr_list.push_back(p_sum_light);
      Win_light_tcut.distr_list.push_back(p_sum_win_light);
      HVP_OS_light_tcut.distr_list.push_back(p_sum_light_OS);

      if(id_stoch != -1) {
	p_sum_light_stoch = p_sum_light_stoch +  4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light_stoch.distr_list[t]*Ker.distr_list[t] ;
	p_sum_win_light_stoch = p_sum_win_light_stoch  +4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_light_stoch.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr) );
	HVP_light_stoch_tcut.distr_list.push_back(p_sum_light_stoch);
	Win_light_stoch_tcut.distr_list.push_back(p_sum_win_light_stoch);
      }
	
      p_sum = p_sum +  4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ.distr_list[t]*Ker.distr_list[t] ;
      p_sum_win = p_sum_win  +4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ.distr_list[t]*Ker.distr_list[t]*( distr_t::f_of_distr(th0, t*a_distr)  - distr_t::f_of_distr(th1, t*a_distr) );
      HVP_tcut.distr_list.push_back(p_sum);
      Win_tcut.distr_list.push_back(p_sum_win);

      p_sum_ls = p_sum_ls +  4.0*w(t,1)*pow(alpha,2)*Corr_disco_JJ_ls.distr_list[t]*Ker.distr_list[t] ;
      HVP_ls_tcut.distr_list.push_back(p_sum_ls);

    }

    Print_To_File({}, {HVP_light_tcut.ave(), HVP_light_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/HVP_light.psum", "", "");
    Print_To_File({}, {Win_light_tcut.ave(), Win_light_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/WIN_light.psum", "", "");
    Print_To_File({}, {HVP_OS_light_tcut.ave(), HVP_OS_light_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/HVP_OS_light.psum", "", "");

    Print_To_File({} , { (HVP_tcut-HVP_light_tcut).ave(), (HVP_tcut-HVP_light_tcut).err()}, "../data/HVP_blinded/"+Ens_sh+"/disco/HVP_s.psum", "", "");
    Print_To_File({} , { (HVP_tcut-HVP_ls_tcut).ave(), (HVP_tcut-HVP_ls_tcut).err()}, "../data/HVP_blinded/"+Ens_sh+"/disco/HVP_c.psum", "", "");
    
  
  

    if(Ens=="cB211b.072.64") {

      Print_To_File({}, {HVP_light_impr_tcut.ave(), HVP_light_impr_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/HVP_light_impr.psum", "", "");
      Print_To_File({}, {Win_light_impr_tcut.ave(), Win_light_impr_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/WIN_light_impr.psum", "", "");

      Print_To_File({}, {HVP_light_impr_stoch_tcut.ave(), HVP_light_impr_stoch_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/HVP_light_impr_stoch.psum", "", "");
      Print_To_File({}, {Win_light_impr_stoch_tcut.ave(), Win_light_impr_stoch_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/WIN_light_impr_stoch.psum", "", "");

    }

    if(Ens=="cB211b.072.64"  || Ens=="cD211a.054.96") {
      //LMA21

      //average Kyriacos and Simones runs

      distr_t_list HVP_light_ave(UseJack);

      for(int t=1; t < Corr.Nt/2 ; t++) {

	double w1= pow( HVP_light_tcut.err(t-1), -2)/( pow(HVP_light_tcut.err(t-1),-2) + pow( HVP_light_LMA21_tcut.err(t-1),-2));
	double w2= pow( HVP_light_LMA21_tcut.err(t-1), -2)/( pow(HVP_light_tcut.err(t-1),-2) + pow( HVP_light_LMA21_tcut.err(t-1),-2));

	HVP_light_ave.distr_list.push_back( w1*HVP_light_tcut.distr_list[t-1] + w2*HVP_light_LMA21_tcut.distr_list[t-1]);

      }
      
      Print_To_File({}, {HVP_light_LMA21_tcut.ave(), HVP_light_LMA21_tcut.err(), HVP_light_ave.ave(), HVP_light_ave.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/HVP_light_LMA21.psum", "", "");
      Print_To_File({}, {Win_light_LMA21_tcut.ave(), Win_light_LMA21_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/WIN_light_LMA21.psum", "", "");

      Print_To_File({}, {HVP_light_stoch_LMA21_tcut.ave(), HVP_light_stoch_LMA21_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/HVP_light_stoch_LMA21.psum", "", "");
      Print_To_File({}, {Win_light_stoch_LMA21_tcut.ave(), Win_light_stoch_LMA21_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/WIN_light_stoch_LMA21.psum", "", "");

      Print_To_File({}, {HVP_light_impr_LMA21_tcut.ave(), HVP_light_impr_LMA21_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/HVP_light_impr_LMA21.psum", "", "");
      Print_To_File({}, {Win_light_impr_LMA21_tcut.ave(), Win_light_impr_LMA21_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/WIN_light_impr_LMA21.psum", "", "");
      
    }

    Print_To_File({}, {HVP_tcut.ave(), HVP_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/HVP.psum", "", "");
    Print_To_File({}, {Win_tcut.ave(), Win_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/WIN.psum", "", "");

    if(id_stoch != -1) {

      Print_To_File({}, {HVP_light_stoch_tcut.ave(), HVP_light_stoch_tcut.err() , ((HVP_light_tcut-HVP_light_stoch_tcut)/HVP_light_tcut).ave(), ((HVP_light_tcut-HVP_light_stoch_tcut)/HVP_light_tcut).err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/HVP_light_stoch.psum", "", "");
      Print_To_File({}, {Win_light_stoch_tcut.ave(), Win_light_stoch_tcut.err()} ,  "../data/HVP_blinded/"+Ens_sh+"/disco/WIN_light_stoch.psum", "", "");
      
    }

   
  }






  return;

}
