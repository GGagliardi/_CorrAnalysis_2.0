#include "../include/axial_WI_disco.h"
#include "Corr_analysis.h"
#include "numerics.h"
#include "random.h"


using namespace std;

const bool UseJack = true;
const int Njacks = 34;
const double ms = 0.01827820;
const double ml = 0.00072;
const double ml2 = 0.00144;
const double ml3 = 0.00216;
const double ml4 = 0.00360;
const double mls = (ms+ml);
const double mll2 = (ml2 + ml);
const double ml2s = (ms + ml2);
const double ml2l3 = (ml2 + ml3);
const double ml3s = (ms + ml3);
const double ml3l4 = (ml3 + ml4);
const double ml4s = (ms + ml4);

const double qs= 1.0/3.0;


const double alpha = 1.0 / 137.035999;
const double fm_to_inv_Gev= 1.0/0.197327;


void axial_WI_disco() {


  auto Sort_disco = [](string A, string B) {


    

			    int conf_length_A= A.length();
			    int conf_length_B= B.length();

			    int pos_a_slash=-1;
			    int pos_b_slash=-1;
			    for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
			    for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;

			    string A_bis = A.substr(pos_a_slash+1, 13);
			    string B_bis = B.substr(pos_b_slash+1, 13);

			    //cout<<"A: "<<A<<endl;
			    //cout<<"A_bis: "<<A_bis<<endl;
			    
			    if(A_bis != B_bis) return A_bis < B_bis;
			    else {

			      string A_new= A.substr(pos_a_slash+14, A.length()-1);
			      string B_new= B.substr(pos_b_slash+14, B.length()-1);

			      //cout<<"new A: "<<A_new<<endl;
			      
			      int pos_und_A=-1; int pos_und_B=-1;
			      for(int i=0;i<(signed)A_new.length();i++) if(A_new.substr(i,1)=="_") { pos_und_A=i; break;}
			      for(int j=0;j<(signed)B_new.length();j++) if(B_new.substr(j,1)=="_")  { pos_und_B=j; break;}
			      int A_copy=stoi(A_new.substr(0,pos_und_A));
			      int B_copy=stoi(B_new.substr(0,pos_und_B));

			      //cout<<"A_copy: "<<A_copy<<endl;

			      int A_hit= stoi(A_new.substr(pos_und_A+1, A_new.length()-1));
			      int B_hit = stoi(B_new.substr(pos_und_B+1, B_new.length()-1));

			      //cout<<"A_hit: "<<A_hit<<endl;


			      if(A_copy != B_copy) return A_copy < B_copy;
			      else {

				return A_hit < B_hit;

			      }
			    }
			    

			    crash("you should not be here");
			    return A < B;
				
  };




  //#######################################################################
  
  //naive
  data_t NAIVE_ll_R0_V1, NAIVE_ll_R1_V1, NAIVE_ss_R0_V1, NAIVE_ss_R1_V1;

  //OET non-optimized
  data_t OET_NO_ls_R0_V1, OET_NO_ls_R1_V1;


  //OET optimized
  data_t OET_ls_R0_V1, OET_ls_R1_V1;
  //further include different masses
  data_t OET_ll2_R0_V1, OET_ll2_R1_V1;
  data_t OET_l2s_R0_V1, OET_l2s_R1_V1;
  //third mass
  data_t OET_l2l3_R0_V1, OET_l2l3_R1_V1;
  data_t OET_l3s_R0_V1, OET_l3s_R1_V1;
  //fourth mass
  data_t OET_l3l4_R0_V1, OET_l3l4_R1_V1;
  data_t OET_l4s_R0_V1, OET_l4s_R1_V1;

 

  //NAIVE
  NAIVE_ll_R0_V1.Read("../disco_ls_WI_pushed", "mes_contr_NAIVE_ll_R0", "S0V1", Sort_disco);
  NAIVE_ll_R1_V1.Read("../disco_ls_WI_pushed", "mes_contr_NAIVE_ll_R1", "S0V1", Sort_disco);
  NAIVE_ss_R0_V1.Read("../disco_ls_WI_pushed", "mes_contr_NAIVE_ss_R0", "S0V1", Sort_disco);
  NAIVE_ss_R1_V1.Read("../disco_ls_WI_pushed", "mes_contr_NAIVE_ss_R1", "S0V1", Sort_disco);

  //OET non-optimized
  OET_NO_ls_R0_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_ls_R0", "P5A1", Sort_disco);
  OET_NO_ls_R1_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_ls_R1", "P5A1", Sort_disco);

  //OET optimized
  OET_ls_R0_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_ls_R0", "A1P5", Sort_disco);
  OET_ls_R1_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_ls_R1", "A1P5", Sort_disco);

  //OET second mass
  OET_l2s_R0_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_l2s_R0", "A1P5", Sort_disco);
  OET_l2s_R1_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_l2s_R1", "A1P5", Sort_disco);

  OET_ll2_R0_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_ll2_R0", "A1P5", Sort_disco);
  OET_ll2_R1_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_ll2_R1", "A1P5", Sort_disco);

  //OET third mass
  OET_l3s_R0_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_l3s_R0", "A1P5", Sort_disco);
  OET_l3s_R1_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_l3s_R1", "A1P5", Sort_disco);
  OET_l2l3_R0_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_l2l3_R0", "A1P5", Sort_disco);
  OET_l2l3_R1_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_l2l3_R1", "A1P5", Sort_disco);

  //OET fourth mass
  OET_l4s_R0_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_l4s_R0", "A1P5", Sort_disco);
  OET_l4s_R1_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_l4s_R1", "A1P5", Sort_disco);
  OET_l3l4_R0_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_l3l4_R0", "A1P5", Sort_disco);
  OET_l3l4_R1_V1.Read("../disco_ls_WI_pushed", "mes_contr_OET_l3l4_R1", "A1P5", Sort_disco);

  //#######################################################################

  
  //#######################################################################

  cout<<"V1 read"<<endl;
  
  //naive
  data_t NAIVE_ll_R0_V2, NAIVE_ll_R1_V2, NAIVE_ss_R0_V2, NAIVE_ss_R1_V2;

  //OET non-optimized
  data_t OET_NO_ls_R0_V2, OET_NO_ls_R1_V2;


  //OET optimized
  data_t OET_ls_R0_V2, OET_ls_R1_V2;
  //further include different masses
  data_t OET_ll2_R0_V2, OET_ll2_R1_V2;
  data_t OET_l2s_R0_V2, OET_l2s_R1_V2;
  //third mass
  data_t OET_l2l3_R0_V2, OET_l2l3_R1_V2;
  data_t OET_l3s_R0_V2, OET_l3s_R1_V2;
  //fourth mass
  data_t OET_l3l4_R0_V2, OET_l3l4_R1_V2;
  data_t OET_l4s_R0_V2, OET_l4s_R1_V2;

  //NAIVE
  NAIVE_ll_R0_V2.Read("../disco_ls_WI_pushed", "mes_contr_NAIVE_ll_R0", "S0V2", Sort_disco);
  NAIVE_ll_R1_V2.Read("../disco_ls_WI_pushed", "mes_contr_NAIVE_ll_R1", "S0V2", Sort_disco);
  NAIVE_ss_R0_V2.Read("../disco_ls_WI_pushed", "mes_contr_NAIVE_ss_R0", "S0V2", Sort_disco);
  NAIVE_ss_R1_V2.Read("../disco_ls_WI_pushed", "mes_contr_NAIVE_ss_R1", "S0V2",Sort_disco);

  //OET non-optimized
  OET_NO_ls_R0_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_ls_R0", "P5A2", Sort_disco);
  OET_NO_ls_R1_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_ls_R1", "P5A2", Sort_disco);

  //OET optimized
  OET_ls_R0_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_ls_R0", "A2P5", Sort_disco);
  OET_ls_R1_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_ls_R1", "A2P5", Sort_disco);

  //OET second mass
  OET_l2s_R0_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_l2s_R0", "A2P5", Sort_disco);
  OET_l2s_R1_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_l2s_R1", "A2P5", Sort_disco);

  OET_ll2_R0_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_ll2_R0", "A2P5", Sort_disco);
  OET_ll2_R1_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_ll2_R1", "A2P5", Sort_disco);

  //OET third mass
  OET_l3s_R0_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_l3s_R0", "A2P5",Sort_disco);
  OET_l3s_R1_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_l3s_R1", "A2P5", Sort_disco);
  OET_l2l3_R0_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_l2l3_R0", "A2P5",Sort_disco);
  OET_l2l3_R1_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_l2l3_R1", "A2P5",Sort_disco);

  //OET fourth mass
  OET_l4s_R0_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_l4s_R0", "A2P5",Sort_disco);
  OET_l4s_R1_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_l4s_R1", "A2P5",Sort_disco);
  OET_l3l4_R0_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_l3l4_R0", "A2P5",Sort_disco);
  OET_l3l4_R1_V2.Read("../disco_ls_WI_pushed", "mes_contr_OET_l3l4_R1", "A2P5",Sort_disco);

  //#######################################################################

  cout<<"V2 read"<<endl;
  
  //#######################################################################
  
  //naive
  data_t NAIVE_ll_R0_V3, NAIVE_ll_R1_V3, NAIVE_ss_R0_V3, NAIVE_ss_R1_V3;

  //OET non-optimized
  data_t OET_NO_ls_R0_V3, OET_NO_ls_R1_V3;


  //OET optimized
  data_t OET_ls_R0_V3, OET_ls_R1_V3;
  //further include different masses
  data_t OET_ll2_R0_V3, OET_ll2_R1_V3;
  data_t OET_l2s_R0_V3, OET_l2s_R1_V3;
  //third mass
  data_t OET_l2l3_R0_V3, OET_l2l3_R1_V3;
  data_t OET_l3s_R0_V3, OET_l3s_R1_V3;
  //fourth mass
  data_t OET_l3l4_R0_V3, OET_l3l4_R1_V3;
  data_t OET_l4s_R0_V3, OET_l4s_R1_V3;

  //NAIVE
  NAIVE_ll_R0_V3.Read("../disco_ls_WI_pushed", "mes_contr_NAIVE_ll_R0", "S0V3",Sort_disco);
  NAIVE_ll_R1_V3.Read("../disco_ls_WI_pushed", "mes_contr_NAIVE_ll_R1", "S0V3",Sort_disco);
  NAIVE_ss_R0_V3.Read("../disco_ls_WI_pushed", "mes_contr_NAIVE_ss_R0", "S0V3",Sort_disco);
  NAIVE_ss_R1_V3.Read("../disco_ls_WI_pushed", "mes_contr_NAIVE_ss_R1", "S0V3",Sort_disco);

  //OET non-optimized
  OET_NO_ls_R0_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_ls_R0", "P5A3",Sort_disco);
  OET_NO_ls_R1_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_ls_R1", "P5A3",Sort_disco);

  //OET optimized
  OET_ls_R0_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_ls_R0", "A3P5",Sort_disco);
  OET_ls_R1_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_ls_R1", "A3P5",Sort_disco);

  //OET second mass
  OET_l2s_R0_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_l2s_R0", "A3P5",Sort_disco);
  OET_l2s_R1_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_l2s_R1", "A3P5",Sort_disco);

  OET_ll2_R0_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_ll2_R0", "A3P5",Sort_disco);
  OET_ll2_R1_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_ll2_R1", "A3P5",Sort_disco);

  //OET third mass
  OET_l3s_R0_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_l3s_R0", "A3P5",Sort_disco);
  OET_l3s_R1_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_l3s_R1", "A3P5",Sort_disco);
  OET_l2l3_R0_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_l2l3_R0", "A3P5",Sort_disco);
  OET_l2l3_R1_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_l2l3_R1", "A3P5",Sort_disco);

  //OET fourth mass
  OET_l4s_R0_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_l4s_R0", "A3P5",Sort_disco);
  OET_l4s_R1_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_l4s_R1", "A3P5",Sort_disco);
  OET_l3l4_R0_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_l3l4_R0", "A3P5",Sort_disco);
  OET_l3l4_R1_V3.Read("../disco_ls_WI_pushed", "mes_contr_OET_l3l4_R1", "A3P5",Sort_disco);

  //#######################################################################

  cout<<"V3 read"<<endl;


  // NON OET2

  data_t NON_OET_V1, NON_OET_V2, NON_OET_V3;

  NON_OET_V1.Read("../disco_ls_WI_non_OET", "bub_J1_ls.txt", "");
  NON_OET_V2.Read("../disco_ls_WI_non_OET", "bub_J2_ls.txt", "");
  NON_OET_V3.Read("../disco_ls_WI_non_OET", "bub_J3_ls.txt", "");

  int T=128;


  /*
  
  //check std identiy
  int Nconfs=NAIVE_ll_R0_V1.col(0)[0][0].size();

  double naive_ave=0;double naive_err=0;
  double oet_ave=0; double oet_err=0;
  double oet2_ave=0.0; double oet2_err=0.0;

  for(int iconf=0;iconf<Nconfs;iconf++) { 

    double n_bub_re=0.0; double o_bub_re=0.0;
    double n_bub_im=0.0; double o_bub_im=0.0;

    double o_bub_2_im= 0.0;
    
    
    for(int t=0;t<T;t++) {

      
      
      n_bub_re += ( 0.0*NAIVE_ll_R0_V1.col(0)[0][t][iconf] + NAIVE_ll_R1_V1.col(0)[0][t][iconf] - NAIVE_ss_R0_V1.col(0)[0][t][iconf] - 0.0*NAIVE_ss_R1_V1.col(0)[0][t][iconf] );
      n_bub_im += 0.5*( 2.0*NAIVE_ll_R0_V1.col(1)[0][t][iconf] + 0.0*NAIVE_ll_R1_V1.col(1)[0][t][iconf] - 0.0*NAIVE_ss_R0_V1.col(1)[0][t][iconf] - 2.0*NAIVE_ss_R1_V1.col(1)[0][t][iconf] );

      o_bub_re += 0.5*mls*( OET_NO_ls_R0_V1.col(0)[0][t][iconf] - OET_NO_ls_R1_V1.col(0)[0][t][iconf]); 
      o_bub_im += 0.5*mls*( OET_NO_ls_R0_V1.col(1)[0][t][iconf] + OET_NO_ls_R1_V1.col(1)[0][t][iconf]);

      o_bub_2_im =  0.5*mll2*(-OET_ll2_R0_V1.col(1)[0][t][iconf]  -OET_ll2_R1_V1.col(1)[0][t][iconf]) +  0.5*ml2s*(OET_l2s_R0_V1.col(1)[0][t][iconf] + OET_l2s_R1_V1.col(1)[0][t][iconf]) ;

   
   
    
    }

  
    cout<<"iconf: "<<iconf<<" NAIVE["<<n_bub_re<<","<<n_bub_im<<"] OET_NO["<<o_bub_re<<","<<o_bub_im<<"]"<<endl;

    naive_ave += n_bub_im;
    naive_err += n_bub_im*n_bub_im;

    oet_ave += o_bub_im;
    oet_err += o_bub_im*o_bub_im;

    oet2_ave += o_bub_2_im;
    oet2_err += o_bub_2_im*o_bub_2_im;

    

    
    
  }
  

  */

  VVfloat BUB_NAIVE_V1 = summ_master( Multiply_Vvector_by_scalar(summ_master(NAIVE_ll_R0_V1.col(1)[0], NAIVE_ll_R1_V1.col(1)[0]),1.0), Multiply_Vvector_by_scalar( summ_master(NAIVE_ss_R0_V1.col(1)[0], NAIVE_ss_R1_V1.col(1)[0]), -0.0));
  VVfloat BUB_OET_V1 = Multiply_Vvector_by_scalar( summ_master(OET_ls_R0_V1.col(1)[0], OET_ls_R1_V1.col(1)[0]), mls);
  VVfloat BUB_OET2_V1 = summ_master( Multiply_Vvector_by_scalar( summ_master(OET_ll2_R0_V1.col(1)[0], OET_ll2_R1_V1.col(1)[0]), -1.0*mll2) , Multiply_Vvector_by_scalar( summ_master(OET_l2s_R0_V1.col(1)[0], OET_l2s_R1_V1.col(1)[0]), ml2s));
  VVfloat BUB_OET3_V1= summ_master( Multiply_Vvector_by_scalar( summ_master(OET_ll2_R0_V1.col(1)[0], OET_ll2_R1_V1.col(1)[0]), -1.0*mll2), Multiply_Vvector_by_scalar( summ_master( OET_l2l3_R0_V1.col(1)[0], OET_l2l3_R1_V1.col(1)[0]), -1.0*ml2l3) , Multiply_Vvector_by_scalar( summ_master( OET_l3s_R0_V1.col(1)[0], OET_l3s_R1_V1.col(1)[0]), ml3s));
  VVfloat BUB_OET4_V1= summ_master( Multiply_Vvector_by_scalar( summ_master(OET_ll2_R0_V1.col(1)[0], OET_ll2_R1_V1.col(1)[0]), -1.0*mll2), Multiply_Vvector_by_scalar( summ_master( OET_l2l3_R0_V1.col(1)[0], OET_l2l3_R1_V1.col(1)[0]), -1.0*ml2l3) ,   Multiply_Vvector_by_scalar( summ_master( OET_l3l4_R0_V1.col(1)[0], OET_l3l4_R1_V1.col(1)[0]), -1.0*ml3l4)          , Multiply_Vvector_by_scalar( summ_master( OET_l4s_R0_V1.col(1)[0], OET_l4s_R1_V1.col(1)[0]), ml4s));

  VVfloat BUB_NAIVE_V2 = summ_master(  Multiply_Vvector_by_scalar(summ_master(NAIVE_ll_R0_V2.col(1)[0], NAIVE_ll_R1_V2.col(1)[0]),1.0), Multiply_Vvector_by_scalar( summ_master(NAIVE_ss_R0_V2.col(1)[0], NAIVE_ss_R1_V2.col(1)[0]), 0.0));
  VVfloat BUB_OET_V2 = Multiply_Vvector_by_scalar( summ_master(OET_ls_R0_V2.col(1)[0], OET_ls_R1_V2.col(1)[0]), mls);
  VVfloat BUB_OET2_V2 = summ_master( Multiply_Vvector_by_scalar( summ_master(OET_ll2_R0_V2.col(1)[0], OET_ll2_R1_V2.col(1)[0]), -1.0*mll2) , Multiply_Vvector_by_scalar( summ_master(OET_l2s_R0_V2.col(1)[0], OET_l2s_R1_V2.col(1)[0]), ml2s));
  VVfloat BUB_OET3_V2= summ_master( Multiply_Vvector_by_scalar( summ_master(OET_ll2_R0_V2.col(1)[0], OET_ll2_R1_V2.col(1)[0]), -1.0*mll2), Multiply_Vvector_by_scalar( summ_master( OET_l2l3_R0_V2.col(1)[0], OET_l2l3_R1_V2.col(1)[0]), -1.0*ml2l3) , Multiply_Vvector_by_scalar( summ_master( OET_l3s_R0_V2.col(1)[0], OET_l3s_R1_V2.col(1)[0]), ml3s));
  VVfloat BUB_OET4_V2= summ_master( Multiply_Vvector_by_scalar( summ_master(OET_ll2_R0_V2.col(1)[0], OET_ll2_R1_V2.col(1)[0]), -1.0*mll2), Multiply_Vvector_by_scalar( summ_master( OET_l2l3_R0_V2.col(1)[0], OET_l2l3_R1_V2.col(1)[0]), -1.0*ml2l3) ,   Multiply_Vvector_by_scalar( summ_master( OET_l3l4_R0_V2.col(1)[0], OET_l3l4_R1_V2.col(1)[0]), -1.0*ml3l4)          , Multiply_Vvector_by_scalar( summ_master( OET_l4s_R0_V2.col(1)[0], OET_l4s_R1_V2.col(1)[0]), ml4s));

  
  VVfloat BUB_NAIVE_V3 = summ_master(  Multiply_Vvector_by_scalar(summ_master(NAIVE_ll_R0_V3.col(1)[0], NAIVE_ll_R1_V3.col(1)[0]),1.0), Multiply_Vvector_by_scalar( summ_master(NAIVE_ss_R0_V3.col(1)[0], NAIVE_ss_R1_V3.col(1)[0]), -0.0));
  VVfloat BUB_OET_V3 = Multiply_Vvector_by_scalar( summ_master(OET_ls_R0_V3.col(1)[0], OET_ls_R1_V3.col(1)[0]), mls);
  VVfloat BUB_OET2_V3 = summ_master( Multiply_Vvector_by_scalar( summ_master(OET_ll2_R0_V3.col(1)[0], OET_ll2_R1_V3.col(1)[0]), -1.0*mll2) , Multiply_Vvector_by_scalar( summ_master(OET_l2s_R0_V3.col(1)[0], OET_l2s_R1_V3.col(1)[0]), ml2s));
  VVfloat BUB_OET3_V3= summ_master( Multiply_Vvector_by_scalar( summ_master(OET_ll2_R0_V3.col(1)[0], OET_ll2_R1_V3.col(1)[0]), -1.0*mll2), Multiply_Vvector_by_scalar( summ_master( OET_l2l3_R0_V3.col(1)[0], OET_l2l3_R1_V3.col(1)[0]), -1.0*ml2l3) , Multiply_Vvector_by_scalar( summ_master( OET_l3s_R0_V3.col(1)[0], OET_l3s_R1_V3.col(1)[0]), ml3s));
  VVfloat BUB_OET4_V3= summ_master( Multiply_Vvector_by_scalar( summ_master(OET_ll2_R0_V3.col(1)[0], OET_ll2_R1_V3.col(1)[0]), -1.0*mll2), Multiply_Vvector_by_scalar( summ_master( OET_l2l3_R0_V3.col(1)[0], OET_l2l3_R1_V3.col(1)[0]), -1.0*ml2l3) ,   Multiply_Vvector_by_scalar( summ_master( OET_l3l4_R0_V3.col(1)[0], OET_l3l4_R1_V3.col(1)[0]), -1.0*ml3l4)          , Multiply_Vvector_by_scalar( summ_master( OET_l4s_R0_V3.col(1)[0], OET_l4s_R1_V3.col(1)[0]), ml4s));



  VVfloat BUB_NON_OET_V1 = NON_OET_V1.col(1)[0];
  VVfloat BUB_NON_OET_V2 = NON_OET_V2.col(1)[0];
  VVfloat BUB_NON_OET_V3 = NON_OET_V3.col(1)[0];
  
  
  //separate stoch sources from gauge

  cout<<"bubbles computed"<<endl;

  int Nc=34;

  int TOT= BUB_NAIVE_V1[0].size();

  int Nhits= TOT/Nc;

  cout<<"Nconfs: "<<Nc<<" Nhits: "<<Nhits<<endl;

  assert( TOT%Nc == 0);

 
  VVVfloat BUB_NAIVE_V1_resh(Nc), BUB_OET_V1_resh(Nc), BUB_OET2_V1_resh(Nc), BUB_OET3_V1_resh(Nc), BUB_OET4_V1_resh(Nc);
  VVVfloat BUB_NAIVE_V2_resh(Nc), BUB_OET_V2_resh(Nc), BUB_OET2_V2_resh(Nc), BUB_OET3_V2_resh(Nc), BUB_OET4_V2_resh(Nc);
  VVVfloat BUB_NAIVE_V3_resh(Nc), BUB_OET_V3_resh(Nc), BUB_OET2_V3_resh(Nc), BUB_OET3_V3_resh(Nc), BUB_OET4_V3_resh(Nc);

  
 

  for(int c=0;c<Nc;c++) {

    BUB_NAIVE_V1_resh[c].resize(Nhits);
    BUB_OET_V1_resh[c].resize(Nhits);
    BUB_OET2_V1_resh[c].resize(Nhits);
    BUB_OET3_V1_resh[c].resize(Nhits);
    BUB_OET4_V1_resh[c].resize(Nhits);
    
    BUB_NAIVE_V2_resh[c].resize(Nhits);
    BUB_OET_V2_resh[c].resize(Nhits);
    BUB_OET2_V2_resh[c].resize(Nhits);
    BUB_OET3_V2_resh[c].resize(Nhits);
    BUB_OET4_V2_resh[c].resize(Nhits);

    BUB_NAIVE_V3_resh[c].resize(Nhits);
    BUB_OET_V3_resh[c].resize(Nhits);
    BUB_OET2_V3_resh[c].resize(Nhits);
    BUB_OET3_V3_resh[c].resize(Nhits);
    BUB_OET4_V3_resh[c].resize(Nhits);
    
    
    for(int ihit=0;ihit<Nhits;ihit++) {

      
      for(int t=0;t<T;t++) {
      
	BUB_NAIVE_V1_resh[c][ihit].push_back( BUB_NAIVE_V1[t][ ihit + c*Nhits]);
	BUB_OET_V1_resh[c][ihit].push_back( BUB_OET_V1[t][ ihit + c*Nhits]);
	BUB_OET2_V1_resh[c][ihit].push_back( BUB_OET2_V1[t][ ihit + c*Nhits]);
	BUB_OET3_V1_resh[c][ihit].push_back( BUB_OET3_V1[t][ ihit + c*Nhits]);
	BUB_OET4_V1_resh[c][ihit].push_back( BUB_OET4_V1[t][ ihit + c*Nhits]);

	BUB_NAIVE_V2_resh[c][ihit].push_back( BUB_NAIVE_V2[t][ ihit + c*Nhits]);
	BUB_OET_V2_resh[c][ihit].push_back( BUB_OET_V2[t][ ihit + c*Nhits]);
	BUB_OET2_V2_resh[c][ihit].push_back( BUB_OET2_V2[t][ ihit + c*Nhits]);
	BUB_OET3_V2_resh[c][ihit].push_back( BUB_OET3_V2[t][ ihit + c*Nhits]);
	BUB_OET4_V2_resh[c][ihit].push_back( BUB_OET4_V2[t][ ihit + c*Nhits]);

	BUB_NAIVE_V3_resh[c][ihit].push_back( BUB_NAIVE_V3[t][ ihit + c*Nhits]);
	BUB_OET_V3_resh[c][ihit].push_back( BUB_OET_V3[t][ ihit + c*Nhits]);
	BUB_OET2_V3_resh[c][ihit].push_back( BUB_OET2_V3[t][ ihit + c*Nhits]);
	BUB_OET3_V3_resh[c][ihit].push_back( BUB_OET3_V3[t][ ihit + c*Nhits]);
	BUB_OET4_V3_resh[c][ihit].push_back( BUB_OET4_V3[t][ ihit + c*Nhits]);
	

      }

    }
  }


  cout<<"bubble reshuffled"<<endl;


  //Read file containing source pos info

  VVint tsou_list(Nc);

  vector<vector<string>> Confs(Nc);

  ifstream Read_sou("../disco_ls_WI_pushed/cB211b.072.64/tsou_list");

  for(int ic=0;ic<Nc;ic++) {
    for(int ihit=0;ihit<Nhits;ihit++) {
      string conf;
      int tsou;
      Read_sou>>conf>>tsou;

      tsou_list[ic].push_back(tsou);
      Confs[ic].push_back(conf);
      if(Read_sou.eof()) crash("Reached EOF before file finished");

    }

  }

  for(int ic=0;ic<Nc;ic++) {
    for(int ihit=0;ihit<Nhits;ihit++) {
      cout<<Confs[ic][ihit]<<", source located at time t="<<tsou_list[ic][ihit]<<endl;
    }
  }

  cout<<"tsource file read"<<endl;

  VVfloat VV_NAIVE_V1 = convolute(BUB_NAIVE_V1_resh, BUB_NAIVE_V1_resh, 1, tsou_list);
  VVfloat VV_OET_V1 = convolute(BUB_OET_V1_resh, BUB_OET_V1_resh, 1,  tsou_list);
  VVfloat VV_OET2_V1 = convolute(BUB_OET2_V1_resh, BUB_OET2_V1_resh, 1,  tsou_list);
  VVfloat VV_OET3_V1 = convolute(BUB_OET3_V1_resh, BUB_OET3_V1_resh, 1,  tsou_list);
  VVfloat VV_OET4_V1 = convolute(BUB_OET4_V1_resh, BUB_OET4_V1_resh, 1,  tsou_list);
    
  VVfloat VV_NAIVE_V2 = convolute(BUB_NAIVE_V2_resh, BUB_NAIVE_V2_resh, 1,  tsou_list);
  VVfloat VV_OET_V2 = convolute(BUB_OET_V2_resh, BUB_OET_V2_resh, 1,  tsou_list);
  VVfloat VV_OET2_V2 = convolute(BUB_OET2_V2_resh, BUB_OET2_V2_resh, 1,  tsou_list);
  VVfloat VV_OET3_V2 = convolute(BUB_OET3_V2_resh, BUB_OET3_V2_resh, 1,  tsou_list);
  VVfloat VV_OET4_V2 = convolute(BUB_OET4_V2_resh, BUB_OET4_V2_resh, 1,  tsou_list);

  
  VVfloat VV_NAIVE_V3 = convolute(BUB_NAIVE_V3_resh, BUB_NAIVE_V3_resh, 1,  tsou_list);
  VVfloat VV_OET_V3 = convolute(BUB_OET_V3_resh, BUB_OET_V3_resh, 1,  tsou_list);
  VVfloat VV_OET2_V3 = convolute(BUB_OET2_V3_resh, BUB_OET2_V3_resh, 1,  tsou_list);
  VVfloat VV_OET3_V3 = convolute(BUB_OET3_V3_resh, BUB_OET3_V3_resh, 1,  tsou_list);
  VVfloat VV_OET4_V3 = convolute(BUB_OET4_V3_resh, BUB_OET4_V3_resh, 1,  tsou_list);


  VVfloat VV_NON_OET_V1 = convolute(BUB_NON_OET_V1, BUB_NON_OET_V1,1);
  VVfloat VV_NON_OET_V2 = convolute(BUB_NON_OET_V2, BUB_NON_OET_V2,1);
  VVfloat VV_NON_OET_V3 = convolute(BUB_NON_OET_V3, BUB_NON_OET_V3,1);
  

  CorrAnalysis Corr(UseJack,Njacks,100);;
  Corr.Perform_Nt_t_average=0;
  Corr.Nt=128;

  

  double L = 64.0;
  double F= pow(L,3)*pow(T,2)*(qs/2)*(qs/2);
  double F_NOET= 1.0/pow(L,3);

  boost::filesystem::create_directory("../data/axial_WI_disco");


  VVfloat TEST_OET=Multiply_Vvector_by_scalar(summ_master(VV_OET_V1, VV_OET_V2, VV_OET_V3),F/3.0);

  for(int c=0;c<Confs.size(); c++) {

    cout<<"conf: "<<Confs[c][0]<<endl;
    for(int t=0;t<T;t++) {
      cout<<"t: "<<t<<" "<<TEST_OET[t][c]<<endl;
    }
    
  }
  

  distr_t_list V_impr_NAIVE= Corr.corr_t(Multiply_Vvector_by_scalar(summ_master(VV_NAIVE_V1, VV_NAIVE_V2, VV_NAIVE_V3),F/3.0), "../data/axial_WI_disco/V_NAIVE_impr");
  distr_t_list V_impr_OET= Corr.corr_t(Multiply_Vvector_by_scalar(summ_master(VV_OET_V1, VV_OET_V2, VV_OET_V3),F/3.0), "../data/axial_WI_disco/V_OET_impr");
  distr_t_list V_impr_OET2= Corr.corr_t(Multiply_Vvector_by_scalar(summ_master(VV_OET2_V1, VV_OET2_V2, VV_OET2_V3),F/3.0), "../data/axial_WI_disco/V_OET2_impr");
  distr_t_list V_impr_OET3= Corr.corr_t(Multiply_Vvector_by_scalar( summ_master(VV_OET3_V1, VV_OET3_V2, VV_OET3_V3),F/3.0), "../data/axial_WI_disco/V_OET3_impr");
  distr_t_list V_impr_OET4= Corr.corr_t(Multiply_Vvector_by_scalar(summ_master(VV_OET4_V1, VV_OET4_V2, VV_OET4_V3),F/3.0), "../data/axial_WI_disco/V_OET4_impr");

  distr_t_list V_non_OET =  Corr.corr_t(Multiply_Vvector_by_scalar(summ_master(VV_NON_OET_V1, VV_NON_OET_V2, VV_NON_OET_V3),F_NOET/3.0), "../data/axial_WI_disco/V_non_OET");




  distr_t a_distr(UseJack);

  double Zv= 0.70637654;

  GaussianMersenne GM(43112);

  double aB_ave= 0.07948*fm_to_inv_Gev;
  double aB_err =0.00011*fm_to_inv_Gev;

  for(int ijack=0;ijack<Njacks;ijack++) a_distr.distr.push_back( aB_ave + GM()*aB_err/sqrt(Njacks-1.0));

  //evaluate amu
  auto K = [&](double Mv, double t, double size) -> double { return kernel_K(t, Mv);};
  distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);


  distr_t_list amu_OET(UseJack), amu_NON_OET(UseJack);


  distr_t p_sum_OET= 0.0*Get_id_jack_distr(Njacks);
  distr_t p_sum_non_OET = 0.0*Get_id_jack_distr(Njacks);


  for(int t=0;t<T/2;t++) {

    p_sum_OET = p_sum_OET + 4.0*Zv*Zv*1e10*w(t,1)*pow(alpha,2)*V_impr_OET4.distr_list[t]*Ker.distr_list[t];

    p_sum_non_OET = p_sum_non_OET + 4.0*Zv*Zv*1e10*w(t,1)*pow(alpha,2)*V_non_OET.distr_list[t]*Ker.distr_list[t];

    amu_OET.distr_list.push_back( p_sum_OET);

    amu_NON_OET.distr_list.push_back( p_sum_non_OET);


  }


  Print_To_File({}, {amu_OET.ave(), amu_OET.err()}, "../data/axial_WI_disco/amu_OET","","");
  Print_To_File({}, {amu_NON_OET.ave(), amu_NON_OET.err()}, "../data/axial_WI_disco/amu_non_OET","","");
  
  
  
  /*
  cout<<"NAIVE: "<<naive_ave/Nconfs<<" +- "<< sqrt( naive_err/(Nconfs*(Nconfs-1)) - (naive_ave*naive_ave/(Nconfs*Nconfs))/(Nconfs-1))<<endl;
  cout<<"OET: "<<oet_ave/Nconfs<<" +- "<< sqrt( oet_err/(Nconfs*(Nconfs-1)) - (oet_ave*oet_ave/(Nconfs*Nconfs))/(Nconfs-1))<<endl;
  cout<<"OET2: "<<oet2_ave/Nconfs<<" +- "<< sqrt( oet2_err/(Nconfs*(Nconfs-1)) - (oet2_ave*oet2_ave/(Nconfs*Nconfs))/(Nconfs-1))<<endl;
  */
  

  
  
    
  



  return;

}
