#include "../include/RC_analysis.h"
#include "input.h"
#include "numerics.h"


const double alpha = 1.0/137.04;
const bool UseJack=1;
const int Njacks=120;
const int Nboots = 800;
const int Nmoms_A = 43;
const int Nmoms_B = 38;
const int Nmoms_C = 41;
const int Nmoms_D = 39;
const int Nmoms_E = 55;
const double fm_to_inv_Gev = 1.0 / 0.197327;
using namespace std;

void Perform_RC_fits() {



  //resample RCs
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack) , a_E(UseJack);

  
  double fmTGeV= 1.0/0.197327;

  LatticeInfo L_info_A, L_info_B, L_info_C, L_info_D, L_info_E;
  L_info_A.LatInfo_new_ens("cA211a.12.48");
  L_info_B.LatInfo_new_ens("cB211b.072.96");
  L_info_C.LatInfo_new_ens("cC211a.06.80");
  L_info_D.LatInfo_new_ens("cD211a.054.96");
  L_info_E.LatInfo_new_ens("cE211a.044.112");

  GaussianMersenne GM(78821);
  

  for(int ijack=0; ijack<Njacks;ijack++) {

  
    a_A.distr.push_back( L_info_A.a_from_afp*fmTGeV + GM()*L_info_A.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_B.distr.push_back( L_info_B.a_from_afp*fmTGeV + GM()*L_info_B.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_C.distr.push_back( L_info_C.a_from_afp*fmTGeV + GM()*L_info_C.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_D.distr.push_back( L_info_D.a_from_afp*fmTGeV + GM()*L_info_D.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
    a_E.distr.push_back( L_info_E.a_from_afp*fmTGeV + GM()*L_info_E.a_from_afp_err*fmTGeV/((UseJack==true)?sqrt(Njacks -1.0):1.0));
      

  }

  
  cout<<"RC generated!"<<endl;


  //read RCs

  distr_t_list Zv_A(UseJack), Za_A(UseJack), ZT_A(UseJack), Zq_A(UseJack), ZS_A(UseJack), ZP_A(UseJack);
  distr_t_list Zv_B(UseJack), Za_B(UseJack), ZT_B(UseJack), Zq_B(UseJack), ZS_B(UseJack), ZP_B(UseJack);
  distr_t_list Zv_C(UseJack), Za_C(UseJack), ZT_C(UseJack), Zq_C(UseJack), ZS_C(UseJack), ZP_C(UseJack);
  distr_t_list Zv_D(UseJack), Za_D(UseJack), ZT_D(UseJack), Zq_D(UseJack), ZS_D(UseJack), ZP_D(UseJack);
  distr_t_list Zv_E(UseJack), Za_E(UseJack), ZT_E(UseJack), Zq_E(UseJack), ZS_E(UseJack), ZP_E(UseJack);


  Vfloat p2_A, p2_B, p2_C, p2_D, p2_E;

  
  //A ensemble
  Vfloat p2_momA = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZV/RCs_Zv_ensA.dat", 0,2,1);
  Vfloat Zv_A_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZV/RCs_Zv_ensA.dat", 1,2,1);
  Vfloat Za_A_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZA/RCs_Za_ensA.dat", 1,2,1);
  Vfloat ZT_A_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZT/RCs_Zt_ensA.dat", 1,2,1);
  Vfloat Zq_A_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZQ/RCs_Zq_ensA.dat", 1,2,1);
  Vfloat ZS_A_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZS/RCs_Zs_ensA.dat", 1,2,1);
  Vfloat ZP_A_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZP/RCs_Zp_ensA.dat", 1,2,1);
  for(int imom=0;imom<Nmoms_A;imom++) {
    p2_A.push_back( p2_momA[imom*120] );
    Vfloat Zv_slice(Zv_A_mom_jack.begin() + imom*120, Zv_A_mom_jack.begin() + (imom+1)*120);
    Zv_A.distr_list.emplace_back( UseJack, Zv_slice);
    Vfloat Za_slice(Za_A_mom_jack.begin() + imom*120, Za_A_mom_jack.begin() + (imom+1)*120);
    Za_A.distr_list.emplace_back( UseJack, Za_slice);
    Vfloat ZT_slice(ZT_A_mom_jack.begin() + imom*120, ZT_A_mom_jack.begin() + (imom+1)*120);
    ZT_A.distr_list.emplace_back( UseJack, ZT_slice);
    Vfloat Zq_slice(Zq_A_mom_jack.begin() + imom*120, Zq_A_mom_jack.begin() + (imom+1)*120);
    Zq_A.distr_list.emplace_back( UseJack, Zq_slice);
    Vfloat ZS_slice(ZS_A_mom_jack.begin() + imom*120, ZS_A_mom_jack.begin() + (imom+1)*120);
    ZS_A.distr_list.emplace_back( UseJack, ZS_slice);
    Vfloat ZP_slice(ZP_A_mom_jack.begin() + imom*120, ZP_A_mom_jack.begin() + (imom+1)*120);
    ZP_A.distr_list.emplace_back( UseJack, ZP_slice);
  }
  //B ensemble
  Vfloat p2_momB = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZV/RCs_Zv_ensB.dat", 0,2,1);
  Vfloat Zv_B_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZV/RCs_Zv_ensB.dat", 1,2,1);
  Vfloat Za_B_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZA/RCs_Za_ensB.dat", 1,2,1);
  Vfloat ZT_B_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZT/RCs_Zt_ensB.dat", 1,2,1);
  Vfloat Zq_B_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZQ/RCs_Zq_ensB.dat", 1,2,1);
  Vfloat ZS_B_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZS/RCs_Zs_ensB.dat", 1,2,1);
  Vfloat ZP_B_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZP/RCs_Zp_ensB.dat", 1,2,1);
  for(int imom=0;imom<Nmoms_B;imom++) {
    p2_B.push_back( p2_momB[imom*120] );
    Vfloat Zv_slice(Zv_B_mom_jack.begin() + imom*120, Zv_B_mom_jack.begin() + (imom+1)*120);
    Zv_B.distr_list.emplace_back( UseJack, Zv_slice);
    Vfloat Za_slice(Za_B_mom_jack.begin() + imom*120, Za_B_mom_jack.begin() + (imom+1)*120);
    Za_B.distr_list.emplace_back( UseJack, Za_slice);
    Vfloat ZT_slice(ZT_B_mom_jack.begin() + imom*120, ZT_B_mom_jack.begin() + (imom+1)*120);
    ZT_B.distr_list.emplace_back( UseJack, ZT_slice);
    Vfloat Zq_slice(Zq_B_mom_jack.begin() + imom*120, Zq_B_mom_jack.begin() + (imom+1)*120);
    Zq_B.distr_list.emplace_back( UseJack, Zq_slice);
    Vfloat ZS_slice(ZS_B_mom_jack.begin() + imom*120, ZS_B_mom_jack.begin() + (imom+1)*120);
    ZS_B.distr_list.emplace_back( UseJack, ZS_slice);
    Vfloat ZP_slice(ZP_B_mom_jack.begin() + imom*120, ZP_B_mom_jack.begin() + (imom+1)*120);
    ZP_B.distr_list.emplace_back( UseJack, ZP_slice);
  }
  //C ensemble
  Vfloat p2_momC = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZV/RCs_Zv_ensC.dat", 0,2,1);
  Vfloat Zv_C_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZV/RCs_Zv_ensC.dat", 1,2,1);
  Vfloat Za_C_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZA/RCs_Za_ensC.dat", 1,2,1);
  Vfloat ZT_C_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZT/RCs_Zt_ensC.dat", 1,2,1);
  Vfloat Zq_C_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZQ/RCs_Zq_ensC.dat", 1,2,1);
  Vfloat ZS_C_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZS/RCs_Zs_ensC.dat", 1,2,1);
  Vfloat ZP_C_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZP/RCs_Zp_ensC.dat", 1,2,1);
  for(int imom=0;imom<Nmoms_C;imom++) {
    p2_C.push_back( p2_momC[imom*120] );
    Vfloat Zv_slice(Zv_C_mom_jack.begin() + imom*120, Zv_C_mom_jack.begin() + (imom+1)*120);
    Zv_C.distr_list.emplace_back( UseJack, Zv_slice);
    Vfloat Za_slice(Za_C_mom_jack.begin() + imom*120, Za_C_mom_jack.begin() + (imom+1)*120);
    Za_C.distr_list.emplace_back( UseJack, Za_slice);
    Vfloat ZT_slice(ZT_C_mom_jack.begin() + imom*120, ZT_C_mom_jack.begin() + (imom+1)*120);
    ZT_C.distr_list.emplace_back( UseJack, ZT_slice);
    Vfloat Zq_slice(Zq_C_mom_jack.begin() + imom*120, Zq_C_mom_jack.begin() + (imom+1)*120);
    Zq_C.distr_list.emplace_back( UseJack, Zq_slice);
    Vfloat ZS_slice(ZS_C_mom_jack.begin() + imom*120, ZS_C_mom_jack.begin() + (imom+1)*120);
    ZS_C.distr_list.emplace_back( UseJack, ZS_slice);
    Vfloat ZP_slice(ZP_C_mom_jack.begin() + imom*120, ZP_C_mom_jack.begin() + (imom+1)*120);
    ZP_C.distr_list.emplace_back( UseJack, ZP_slice);
  }
  //D ensemble
  Vfloat p2_momD = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZV/RCs_Zv_ensD.dat", 0,2,1);
  Vfloat Zv_D_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZV/RCs_Zv_ensD.dat", 1,2,1);
  Vfloat Za_D_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZA/RCs_Za_ensD.dat", 1,2,1);
  Vfloat ZT_D_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZT/RCs_Zt_ensD.dat", 1,2,1);
  Vfloat Zq_D_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZQ/RCs_Zq_ensD.dat", 1,2,1);
  Vfloat ZS_D_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZS/RCs_Zs_ensD.dat", 1,2,1);
  Vfloat ZP_D_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZP/RCs_Zp_ensD.dat", 1,2,1);
  for(int imom=0;imom<Nmoms_D;imom++) {
    p2_D.push_back( p2_momD[imom*120] );
    Vfloat Zv_slice(Zv_D_mom_jack.begin() + imom*120, Zv_D_mom_jack.begin() + (imom+1)*120);
    Zv_D.distr_list.emplace_back( UseJack, Zv_slice);
    Vfloat Za_slice(Za_D_mom_jack.begin() + imom*120, Za_D_mom_jack.begin() + (imom+1)*120);
    Za_D.distr_list.emplace_back( UseJack, Za_slice);
    Vfloat ZT_slice(ZT_D_mom_jack.begin() + imom*120, ZT_D_mom_jack.begin() + (imom+1)*120);
    ZT_D.distr_list.emplace_back( UseJack, ZT_slice);
    Vfloat Zq_slice(Zq_D_mom_jack.begin() + imom*120, Zq_D_mom_jack.begin() + (imom+1)*120);
    Zq_D.distr_list.emplace_back( UseJack, Zq_slice);
    Vfloat ZS_slice(ZS_D_mom_jack.begin() + imom*120, ZS_D_mom_jack.begin() + (imom+1)*120);
    ZS_D.distr_list.emplace_back( UseJack, ZS_slice);
    Vfloat ZP_slice(ZP_D_mom_jack.begin() + imom*120, ZP_D_mom_jack.begin() + (imom+1)*120);
    ZP_D.distr_list.emplace_back( UseJack, ZP_slice);
  }
  //E ensemble
  Vfloat p2_momE = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZV/RCs_Zv_ensE.dat", 0,2,1);
  Vfloat Zv_E_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZV/RCs_Zv_ensE.dat", 1,2,1);
  Vfloat Za_E_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZA/RCs_Za_ensE.dat", 1,2,1);
  Vfloat ZT_E_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZT/RCs_Zt_ensE.dat", 1,2,1);
  Vfloat Zq_E_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZQ/RCs_Zq_ensE.dat", 1,2,1);
  Vfloat ZS_E_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZS/RCs_Zs_ensE.dat", 1,2,1);
  Vfloat ZP_E_mom_jack = Read_From_File("../data_RC_Gregoris/RCs_after_evolution_21GeV2/ZP/RCs_Zp_ensE.dat", 1,2,1);
  for(int imom=0;imom<Nmoms_E;imom++) {
    p2_E.push_back( p2_momE[imom*120] );
    Vfloat Zv_slice(Zv_E_mom_jack.begin() + imom*120, Zv_E_mom_jack.begin() + (imom+1)*120);
    Zv_E.distr_list.emplace_back( UseJack, Zv_slice);
    Vfloat Za_slice(Za_E_mom_jack.begin() + imom*120, Za_E_mom_jack.begin() + (imom+1)*120);
    Za_E.distr_list.emplace_back( UseJack, Za_slice);
    Vfloat ZT_slice(ZT_E_mom_jack.begin() + imom*120, ZT_E_mom_jack.begin() + (imom+1)*120);
    ZT_E.distr_list.emplace_back( UseJack, ZT_slice);
    Vfloat Zq_slice(Zq_E_mom_jack.begin() + imom*120, Zq_E_mom_jack.begin() + (imom+1)*120);
    Zq_E.distr_list.emplace_back( UseJack, Zq_slice);
    Vfloat ZS_slice(ZS_E_mom_jack.begin() + imom*120, ZS_E_mom_jack.begin() + (imom+1)*120);
    ZS_E.distr_list.emplace_back( UseJack, ZS_slice);
    Vfloat ZP_slice(ZP_E_mom_jack.begin() + imom*120, ZP_E_mom_jack.begin() + (imom+1)*120);
    ZP_E.distr_list.emplace_back( UseJack, ZP_slice);
  }


  
  

  return;
}
