#ifndef __PINGU__
#define __PINGU__


#include "Meson_mass_extrapolation.h"
#include "Bootstrap_fit.h"
#include "scale_setting.h"
#include "stat.h"
#include "binary_io.h"
#include "RC_WI_analysis.h"
#include "numerics.h"
#include "Corr_analysis.h"
#include "complex_stat.h"
#include "highPrec.h"
#include "Spectral.h"


using namespace std;





class FERM4_T {

public:
  FERM4_T() {}
  FERM4_T(const vector<int>& s) : strp(s) , dims(s.size()) {}
  FERM4_T(const vector<int>& s, int N) :  T_list(N), strp(s) , dims(s.size()) {}
  

  int size() const { return T_list.size() ; };
  int coord(int pos, int dir) const {
    int out=pos;
    for(int id=0;id<dir;id++)  out /= strp[id];
    return out%strp[dir];
  }
 
  int pos(const vector<int>& p) {
    int r=0; int N=1;
    for(int i=0;i<dims;i++) { assert(p[i] < strp[i]);   r += p[i]*N; N*=strp[i];}
    return r;
  }
  int t_size() const { return T_list[0].size(); }

  void free_mem() { dims=0; strp.clear(), T_list.clear(); T_list.shrink_to_fit(); } 

  vector<complex_distr_t_list> T_list;
  vector<int> strp;
  int dims;
};


void summ_FERM4_T(const FERM4_T &in, const FERM4_T& b) ;
     



void FERM4_T_gprod(const FERM4_T& in, FERM4_T& out, const C_MATRIX& g, int str_id) ;

complex_distr_t_list glb_red_FERM4_T(const FERM4_T &A, const FERM4_T &B,
                                     const vector<pair<int, int>> &rules);

complex_distr_t_list glb_red_FERM4_T( const FERM4_T &A,  const vector<pair<int,int>>& rules ); 


class NISSA_GAMMA {

 public:
 NISSA_GAMMA() : G(5) {Init_gamma();}


  void Init_gamma() { for(int i=0;i<5;i++) G[i] = Get_nissa_gamma(i); };
  C_MATRIX PROD(int i, int j) {
    if(i>4 || j>4) crash("PROD called with i,j > 4");
    return complex_prod_matr( G[i], G[j]);
  };
  C_MATRIX G5PROD(int i) {
    if(i>4 ) crash("G5PROD called with i > 4");
    return complex_prod_matr( G[i], G[4]);
  };
  C_MATRIX T(int i) { return TRANSPOSE(G[i]);};
  C_MATRIX DAG(int i) { return DAGGER(G[i]);};
  C_MATRIX ID() { return IDENTITY(4);}
  
  vector<C_MATRIX> G;

  
  

};


void Get_PINGU();
void penguin_spectral_reco();
void TEST_PI_Q2();



#endif
