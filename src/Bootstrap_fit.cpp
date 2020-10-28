#include "../include/Bootstrap_fit.h"


using namespace std;





double Boot_ave(Vfloat& A) {

  if(A.size() ==0) crash(" Boot_average called with an empty vector");
  
  double res=0;;

  for(auto &el: A) res+= el;

  return res/A.size();


}

Pfloat Boot_ave_err(Vfloat& A) {

  if(A.size() ==0) crash("Boot_average called with an empty vector");

  double res=0;
  double err=0;


  for(auto &el: A) { res+=el/A.size(); err+= el*el/A.size();}


  err = (err - res*res)*(A.size())/(A.size() -1);

  return make_pair(res, sqrt(err));

}

double Boot_err(Vfloat& A) {

  if(A.size() ==0) crash("Boot_average called with an empty vector");

  double res=0;
  double err=0;


  for(auto &el: A) { res+=el/A.size(); err+= el*el/A.size();}


  err = (err - res*res)*(A.size())/(A.size() -1);

  return sqrt(err);

}



Pfloat Boot_ave_err(VVfloat& A) {

  if(A.size() ==0) crash("Boot_average called with an empty vector<vector>> ");

  int Nbr = A.size();

  Vfloat Res_boot(Nbr,0) ;
  Vfloat Err_boot(Nbr,0);

   double Tot_result=0;
   double Tot_err=0;

  for(int ibranch=0; ibranch<Nbr;ibranch++) {


    if(A[ibranch].size() == 0) crash("Boot_average called with a vector<vector> A containing an empty <vector>");
    
    for(auto &i_boot: A[ibranch]) {
    
      Res_boot[ibranch] += i_boot/A[ibranch].size();
      Err_boot[ibranch] += pow(i_boot, 2)/A[ibranch].size();
    }

    Err_boot[ibranch] = (Err_boot[ibranch] - pow(Res_boot[ibranch],2))*A[ibranch].size()/(A[ibranch].size() -1);
    Tot_result += Res_boot[ibranch]/Nbr;
    Tot_err += Err_boot[ibranch];
  }

  for(int ibranch=0; ibranch<Nbr;ibranch++) Tot_err += pow(Res_boot[ibranch]- Tot_result,2)/Nbr;

  Tot_err = sqrt(Tot_err);


  return make_pair(Tot_result, Tot_err);
  
}



double Boot_err(VVfloat& A) {

  if(A.size() ==0) crash("Boot_average called with an empty vector<vector>> ");

  int Nbr = A.size();

  Vfloat Res_boot(Nbr,0) ;
  Vfloat Err_boot(Nbr,0);

   double Tot_result=0;
   double Tot_err=0;

  for(int ibranch=0; ibranch<Nbr;ibranch++) {


    if(A[ibranch].size() == 0) crash("Boot_average called with a vector<vector> A containing an empty <vector>");
    
    for(auto &i_boot: A[ibranch]) {
    
      Res_boot[ibranch] += i_boot/A[ibranch].size();
      Err_boot[ibranch] += pow(i_boot, 2)/A[ibranch].size();
    }

    Err_boot[ibranch] = (Err_boot[ibranch] - pow(Res_boot[ibranch],2))*A[ibranch].size()/(A[ibranch].size() -1);
    Tot_result += Res_boot[ibranch]/Nbr;
    Tot_err += Err_boot[ibranch];
  }

  for(int ibranch=0; ibranch<Nbr;ibranch++) Tot_err += pow(Res_boot[ibranch]- Tot_result,2)/Nbr;

  Tot_err = sqrt(Tot_err);


  return Tot_err;
  
}


double Boot_ave(VVfloat& A) {

  if(A.size() ==0) crash("Boot_average called with an empty vector<vector>> ");

  int Nbr = A.size();

  Vfloat Res_boot(Nbr,0) ;

   double Tot_result=0;

  for(int ibranch=0; ibranch<Nbr;ibranch++) {


    if(A[ibranch].size() == 0) crash("Boot_average called with a vector<vector> A containing an empty <vector>");
    
    for(auto &i_boot: A[ibranch]) {
    
      Res_boot[ibranch] += i_boot/A[ibranch].size();
    }

    Tot_result += Res_boot[ibranch]/Nbr;
  }




  return Tot_result;
  
}
