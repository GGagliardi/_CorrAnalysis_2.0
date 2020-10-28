#ifndef __stat__
#define __stat__

#include "numerics.h"
#include "random.h"


using namespace std;



Pfloat BootAve(const Vfloat &BootDistr);
double Compute_boot_cov(const Vfloat& A,const Vfloat& B);
Pfloat JackAve(const Vfloat &JackDistr);
double Compute_jack_cov(const Vfloat&,const Vfloat&);
void Compute_covariance_matrix(bool UseJack, Eigen::MatrixXd& Cov, int nargs,...);
void Compute_correlation_matrix(bool UseJack, Eigen::MatrixXd& Corr, int nargs,...);
  





class distr_t {

 public:
 distr_t()  {}
 distr_t(bool a) : UseJack(a) {}
 distr_t(bool a, Vfloat b) : UseJack(a), distr(b) {}
 distr_t(bool a, int size) : UseJack(a), distr(size,0.0) {}

 //////////////////////////////////////////

 ////////////////////////////////////

 double ave() const;
 double err() const;
 Pfloat ave_err() const;
 int size() const;
 

 bool UseJack;
 Vfloat distr;
};


//operator overloading

distr_t operator*(const distr_t& A, const distr_t& B);
 
distr_t operator+(const distr_t& A,const  distr_t& B);
 
distr_t operator-(const distr_t& A, const distr_t& B);

distr_t operator/(const distr_t& A, const distr_t& B);

double operator%(const distr_t& A,const  distr_t& B) ;

distr_t operator*(const distr_t& A, double B);

distr_t operator*(double B,  const distr_t& A);

distr_t operator+(const distr_t& A, double B);

distr_t operator+(double B,const distr_t& A);
distr_t operator-(const distr_t& A, double B) ;

distr_t operator-(double B,const distr_t& A);
distr_t operator/(const distr_t& A, double B) ;

distr_t operator/(double B,const distr_t& A);

distr_t operator*(const distr_t& A,const Vfloat& B);
 
distr_t operator+(const distr_t& A, const Vfloat& B) ;
 
distr_t operator-(const distr_t& A, const Vfloat& B) ;
distr_t operator/(const distr_t& A,const Vfloat& B) ;

distr_t operator*(const Vfloat &B,const distr_t& A) ;
distr_t operator+(const Vfloat &B, const distr_t& A);
distr_t operator-(const Vfloat& B, const distr_t& A) ;
distr_t operator/(const Vfloat& B, const distr_t& A);

// end operator overloading



class distr_t_list {

 public:
  distr_t_list() {}
  distr_t_list( bool sampling_type, int size) :  UseJack(sampling_type) {
    for(int i_distr=0; i_distr<size;i_distr++) distr_list.emplace_back(sampling_type);
  }
  distr_t_list( bool sampling_type, VVfloat raw_distr_list) : UseJack(sampling_type) {
    for (unsigned int i_distr=0;i_distr<raw_distr_list.size(); i_distr++) this->distr_list.emplace_back(sampling_type, raw_distr_list[i_distr]);
  }
  distr_t_list(bool sampling_type) : UseJack(sampling_type) {}
  distr_t_list(int size, distr_t A) : UseJack(A.UseJack),  distr_list(size,A) {}
  distr_t_list(int sampling_type, int size, int sample_size) : UseJack(sampling_type) {
    for(int i_distr=0; i_distr<size;i_distr++) this->distr_list.emplace_back(sampling_type,sample_size);}
  

 distr_t operator[](int k) {
    if(k>= this->size()) crash("In distr_t_list::operator[](int k) k >= size");
    return this->distr_list[k];
  }

  //////////////////////////////////

 
 Vfloat ave() const ; 
 double ave(int i_distr) const;
 Vfloat err() const;
 double err(int i_distr) const;
 vector<Pfloat> ave_err() const;
 Pfloat ave_err(int i_distr) const ;
 int size() const;
 Vfloat Get_distr_index(int k) const;

  bool UseJack;
  vector<distr_t> distr_list;
};

//operator overloading


  
distr_t_list operator*(const distr_t_list& A, const distr_t_list& B);

distr_t_list operator+(const distr_t_list& A,const distr_t_list& B);

distr_t_list operator-(const distr_t_list& A, const distr_t_list& B);

distr_t_list operator/(const distr_t_list& A,const distr_t_list& B) ;

Vfloat operator%(const distr_t_list& A,const distr_t_list& B);

////////////////////////////////

distr_t_list operator*(const distr_t_list& A, const distr_t& B);
distr_t_list operator*(const distr_t& B, const distr_t_list& A);

distr_t_list operator+(const distr_t_list& A,const distr_t& B);
distr_t_list operator+(const distr_t& B, const distr_t_list& A);

distr_t_list operator-(const distr_t_list& A,const distr_t& B);

distr_t_list operator-(const distr_t& B,const distr_t_list& A);

distr_t_list operator/(const distr_t_list& A,const distr_t& B);

distr_t_list operator/(const distr_t& B,const distr_t_list& A);
  

Vfloat operator%(const distr_t_list& A,const distr_t& B);

Vfloat operator%(const distr_t& B, const distr_t_list& A);
///////////////////////////////


distr_t_list operator*(const distr_t_list& A, double B);
distr_t_list operator*(double B, const distr_t_list& A);

distr_t_list operator+(const distr_t_list& A, double B);

distr_t_list operator+(double B, const distr_t_list& A);

distr_t_list operator-(const distr_t_list& A, double B);

distr_t_list operator-(double B, const distr_t_list& A);

distr_t_list operator/(const distr_t_list& A, double B);
  
distr_t_list operator/(double B, const distr_t_list& A);

distr_t_list operator+(const distr_t_list &A, const Vfloat& B);

distr_t_list operator+(const Vfloat& B, const distr_t_list& A) ;


distr_t_list operator*(const distr_t_list &A, const Vfloat& B);

distr_t_list operator*(const Vfloat& B, const distr_t_list& A);


//end operator overloading

class Jackknife {

 public:
  Jackknife() {}
  Jackknife(int block_size_max, int Njacks) {
    ThermalMeasToErase=0;
    this->block_size_max = block_size_max;
    this->Njacks= Njacks;
    Verbose_jack=0;
    ReturnBlockSizeMax=0;
    }
  distr_t DoJack(int , ...);

  


  int ThermalMeasToErase;
  int block_size_max;
  int Verbose_jack;
  int Njacks;
  bool ReturnBlockSizeMax;


} ;




class Bootstrap {

 public:
  Bootstrap() {}
  Bootstrap(int Nboots, int seed, int size) : RN(seed, size-1) { 
    this->Nboots= Nboots;
    this->size=size;
    this->generate_data();
    
  }
  distr_t DoBoot(int Nobs, ...);
 
  
  int Nboots;
  int size;
  vector<vector<int>> bootstrap_data;

  private:
  RandomMersenne RN;
  void generate_data() {
    bootstrap_data.resize(Nboots);
    for(int iboot=0; iboot<Nboots;iboot++)
      for(int el=0; el < size; el ++) bootstrap_data[iboot].push_back( RN());
  }

} ;







#endif
