#ifndef __complex_stat__
#define __complex_stat__

#include "numerics.h"
#include "random.h"
#include "stat.h"

using namespace std;




  





class complex_distr_t {

 public:
  complex_distr_t() : RE(1), IM(1)  {UseJack = 1;}
  complex_distr_t(bool a) : UseJack(a) , RE(a), IM(a) {}
  complex_distr_t(bool a, Vfloat b) : UseJack(a), RE(a,b), IM(a,b) {}
  complex_distr_t(bool a, int size) : UseJack(a), RE(a,size), IM(a,size) {}
  complex_distr_t(bool a, Vfloat r, Vfloat i) : UseJack(a), RE(a,r) , IM(a,i) {}
  complex_distr_t( distr_t r, distr_t i) : UseJack(r.UseJack) { assert(r.UseJack == i.UseJack); assert(r.size() == i.size()); this->RE=r; this->IM=i;}
  
  //////////////////////////////////////////
  
  ////////////////////////////////////

  double ave(int i) const;
  double err(int i) const;
  complex<double> ave() const;
  complex<double> err() const;
  int size() const;
  complex_distr_t dagger() const { complex_distr_t ret(this->UseJack); ret.RE= this->RE; ret.IM = -1.0*this->IM; return ret; }
  
  bool UseJack;
  distr_t RE;
  distr_t IM;
};


//operator overloading

complex_distr_t operator*(const complex_distr_t& A, const complex_distr_t& B);
 
complex_distr_t operator+(const complex_distr_t& A,const  complex_distr_t& B);

complex_distr_t operator-(const complex_distr_t &A, const complex_distr_t &B);


complex_distr_t operator/(const complex_distr_t &A, const complex_distr_t &B);



complex<double> operator%(const complex_distr_t &A, const complex_distr_t &B);
complex<double> operator%(const complex_distr_t &A, const distr_t &B);
complex<double> operator%(const distr_t& A,const  complex_distr_t& B) ;


complex_distr_t operator*(const complex_distr_t &A, double B);
complex_distr_t operator*(const complex_distr_t& A, complex<double> B);

complex_distr_t operator*(double B, const complex_distr_t &A);
complex_distr_t operator*(complex<double> B,  const complex_distr_t& A);

complex_distr_t operator+(const complex_distr_t &A, double B);
complex_distr_t operator+(const complex_distr_t& A, complex<double> B);

complex_distr_t operator+(double B, const complex_distr_t &A);
complex_distr_t operator+(complex<double> B,const complex_distr_t& A);
complex_distr_t operator-(const complex_distr_t &A, double B);
complex_distr_t operator-(const complex_distr_t& A, complex<double> B) ;

complex_distr_t operator-(double B, const complex_distr_t &A);
complex_distr_t operator-(complex<double> B,const complex_distr_t& A);
complex_distr_t operator/(const complex_distr_t &A, double B);
complex_distr_t operator/(const complex_distr_t& A, complex<double> B) ;

complex_distr_t operator/(double B, const complex_distr_t &A);
complex_distr_t operator/(complex<double> B,const complex_distr_t& A);

complex_distr_t operator*(const complex_distr_t& A,const distr_t& B);
 
complex_distr_t operator+(const complex_distr_t& A, const distr_t& B) ;
 
complex_distr_t operator-(const complex_distr_t& A, const distr_t& B) ;
complex_distr_t operator/(const complex_distr_t& A,const distr_t& B) ;

complex_distr_t operator*(const distr_t &B,const complex_distr_t& A) ;
complex_distr_t operator+(const distr_t &B, const complex_distr_t& A);
complex_distr_t operator-(const distr_t &B, const complex_distr_t &A);
complex_distr_t operator/(const distr_t& B, const complex_distr_t& A);

// end operator overloading



class complex_distr_t_list {

public:
  complex_distr_t_list() {UseJack=1;}
  complex_distr_t_list( bool sampling_type, int size) :  UseJack(sampling_type) {
    for(int i_distr=0; i_distr<size;i_distr++) distr_list.emplace_back(sampling_type);
  }
  complex_distr_t_list(bool sampling_type) : UseJack(sampling_type) {}
  complex_distr_t_list(int size, complex_distr_t A) : UseJack(A.UseJack),  distr_list(size,A) {}
  complex_distr_t_list(int sampling_type, int size, int sample_size) : UseJack(sampling_type) {
    for(int i_distr=0; i_distr<size;i_distr++) this->distr_list.emplace_back(sampling_type,sample_size);}
  complex_distr_t_list( distr_t_list RE, distr_t_list IM) : UseJack(RE.UseJack) {
    assert(RE.UseJack == IM.UseJack);
    assert(RE.size() == IM.size());
    for(int t=0;t<(signed)RE.size();t++)  this->distr_list.emplace_back( RE[t], IM[t]);
  }
  complex_distr_t_list ( bool a, VVfloat r, VVfloat s) : UseJack(a) {
    assert(r.size() == s.size());
    for(int t=0;t<(signed)r.size();t++) {
      assert(r[t].size() == s[t].size());
      this->distr_list.emplace_back( a, r[t], s[t]);
    }
  }
  
  
  complex_distr_t operator[](int k) {
    if(k>= this->size()) crash("In distr_t_list::operator[](int k) k >= size");
    return this->distr_list[k];
  }

  VVfloat Get_vvector(int k) {

    VVfloat ret;
    for(int t=0;t<this->size();t++) {
      if(k==0)  ret.push_back( this->distr_list[t].RE.distr );
      else ret.push_back( this->distr_list[t].IM.distr );
							      
    }

    return ret;
  }

  void Add(const complex_distr_t_list &A, int c) {

    assert( this->size() == A.size());
    for(int t=0; t<this->size();t++) {
      assert( c < this->distr_list[t].size());
      assert( A.distr_list[t].size() == 1);
      this->distr_list[t].RE.distr[c] += A.distr_list[t].RE.distr[0];
      this->distr_list[t].IM.distr[c] += A.distr_list[t].IM.distr[0];
    }
    
    return;
  }

  static complex_distr_t_list external_prod(const complex_distr_t_list& A, const complex_distr_t_list& B) {

    assert(A.UseJack == B.UseJack);
    complex_distr_t_list ret(A.UseJack);
    int TT= A.size()*B.size();

    for(int t=0;t<TT;t++) {

      int t1= t%A.size();
      int t2= (t/A.size())%B.size();

      ret.distr_list.push_back( A.distr_list[t1]*B.distr_list[t2]);

    }

    return ret;
  }
  //////////////////////////////////
  
  
  Vfloat RE_ave() const ;
  Vfloat IM_ave() const;
  double RE_ave(int i_distr) const;
  double IM_ave(int i_distr) const;
  Vfloat RE_err() const;
  Vfloat IM_err() const;
  double RE_err(int i_distr) const;
  double IM_err(int i_distr) const;
  int size() const;
  distr_t_list RE() const {
    distr_t_list ret(this->UseJack);
    for(int t=0;t<(signed)this->distr_list.size();t++) ret.distr_list.push_back( this->distr_list[t].RE );
    return ret;
  }
  distr_t_list IM() const {
    distr_t_list ret(this->UseJack);
    for(int t=0;t<(signed)this->distr_list.size();t++) ret.distr_list.push_back( this->distr_list[t].IM );
    return ret;
  }
  complex_distr_t_list dagger() const {
    complex_distr_t_list ret(this->UseJack);
    for(int t=0;t<this->size();t++) ret.distr_list.push_back( this->distr_list[t].dagger());
    return ret;
  }
 

  

  bool UseJack;
  vector<complex_distr_t> distr_list;
};

//operator overloading




  
complex_distr_t_list operator*(const complex_distr_t_list& A, const complex_distr_t_list& B);

complex_distr_t_list operator+(const complex_distr_t_list& A,const complex_distr_t_list& B);

complex_distr_t_list operator-(const complex_distr_t_list& A, const complex_distr_t_list& B);

complex_distr_t_list operator/(const complex_distr_t_list &A,
                               const complex_distr_t_list &B);

complex_distr_t_list operator/(const complex_distr_t_list& A, const distr_t_list& B);

vector<complex<double>> operator%(const complex_distr_t_list& A,const complex_distr_t_list& B);

////////////////////////////////

complex_distr_t_list operator*(const complex_distr_t_list& A, const distr_t& B);
complex_distr_t_list operator*(const distr_t& B, const complex_distr_t_list& A);

complex_distr_t_list operator+(const complex_distr_t_list& A,const distr_t& B);
complex_distr_t_list operator+(const distr_t& B, const complex_distr_t_list& A);

complex_distr_t_list operator-(const complex_distr_t_list& A,const distr_t& B);

complex_distr_t_list operator-(const distr_t& B,const complex_distr_t_list& A);

complex_distr_t_list operator/(const complex_distr_t_list& A,const distr_t& B);

complex_distr_t_list operator/(const distr_t& B,const complex_distr_t_list& A);
  

vector<complex<double>> operator%(const complex_distr_t_list& A,const distr_t& B);

vector<complex<double>> operator%(const distr_t& B, const complex_distr_t_list& A);
///////////////////////////////


complex_distr_t_list operator*(const complex_distr_t_list& A, const complex_distr_t& B);
complex_distr_t_list operator*(const complex_distr_t& B, const complex_distr_t_list& A);

complex_distr_t_list operator+(const complex_distr_t_list& A,const complex_distr_t& B);
complex_distr_t_list operator+(const complex_distr_t& B, const complex_distr_t_list& A);

complex_distr_t_list operator-(const complex_distr_t_list& A,const complex_distr_t& B);

complex_distr_t_list operator-(const complex_distr_t& B,const complex_distr_t_list& A);

complex_distr_t_list operator/(const complex_distr_t_list& A,const complex_distr_t& B);

complex_distr_t_list operator/(const complex_distr_t& B,const complex_distr_t_list& A);
  

vector<complex<double>> operator%(const complex_distr_t_list& A,const complex_distr_t& B);

vector<complex<double>> operator%(const complex_distr_t& B, const complex_distr_t_list& A);


////////////////////////////////

complex_distr_t_list operator*(const complex_distr_t_list& A, double B);
complex_distr_t_list operator*(double B, const complex_distr_t_list& A);

complex_distr_t_list operator+(const complex_distr_t_list& A, double B);

complex_distr_t_list operator+(double B, const complex_distr_t_list& A);

complex_distr_t_list operator-(const complex_distr_t_list& A, double B);

complex_distr_t_list operator-(double B, const complex_distr_t_list& A);

complex_distr_t_list operator/(const complex_distr_t_list& A, double B);

complex_distr_t_list operator/(double B, const complex_distr_t_list &A);

//////////////////////////////


complex_distr_t_list operator*(const complex_distr_t_list& A, complex<double> B);
complex_distr_t_list operator*(complex<double> B, const complex_distr_t_list& A);

complex_distr_t_list operator+(const complex_distr_t_list& A, complex<double> B);

complex_distr_t_list operator+(complex<double> B, const complex_distr_t_list& A);

complex_distr_t_list operator-(const complex_distr_t_list& A, complex<double> B);

complex_distr_t_list operator-(complex<double> B, const complex_distr_t_list& A);

complex_distr_t_list operator/(const complex_distr_t_list& A, complex<double> B);

complex_distr_t_list operator/(complex<double> B, const complex_distr_t_list &A);


/////////////////////////////




complex_distr_t_list operator+(const complex_distr_t_list &A, const Vfloat &B);
complex_distr_t_list operator-(const complex_distr_t_list &A, const Vfloat& B);

complex_distr_t_list operator+(const Vfloat &B, const complex_distr_t_list &A);
complex_distr_t_list operator-(const Vfloat& B, const complex_distr_t_list& A);


complex_distr_t_list operator*(const complex_distr_t_list &A, const Vfloat& B);

complex_distr_t_list operator*(const Vfloat& B, const complex_distr_t_list& A);


//end operator overloading


#endif
