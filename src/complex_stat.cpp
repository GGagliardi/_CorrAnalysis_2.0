#include "../include/complex_stat.h"

using namespace std;






double complex_distr_t::ave(int i) const {

  if(i==0) return this->RE.ave();

  return this->IM.ave();
}

double complex_distr_t::err(int i) const {

  if(i==0) return this->RE.err();

  return this->IM.err();
}

complex<double> complex_distr_t::ave() const {
  return this->RE.ave() + 1i*this->IM.ave();
}

complex<double> complex_distr_t::err() const {
  return this->RE.err() + 1i*this->IM.err();
}

 
int complex_distr_t::size() const { assert(RE.size() == IM.size());  return this->RE.size();}


Vfloat complex_distr_t_list::RE_ave() const {
    Vfloat res;
    for(int i=0; i < this->size();i++) res.push_back(distr_list[i].ave(0));
    return res;
}

Vfloat complex_distr_t_list::IM_ave() const {
    Vfloat res;
    for(int i=0; i < this->size();i++) res.push_back(distr_list[i].ave(1));
    return res;
  }

double complex_distr_t_list::RE_ave(int i_distr) const {
    if(i_distr >= (signed)distr_list.size()) crash("In distr_t_list function ave called with positional argument greater than distr_list size");
    return distr_list[i_distr].ave(0);
}

double complex_distr_t_list::IM_ave(int i_distr) const {
    if(i_distr >= (signed)distr_list.size()) crash("In distr_t_list function ave called with positional argument greater than distr_list size");
    return distr_list[i_distr].ave(1);
}

Vfloat complex_distr_t_list::RE_err()  const{
    Vfloat res;
    for(int i=0; i < this->size();i++) res.push_back(this->distr_list[i].err(0));
    return res;
}

Vfloat complex_distr_t_list::IM_err()  const{
    Vfloat res;
    for(int i=0; i < this->size();i++) res.push_back(this->distr_list[i].err(1));
    return res;
  }

double complex_distr_t_list::RE_err(int i_distr) const {
    if(i_distr >=  (signed)distr_list.size()) crash("In distr_t_list function ave called with positional argument greater than distr_list size");
    return this->distr_list[i_distr].err(0);
}

double complex_distr_t_list::IM_err(int i_distr) const {
    if(i_distr >=  (signed)distr_list.size()) crash("In distr_t_list function ave called with positional argument greater than distr_list size");
    return this->distr_list[i_distr].err(1);
  }


int complex_distr_t_list::size() const { return distr_list.size();}



//operator overloading

complex_distr_t operator*(const complex_distr_t& A, const complex_distr_t& B) {
  complex_distr_t ret = A;
  ret.RE = A.RE*B.RE - A.IM*B.IM;
  ret.IM =  A.RE*B.IM + A.IM*B.RE; 
  return ret;
 }
 
complex_distr_t operator+(const complex_distr_t& A,const  complex_distr_t& B) {
  complex_distr_t ret = A;
  ret.RE = A.RE + B.RE;
  ret.IM = A.IM + B.IM;
  return ret;
}

 
complex_distr_t operator-(const complex_distr_t& A, const complex_distr_t& B) {
  complex_distr_t ret = A;
  ret.RE = A.RE - B.RE;
  ret.IM = A.IM - B.IM;
  return ret;
}

complex_distr_t operator/(const complex_distr_t& A, const distr_t& B) {
  complex_distr_t ret=A;

  ret.RE = A.RE/B;
  ret.IM = A.IM/B;

  return ret;
  
}

complex_distr_t operator/(const complex_distr_t& A, const complex_distr_t& B) {
  complex_distr_t ret;
  return  A*B.dagger()/( B.RE*B.RE + B.IM*B.IM) ;
}

complex<double> operator%(const complex_distr_t& A,const  complex_distr_t& B) {

  return (A.RE%B.RE) + 1i*(A.IM%B.IM);
  
}

complex<double> operator%(const complex_distr_t& A,const  distr_t& B) {

  return (A.RE%B) + 1i*(A.IM%B);
}

complex<double> operator%(const distr_t& A,const  complex_distr_t& B) {

  return B%A;
}
  
  

complex_distr_t operator*(const complex_distr_t& A, double B) {
  complex_distr_t ret=A;
  
  ret.RE= B*ret.RE;
  ret.IM= B*ret.IM;

  return ret;
}

complex_distr_t operator*(const complex_distr_t& A, complex<double> B) {
  complex_distr_t ret=A;

 
  ret.RE= B.real()*A.RE -1.0*B.imag()*A.IM;
  ret.IM= A.RE*B.imag() + A.IM*B.real();
 
  return ret;
  
 }


complex_distr_t operator*(double B, const complex_distr_t &A) { return A * B; }

complex_distr_t operator*(complex<double> B,  const complex_distr_t& A) {
   return A*B;
}

complex_distr_t operator+(const complex_distr_t& A, double B) {
   complex_distr_t ret=A;

   ret.RE = A.RE + B;
   ret.IM = A.IM + B;

   return ret;  
 }

complex_distr_t operator+(double B, const complex_distr_t &A) { return A + B; }

complex_distr_t operator+(const complex_distr_t& A, complex<double> B) {
  complex_distr_t ret=A;

  ret.RE = A.RE + B.real();
  ret.IM = A.IM + B.imag();

  return ret;  
 }

 complex_distr_t operator+(complex<double> B,const complex_distr_t& A) {
   return A+B;
 }

 complex_distr_t operator-(const complex_distr_t& A, double B) {

   complex_distr_t ret=A;

   ret.RE = A.RE - B;
   ret.IM = A.IM - B;
   
   return ret;
 }

 complex_distr_t operator-(const complex_distr_t& A, complex<double> B) {

   complex_distr_t ret=A;

   ret.RE = A.RE - B.real();
   ret.IM = A.IM - B.imag();
   
   return ret;  
 }


 

 complex_distr_t operator-(double B,const complex_distr_t& A) {

   return (-1.0)*A + B;
 }

complex_distr_t operator-(complex<double> B,const complex_distr_t& A) {

   return (-1.0)*A + B;
 }

 complex_distr_t operator/(const complex_distr_t& A, double B) {

   complex_distr_t ret=A;

   ret.RE = A.RE/B;
   ret.IM = A.IM/B;
   
   return ret;
 }

complex_distr_t operator/(const complex_distr_t& A, complex<double> B) {
    return A*(1.0/(B));
  }
 

 complex_distr_t operator/(double B,const complex_distr_t& A) {
   return B*A.dagger()/(A.RE*A.RE + A.IM*A.IM);
 }

complex_distr_t operator/(complex<double> B,const complex_distr_t& A) {
  return B*( 1.0/A );
 }

 complex_distr_t operator*(const complex_distr_t& A,const distr_t& B) {
   complex_distr_t ret=A;
   ret.RE = A.RE*B;
   ret.IM = A.IM*B;
   return ret;
 }
 
 complex_distr_t operator+(const complex_distr_t& A, const distr_t& B) {
   complex_distr_t ret=A;
   ret.RE = A.RE +B;
   ret.IM = A.IM +B;
   return ret;
 }
 
 complex_distr_t operator-(const complex_distr_t& A, const distr_t& B) {
   complex_distr_t ret=A;
   ret.RE = A.RE -B;
   ret.IM = A.IM -B;
   return ret;
 }
 
 
complex_distr_t operator*(const distr_t &B,const complex_distr_t& A) {
  return A*B;
}
complex_distr_t operator+(const distr_t &B, const complex_distr_t& A) {
  return A+B;
}
complex_distr_t operator-(const distr_t& B, const complex_distr_t& A) {
  return -1.0*A + B;
}


complex_distr_t operator/(const distr_t& B, const complex_distr_t& A) {

  return B*A.dagger()/(A.RE*A.RE + A.IM*A.IM);
}





  
complex_distr_t_list operator*(const complex_distr_t_list& A, const complex_distr_t_list& B) {
  if(A.size() != B.size()) crash("In complex_distr_t_list, call to A*B is invalid. A and B have different sizes");
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] =A.distr_list[i]*B.distr_list[i];
    return res;
}

complex_distr_t_list operator+(const complex_distr_t_list& A,const complex_distr_t_list& B) {

  if(A.size() != B.size()) crash("In complex_distr_t_list, call to A+B is invalid. A and B have different sizes");
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]+B.distr_list[i];
  return res;
}

complex_distr_t_list operator-(const complex_distr_t_list& A, const complex_distr_t_list& B) {

  if(A.size() != B.size()) crash("In complex_distr_t_list, call to A-B is invalid. A and B have different sizes");
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]-B.distr_list[i];
  return res;
}

complex_distr_t_list operator/(const complex_distr_t_list& A,const complex_distr_t_list& B) {
  if(A.size() != B.size()) crash("In complex_distr_t_list, call to A/B is invalid. A and B have different sizes");
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]/B.distr_list[i];
    
  return res;
}

complex_distr_t_list operator/(const complex_distr_t_list& A,const distr_t_list& B) {
  if(A.size() != B.size()) crash("In complex_distr_t_list, call to A/B is invalid. A and B have different sizes");
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]/B.distr_list[i];
    
  return res;
}

vector<complex<double>> operator%(const complex_distr_t_list& A,const complex_distr_t_list& B) {
  if(A.size() != B.size()) crash("In complex_distr_t_list, call to A%B is invalid. A and B have different sizes");
  vector<complex<double>> cov(A.size());
  for(int i=0;i<A.size();i++) {
    if(A.distr_list[i].size() != B.distr_list[i].size()) crash("In complex_distr_t_list, call to A%B invalid, the two operands have different distr_t sizes");
    cov[i] = A.distr_list[i]%B.distr_list[i];
  }
  return cov;
}

////////////////////////////////

complex_distr_t_list operator*(const complex_distr_t_list& A, const complex_distr_t& B) {

  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] =A.distr_list[i]*B;
  return res;
}
complex_distr_t_list operator*(const complex_distr_t& B, const complex_distr_t_list& A) {return A*B;}

complex_distr_t_list operator+(const complex_distr_t_list& A,const complex_distr_t& B) {
  complex_distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]+B;
  return res;
}
complex_distr_t_list operator+(const complex_distr_t& B, const complex_distr_t_list& A) { return A+B;}

complex_distr_t_list operator-(const complex_distr_t_list& A,const complex_distr_t& B) {
  complex_distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]-B;
  return res;
}

complex_distr_t_list operator-(const complex_distr_t& B,const complex_distr_t_list& A) {
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = B- A.distr_list[i];
  return res;
}

complex_distr_t_list operator/(const complex_distr_t_list& A,const complex_distr_t& B) {
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]/B;
  return res;
}

complex_distr_t_list operator/(const complex_distr_t& B,const complex_distr_t_list& A) {
  complex_distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = B/A.distr_list[i];
  return res;
}
  

vector<complex<double>> operator%(const complex_distr_t_list& A,const complex_distr_t& B) {
  vector<complex<double>> cov(A.size());
  for( int i=0;i<A.size();i++) {
    if(A.distr_list[i].size() != B.size()) crash("In complex_distr_t_list, call to A%B invalid, the two operands have different complex_distr_t sizes");
    cov[i] = A.distr_list[i]%B;
  }
  return cov;
}

vector<complex<double>> operator%(const complex_distr_t& B, const complex_distr_t_list& A) {
  return A%B;
}


/////////////////////////////////////////////////////////




complex_distr_t_list operator*(const complex_distr_t_list& A, const distr_t& B) {

  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] =A.distr_list[i]*B;
  return res;
}
complex_distr_t_list operator*(const distr_t& B, const complex_distr_t_list& A) {return A*B;}

complex_distr_t_list operator+(const complex_distr_t_list& A,const distr_t& B) {
  complex_distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]+B;
  return res;
}
complex_distr_t_list operator+(const distr_t& B, const complex_distr_t_list& A) { return A+B;}

complex_distr_t_list operator-(const complex_distr_t_list& A,const distr_t& B) {
  complex_distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]-B;
  return res;
}

complex_distr_t_list operator-(const distr_t& B,const complex_distr_t_list& A) {
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = B- A.distr_list[i];
  return res;
}

complex_distr_t_list operator/(const complex_distr_t_list& A,const distr_t& B) {
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]/B;
  return res;
}

complex_distr_t_list operator/(const distr_t& B,const complex_distr_t_list& A) {
  complex_distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = B/A.distr_list[i];
  return res;
}
  

vector<complex<double>> operator%(const complex_distr_t_list& A,const distr_t& B) {
  vector<complex<double>> cov(A.size());
  for( int i=0;i<A.size();i++) {
    if(A.distr_list[i].size() != B.size()) crash("In complex_distr_t_list, call to A%B invalid, the two operands have different complex_distr_t sizes");
    cov[i] = A.distr_list[i]%B;
  }
  return cov;
}

vector<complex<double>> operator%(const distr_t& B, const complex_distr_t_list& A) {
  return A%B;
}



///////////////////////////////


complex_distr_t_list operator*(const complex_distr_t_list& A, double B) {
  complex_distr_t_list res(A.UseJack,A.size());
    for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]*B;
    return res;
}

complex_distr_t_list operator*(double B, const complex_distr_t_list& A) {return A*B;}

complex_distr_t_list operator+(const complex_distr_t_list& A, double B) {
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]+B;
     return res;
}

complex_distr_t_list operator+(double B, const complex_distr_t_list& A) { return A+B;}

complex_distr_t_list operator-(const complex_distr_t_list& A, double B) {
  complex_distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]-B;
  return res;
  }

complex_distr_t_list operator-(double B, const complex_distr_t_list& A) {
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = B- A.distr_list[i];
  return res;
}

complex_distr_t_list operator/(const complex_distr_t_list& A, double B) {
  complex_distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]/B;
  return res;
}
  
complex_distr_t_list operator/(double B, const complex_distr_t_list& A) {
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = B/A.distr_list[i];
  return res;
}

//############################################################################################


complex_distr_t_list operator*(const complex_distr_t_list& A, complex<double> B) {
  complex_distr_t_list res(A.UseJack,A.size());
    for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]*B;
    return res;
}

complex_distr_t_list operator*(complex<double> B, const complex_distr_t_list& A) {return A*B;}

complex_distr_t_list operator+(const complex_distr_t_list& A, complex<double> B) {
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]+B;
     return res;
}

complex_distr_t_list operator+(complex<double> B, const complex_distr_t_list& A) { return A+B;}

complex_distr_t_list operator-(const complex_distr_t_list& A, complex<double> B) {
  complex_distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]-B;
  return res;
  }

complex_distr_t_list operator-(complex<double> B, const complex_distr_t_list& A) {
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = B- A.distr_list[i];
  return res;
}

complex_distr_t_list operator/(const complex_distr_t_list& A, complex<double> B) {
  complex_distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]/B;
  return res;
}
  
complex_distr_t_list operator/(complex<double> B, const complex_distr_t_list& A) {
  complex_distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = B/A.distr_list[i];
  return res;
}



//#############################################################################################

complex_distr_t_list operator+(const complex_distr_t_list &A, const Vfloat& B) {
  complex_distr_t_list res(A.UseJack, A.size());
  if(A.size() != (signed)B.size()) crash("Call to operator complex_distr_t_list*Vfloat is invalid, sizeof(Vfloat) and size(complex_distr_t_list) do not coincide");

  for(int i=0; i<A.size();i++) res.distr_list[i] = A.distr_list[i]+B[i];

  return res;
}

complex_distr_t_list operator-(const complex_distr_t_list &A, const Vfloat& B) {
  complex_distr_t_list res(A.UseJack, A.size());
  if(A.size() != (signed)B.size()) crash("Call to operator complex_distr_t_list*Vfloat is invalid, sizeof(Vfloat) and size(complex_distr_t_list) do not coincide");

  for(int i=0; i<A.size();i++) res.distr_list[i] = A.distr_list[i]-B[i];

  return res;
}

complex_distr_t_list operator+(const Vfloat &B, const complex_distr_t_list &A) { return A + B; }
complex_distr_t_list operator-(const Vfloat& B, const complex_distr_t_list& A) { return -1.0*(A-B); }


complex_distr_t_list operator*(const complex_distr_t_list &A, const Vfloat& B) {
  complex_distr_t_list res(A.UseJack, A.size());
  if(A.size() != (signed)B.size()) crash("Call to operator complex_distr_t_list*Vfloat is invalid, sizeof(Vfloat) and size(complex_distr_t_list) do not coincide");

  for(int i=0; i<A.size();i++) res.distr_list[i] = A.distr_list[i]*B[i];

  return res;
}

complex_distr_t_list operator*(const Vfloat& B, const complex_distr_t_list& A) { return A*B;}





