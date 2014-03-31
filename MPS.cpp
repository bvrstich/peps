#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;

#include "include.h"

/**
 * construct constructs a standard MPS object, by creating a double layer peps object from the top row
 */
template<typename T>
MPS<T>::MPS(const PEPS<T> &peps_1,const PEPS<T> &peps_2) : vector< TArray<T,3> >(PEPS<T>::lat.gLx()) {

   D = peps_1.gD() * peps_2.gD();

   //c == 0
   (*this)[0].resize(1,D,D);

   for(int c = 1;c < this->size() - 1;++c)
      (*this)[c].resize(D,D,D);

   (*this)[this->size() - 1].resize(D,D,1);

   enum {i,j,k,l,m,n,o,p,q};

   TArray<T,8> tmp;

   for(int c = 0;c < this->size();++c){

      tmp.clear();

      Contract((T)1.0,peps_1(0,c),shape(i,j,k,l,m),peps_2(0,c),shape(n,o,k,p,q),(T)0.0,tmp,shape(i,n,j,o,l,p,m,q));

      CopyR(tmp,(*this)[c]);

   }

}

/**
 * copy constructor
 */
template<typename T>
MPS<T>::MPS(const MPS<T> &mps_copy) : vector< TArray<T,3> >(mps_copy) {

   D = mps_copy.gD();

}

/**
 * empty destructor
 */
template<typename T>
MPS<T>::~MPS(){ }

/**
 * @return virtual dimension of the MPS
 */
template<typename T>
int MPS<T>::gD() const {

   return D;

}

/**
 * act with an MPO on this MPS, resulting MPS is returned as *this object
 * @param uplo if == 'U' contract with the upper physical index of the MPO, if == 'L', contract with the lower
 * @param mpo the MPO
 */
template<typename T>
void MPS<T>::gemv(char uplo,const MPO<T> &mpo){

   int L = this->size();

   if(uplo == 'U'){

      //first site
      TArray<T,5> tmp;

      enum {i,j,k,l,m,n,o,p,q};

      Contract((T)1.0,mpo[0],shape(i,j,k,l),(*this)[0],shape(m,j,n),(T)0.0,tmp,shape(m,i,k,n,l));

      (*this)[0].resize(1,D,D*D);

      CopyR(tmp,(*this)[0]);

      //middle sites
      for(int c = 1;c < L - 1;++c){

         tmp.clear();

         Contract((T)1.0,mpo[c],shape(i,j,k,l),(*this)[c],shape(m,j,n),(T)0.0,tmp,shape(m,i,k,n,l));

         (*this)[c].resize(D*D,D,D*D);

         CopyR(tmp,(*this)[c]);

      }

      //last site
      tmp.clear();

      Contract((T)1.0,mpo[L - 1],shape(i,j,k,l),(*this)[L - 1],shape(m,j,n),(T)0.0,tmp,shape(m,i,k,n,l));

      (*this)[L - 1].resize(D*D,D,1);

      CopyR(tmp,(*this)[L - 1]);

   }
   else{//L

      //first site
      TArray<T,5> tmp;

      enum {i,j,k,l,m,n,o,p,q};

      Contract((T)1.0,mpo[0],shape(i,j,k,l),(*this)[0],shape(m,k,n),(T)0.0,tmp,shape(m,i,j,n,l));

      (*this)[0].resize(1,D,D*D);

      CopyR(tmp,(*this)[0]);

      //middle sites
      for(int c = 1;c < L - 1;++c){

         tmp.clear();

         Contract((T)1.0,mpo[c],shape(i,j,k,l),(*this)[c],shape(m,k,n),(T)0.0,tmp,shape(m,i,j,n,l));

         (*this)[c].resize(D*D,D,D*D);

         CopyR(tmp,(*this)[c]);

      }

      //last site
      tmp.clear();

      Contract((T)1.0,mpo[L - 1],shape(i,j,k,l),(*this)[L - 1],shape(m,k,n),(T)0.0,tmp,shape(m,i,j,n,l));

      (*this)[L - 1].resize(D*D,D,1);

      CopyR(tmp,(*this)[L - 1]);

   }

   D *= mpo.gD();

}

template MPS<double>::MPS(const PEPS<double> &,const PEPS<double> &);
template MPS< complex<double> >::MPS(const PEPS< complex<double> > &,const PEPS< complex<double> > &);

template MPS<double>::MPS(const MPS<double> &);
template MPS< complex<double> >::MPS(const MPS< complex<double> > &);

template MPS<double>::~MPS();
template MPS< complex<double> >::~MPS();

template int MPS<double>::gD() const;
template int MPS< complex<double> >::gD() const;

template void MPS<double>::gemv(char uplo,const MPO<double> &mpo);
template void MPS< complex<double> >::gemv(char uplo,const MPO< complex<double> > &mpo);
