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
 * construct constructs a standard MPO object, by creating a double layer peps object from the top row
 */
template<typename T>
MPO<T>::MPO(int row,const PEPS<T> &peps_1,const PEPS<T> &peps_2) : vector< TArray<T,4> >(PEPS<T>::lat.gLx()) {

   D = peps_1.gD() * peps_2.gD();

   //c == 0
   (*this)[0].resize(1,D,D,D);

   for(int c = 1;c < this->size() - 1;++c)
      (*this)[c].resize(D,D,D,D);

   (*this)[this->size() - 1].resize(D,D,D,1);

   enum {i,j,k,l,m,n,o,p,q};

   TArray<T,8> tmp;

   for(int c = 0;c < this->size();++c){

      tmp.clear();

      Contract((T)1.0,peps_1(row,c),shape(i,j,k,l,m),peps_2(row,c),shape(n,o,k,p,q),(T)0.0,tmp,shape(i,n,j,o,l,p,m,q));

      CopyR(tmp,(*this)[c]);

   }

}

/**
 * copy constructor
 */
template<typename T>
MPO<T>::MPO(const MPO<T> &mps_copy) : vector< TArray<T,4> >(mps_copy) {

   D = mps_copy.gD();

}

/**
 * empty destructor
 */
template<typename T>
MPO<T>::~MPO(){ }

/**
 * @return virtual dimension of the MPO
 */
template<typename T>
int MPO<T>::gD() const {

   return D;

}

template MPO<double>::MPO(int,const PEPS<double> &,const PEPS<double> &);
template MPO< complex<double> >::MPO(int,const PEPS< complex<double> > &,const PEPS< complex<double> > &);

template MPO<double>::MPO(const MPO<double> &);
template MPO< complex<double> >::MPO(const MPO< complex<double> > &);

template MPO<double>::~MPO();
template MPO< complex<double> >::~MPO();

template int MPO<double>::gD() const;
template int MPO< complex<double> >::gD() const;
