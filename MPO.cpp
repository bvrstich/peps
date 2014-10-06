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

using namespace global;

/**
 * construct constructs a standard MPO object, by creating a double layer peps object
 * @param option 'V' == vertical stripe, 'H' == horizontal stripe
 * @param rc row or column index
 */
template<typename T>
MPO<T>::MPO(char option,int rc,const PEPS<T> &peps_1,const PEPS<T> &peps_2) : vector< TArray<T,4> >() {

   D = peps_1.gD() * peps_2.gD();

   if(option == 'H'){//horizontal

      this->resize(Lx);

      enum {i,j,k,l,m,n,o,p,q};

      TArray<T,8> tmp;

      for(int c = 0;c < Lx;++c){

         int DL = peps_1(rc,c).shape(0) * peps_2(rc,c).shape(0);
         int DU = peps_1(rc,c).shape(1) * peps_2(rc,c).shape(1);
         int DD = peps_1(rc,c).shape(3) * peps_2(rc,c).shape(3);
         int DR = peps_1(rc,c).shape(4) * peps_2(rc,c).shape(4);

         Contract((T)1.0,peps_1(rc,c),shape(i,j,k,l,m),peps_2(rc,c),shape(n,o,k,p,q),(T)0.0,tmp,shape(i,n,j,o,l,p,m,q));

         (*this)[c] = tmp.reshape_clear(shape(DL,DU,DD,DR));

      }

   }
   else{//vertical!

      this->resize(Ly);

      enum {i,j,k,l,m,n,o,p,q};

      TArray<T,8> tmp;

      for(int r = 0;r < Ly;++r){

         int DL = peps_1(r,rc).shape(3) * peps_2(r,rc).shape(3);
         int DU = peps_1(r,rc).shape(4) * peps_2(r,rc).shape(4);
         int DD = peps_1(r,rc).shape(0) * peps_2(r,rc).shape(0);
         int DR = peps_1(r,rc).shape(1) * peps_2(r,rc).shape(1);

         Contract((T)1.0,peps_1(r,rc),shape(i,j,k,l,m),peps_2(r,rc),shape(n,o,k,p,q),(T)0.0,tmp,shape(l,p,m,q,i,n,j,o));

         (*this)[r] = tmp.reshape_clear(shape(DL,DU,DD,DR));

      }

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

template MPO<double>::MPO(char,int,const PEPS<double> &,const PEPS<double> &);
template MPO< complex<double> >::MPO(char,int,const PEPS< complex<double> > &,const PEPS< complex<double> > &);

template MPO<double>::MPO(const MPO<double> &);
template MPO< complex<double> >::MPO(const MPO< complex<double> > &);

template MPO<double>::~MPO();
template MPO< complex<double> >::~MPO();

template int MPO<double>::gD() const;
template int MPO< complex<double> >::gD() const;
