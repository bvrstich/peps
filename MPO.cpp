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
 * construct constructs a standard MPO object, by creating a double layer peps object
 * @param option 'V' == vertical stripe, 'H' == horizontal stripe
 * @param rc row or column index
 */
template<typename T>
MPO<T>::MPO(char option,int rc,const PEPS<T> &peps_1,const PEPS<T> &peps_2) : vector< TArray<T,4> >() {

   D = peps_1.gD() * peps_2.gD();

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   if(option == 'H'){//horizontal

      this->resize(Lx);

      enum {i,j,k,l,m,n,o,p,q};

      TArray<T,8> tmp;

      //c == 0
      Contract((T)1.0,peps_1(rc,0),shape(i,j,k,l,m),peps_2(rc,0),shape(n,o,k,p,q),(T)0.0,tmp,shape(i,n,j,o,l,p,m,q));

      (*this)[0] = tmp.reshape_clear(shape(1,D,D,D));

      //c == 1 -> L - 2
      for(int c = 1;c < Lx - 1;++c){

         Contract((T)1.0,peps_1(rc,c),shape(i,j,k,l,m),peps_2(rc,c),shape(n,o,k,p,q),(T)0.0,tmp,shape(i,n,j,o,l,p,m,q));

         (*this)[c] = tmp.reshape_clear(shape(D,D,D,D));

      }

      Contract((T)1.0,peps_1(rc,Lx-1),shape(i,j,k,l,m),peps_2(rc,Lx-1),shape(n,o,k,p,q),(T)0.0,tmp,shape(i,n,j,o,l,p,m,q));

      (*this)[Lx-1] = tmp.reshape_clear(shape(D,D,D,1));

   }
   else{//vertical!

      this->resize(Ly);

      enum {i,j,k,l,m,n,o,p,q};

      TArray<T,8> tmp;

      //r == 0
      Contract((T)1.0,peps_1(0,rc),shape(i,j,k,l,m),peps_2(0,rc),shape(n,o,k,p,q),(T)0.0,tmp,shape(l,p,i,n,m,q,j,o));

      (*this)[0] = tmp.reshape_clear(shape(1,D,D,D));

      //r == 1 -> L - 2
      for(int r = 1;r < Ly - 1;++r){

         Contract((T)1.0,peps_1(r,rc),shape(i,j,k,l,m),peps_2(r,rc),shape(n,o,k,p,q),(T)0.0,tmp,shape(l,p,i,n,m,q,j,o));

         (*this)[r] = tmp.reshape_clear(shape(D,D,D,D));

      }

      Contract((T)1.0,peps_1(Ly-1,rc),shape(i,j,k,l,m),peps_2(Ly-1,rc),shape(n,o,k,p,q),(T)0.0,tmp,shape(l,p,i,n,m,q,j,o));

      (*this)[Ly-1] = tmp.reshape_clear(shape(D,D,D,1));

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
