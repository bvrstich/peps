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
 * construct constructs DLPEPS object by contracting the physical dimension of two PEPS,
 */
template<typename T>
DLPEPS<T>::DLPEPS(const PEPS<T> &peps_1,const PEPS<T> &peps_2) : vector< TArray<T,4> >(PEPS<T>::lat.gLx() * PEPS<T>::lat.gLy()) {

   D = peps_1.gD() * peps_2.gD();

   int Lx = PEPS<T>::lat.gLx();
   int Ly = PEPS<T>::lat.gLy();

   //corners first

   //r == 0 : c == 0
   (*this)[ PEPS<T>::lat.grc2i(0,0) ].resize(1,D,1,D);

   //r == 0 : c == L - 1
   (*this)[ PEPS<T>::lat.grc2i(0,Lx-1) ].resize(D,D,1,1);

   //r == L - 1 : c == 0
   (*this)[ PEPS<T>::lat.grc2i(Ly-1,0) ].resize(1,1,D,D);

   //r == L - 1 : c == L - 1
   (*this)[ PEPS<T>::lat.grc2i(Ly-1,Lx-1) ].resize(D,1,D,1);

   //sides:

   //r == 0
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ PEPS<T>::lat.grc2i(0,c) ].resize(D,D,1,D);

   //r == Ly - 1
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ PEPS<T>::lat.grc2i(Ly-1,c) ].resize(D,1,D,D);

   //c == 0
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ PEPS<T>::lat.grc2i(r,0) ].resize(1,D,D,D);

   //c == Lx - 1
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ PEPS<T>::lat.grc2i(r,Lx - 1) ].resize(D,D,D,1);

   //the rest is full
   for(int r = 1;r < Ly - 1;++r)
      for(int c = 1;c < Lx - 1;++c)
         (*this)[ PEPS<T>::lat.grc2i(r,c) ].resize(D,D,D,D);

   //now fill it!
   TArray<T,8> tmp;

   enum {i,j,k,l,m,n,o,p,q};

   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         tmp.clear();

         Contract((T)1.0,peps_1(r,c),shape(i,j,k,l,m),peps_2(r,c),shape(n,o,k,p,q),(T)0.0,tmp,shape(i,n,j,o,l,p,m,q));

         CopyR(tmp,(*this)(r,c));

      }

}

/**
 * copy constructor
 */
template<typename T>
DLPEPS<T>::DLPEPS(const DLPEPS<T> &peps_copy) : vector< TArray<T,4> >(peps_copy) { }

/**
 * empty destructor
 */
template<typename T>
DLPEPS<T>::~DLPEPS(){ }

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
const TArray<T,4> &DLPEPS<T>::operator()(int r,int c) const {

   return (*this)[PEPS<T>::lat.grc2i(r,c)];

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
TArray<T,4> &DLPEPS<T>::operator()(int r,int c) {

   return (*this)[PEPS<T>::lat.grc2i(r,c)];

}

//forward declarations for types to be used!
template DLPEPS<double>::DLPEPS(const PEPS<double> &,const PEPS<double> &);
template DLPEPS< complex<double> >::DLPEPS(const PEPS< complex<double> > &,const PEPS< complex<double> > &);

template DLPEPS<double>::DLPEPS(const DLPEPS<double> &);
template DLPEPS< complex<double> >::DLPEPS(const DLPEPS< complex<double> > &);

template DLPEPS<double>::~DLPEPS();
template DLPEPS< complex<double> >::~DLPEPS();

template TArray<double,4> &DLPEPS<double>::operator()(int r,int c);
template TArray<complex<double>,4> &DLPEPS< complex<double> >::operator()(int r,int c);

template const TArray<double,4> &DLPEPS<double>::operator()(int r,int c) const;
template const TArray<complex<double>,4> &DLPEPS< complex<double> >::operator()(int r,int c) const;
