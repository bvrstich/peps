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

template<typename T>
Lattice PEPS<T>::lat;

template<typename T>
Random PEPS<T>::RN;

/**
 * construct constructs a standard PEPS object, note: be sure to initialize the Lattice object before calling the constructor
 * @param D_in cutoff virtual dimension
 */
template<typename T>
PEPS<T>::PEPS(int D_in) : vector< TArray<T,5> >(lat.gLx() * lat.gLy()) {

   D = D_in;

   int Lx = lat.gLx();
   int Ly = lat.gLy();
   int d = lat.gd();

   //corners first

   //r == 0 : c == 0
   (*this)[ lat.grc2i(0,0) ].resize(1,D,d,1,D);

   //r == 0 : c == L - 1
   (*this)[ lat.grc2i(0,Lx-1) ].resize(D,D,d,1,1);

   //r == L - 1 : c == 0
   (*this)[ lat.grc2i(Ly-1,0) ].resize(1,1,d,D,D);

   //r == L - 1 : c == L - 1
   (*this)[ lat.grc2i(Ly-1,Lx-1) ].resize(D,1,d,D,1);

   //sides:

   //r == 0
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ lat.grc2i(0,c) ].resize(D,D,d,1,D);

   //r == Ly - 1
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ lat.grc2i(Ly-1,c) ].resize(D,1,d,D,D);

   //c == 0
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ lat.grc2i(r,0) ].resize(1,D,d,D,D);

   //c == Lx - 1
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ lat.grc2i(r,Lx - 1) ].resize(D,D,d,D,1);

   //the rest is full
   for(int r = 1;r < Ly - 1;++r)
      for(int c = 1;c < Lx - 1;++c)
         (*this)[ lat.grc2i(r,c) ].resize(D,D,d,D,D);

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c)
         for(typename TArray<T,5>::iterator it = (*this)[ lat.grc2i(r,c) ].begin();it != (*this)[lat.grc2i(r,c) ].end();++it)
            *it = rgen();

}

/**
 * copy constructor
 */
template<typename T>
PEPS<T>::PEPS(const PEPS<T> &peps_copy) : vector< TArray<T,5> >(peps_copy) {

   D = peps_copy.gD();

}

/**
 * empty destructor
 */
template<typename T>
PEPS<T>::~PEPS(){ }

//!function which generates random complex numbers uniformly on a square of side 2
template<>
complex<double> PEPS< complex<double> >::rgen(){ 

   return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); 

}

//!function which generates uniform random numbers between [-1:1]
template<>
double PEPS<double>::rgen(){ 

   return 2.0*RN() - 1.0;

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
const TArray<T,5> &PEPS<T>::operator()(int r,int c) const {

   return (*this)[lat.grc2i(r,c)];

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
TArray<T,5> &PEPS<T>::operator()(int r,int c) {

   return (*this)[lat.grc2i(r,c)];

}

/**
 * @return the cutoff virutal dimension
 */
template<typename T>
int PEPS<T>::gD() const {

   return D;

}

/**
 * @param peps_i peps to take the overlap with
 * @param D_aux auxiliary dimension of the contraction (determines the accuracy of the contraction)
 * @return the inner product of two PEPS <psi1|psi2> 
 */
template<typename T>
T PEPS<T>::dot(const PEPS<T> &peps_i,int D_aux) const {

   //start from bottom
   MPS<T> mps_b('b',*this,peps_i);

   //make it left canonicalized
   mps_b.canonicalize(Left,false);

   for(int i = 1;i < PEPS<T>::lat.gLy()/2;++i){

      //i'th row as MPO
      MPO<T> mpo(i,*this,peps_i);

      //apply to form MPS with bond dimension D^4
      mps_b.gemv('L',mpo);

      //reduce the dimensions of the edge states using thin svd
      mps_b.cut_edges();

      MPS<T> mps_c(mps_b.size());

      //compress using svd only
      mps_c.compress(D_aux,mps_b,5);

      mps_b = std::move(mps_c);

   }

   //then from top 
   MPS<T> mps_t('t',*this,peps_i);

   for(int i = PEPS<T>::lat.gLy() - 2;i >= PEPS<T>::lat.gLy()/2;--i){

      //i'th row as MPO
      MPO<T> mpo(i,*this,peps_i);

      //apply to form MPS with bond dimension D^4
      mps_t.gemv('U',mpo);

      //reduce the dimensions of the edge states using thin svd
      mps_t.cut_edges();

      MPS<T> mps_c(mps_t.size());

      //compress using svd only
      mps_c.compress(D_aux,mps_t,5);

      mps_t = std::move(mps_c);

   }

   return mps_b.dot(mps_t);

}

//forward declarations for types to be used!
template PEPS<double>::PEPS(int);
template PEPS< complex<double> >::PEPS(int);

template PEPS<double>::PEPS(const PEPS<double> &);
template PEPS< complex<double> >::PEPS(const PEPS< complex<double> > &);

template PEPS<double>::~PEPS();
template PEPS< complex<double> >::~PEPS();

template TArray<double,5> &PEPS<double>::operator()(int r,int c);
template TArray<complex<double>,5> &PEPS< complex<double> >::operator()(int r,int c);

template const TArray<double,5> &PEPS<double>::operator()(int r,int c) const;
template const TArray<complex<double>,5> &PEPS< complex<double> >::operator()(int r,int c) const;

template double PEPS<double>::dot(const PEPS<double> &peps_i,int D_aux) const;

template int PEPS<double>::gD() const;
template int PEPS< complex<double> >::gD() const;
