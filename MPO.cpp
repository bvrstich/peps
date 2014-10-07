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
 * empty constructor
 */
template<typename T>
MPO<T>::MPO() : vector< TArray<T,4> >() { }

/** 
 * constructor: just sets the length of the vector, nothing is allocates or initialized
 * @param L_in length of the chain
 */
template<typename T>
MPO<T>::MPO(int L_in) : vector< TArray<T,4> >(L_in) { }

/** 
 * standard constructor:
 * @param L_in length of the chain
 * @param d_phys_in physical dimension
 * @param D_in virtual max bond dimension
 * allocates the tensors with correct dimensions
 */
template<typename T>
MPO<T>::MPO(int L_in,int d_phys_in,int D_in) : vector< TArray<T,4> >(L_in) {

   D = D_in;
   d_phys = d_phys_in;

   vector<int> vdim(L_in + 1);

   vdim[0] = 1;

   for(int i = 1;i < L_in;++i){

      int tmp = vdim[i - 1] * d_phys * d_phys;

      if(tmp < D)
         vdim[i] = tmp;
      else 
         vdim[i] = D;

   }

   vdim[L_in] = 1;

   for(int i = L_in - 1;i > 0;--i){

      int tmp = vdim[i + 1] * d_phys * d_phys;

      if(tmp < vdim[i])
         vdim[i] = tmp;

   }

   for(int i = 0;i < this->size();++i)
      (*this)[i].resize(vdim[i],d_phys,d_phys,vdim[i+1]);

}

/**
 * copy constructor
 */
template<typename T>
MPO<T>::MPO(const MPO<T> &mpo_copy) : vector< TArray<T,4> >(mpo_copy) {

   D = mpo_copy.gD();
   d_phys = mpo_copy.gd_phys();

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

/**
 * @return physical dimension of the MPO
 */
template<typename T>
int MPO<T>::gd_phys() const {

   return d_phys;

}

/**
 * @param bra the bra of the inner product
 * @return the inner product of two MPS's, with *this being the ket
 */
template<typename T>
T MPO<T>::dot(const MPO<T> &bra) const {
/*
   TArray<T,2> E;

   Contract((T)1.0,bra[0],shape(0,1),(*this)[0],shape(0,1),(T)0.0,E);

   TArray<T,3> I;

   for(int i = 1;i < this->size();++i){

      I.clear();

      Contract((T)1.0,bra[i],shape(0),E,shape(0),(T)0.0,I);

      E.clear();

      Contract((T)1.0,I,shape(2,0),(*this)[i],shape(0,1),(T)0.0,E);

   }

   return E(0,0);
*/
   return (T) 0.0;
}

/**
 * Fill the MPO by contracting a the physical dimensions of a peps object
 * @param option 'b'ottom, 't'op, 'l'eft or 'r'ight
 * @param peps input PEPS<double> object
 */
template<>
void MPO<double>::fill(const char option,const PEPS<double> &peps){

   if(option == 'b'){

      enum {i,j,k,l,m,n,o,p,q};

      DArray<8> tmp;

      //share the pointer
      for(int col = 0;col < Lx;col++){

         tmp.share_mem( (*this)[col] );

         Contract(1.0,peps(0,col),shape(i,j,k,l,m),peps(0,col),shape(n,o,k,p,q),0.0,tmp,shape(i,n,j,o,l,p,m,q));

      }

   }
   else if(option == 't'){


   }
   else if(option == 'l'){


   }
   else{

   }

}

template MPO<double>::MPO();
template MPO< complex<double> >::MPO();

template MPO<double>::MPO(int);
template MPO< complex<double> >::MPO(int);

template MPO<double>::MPO(int,int,int);
template MPO< complex<double> >::MPO(int,int,int);

template MPO<double>::MPO(const MPO<double> &);
template MPO< complex<double> >::MPO(const MPO< complex<double> > &);

template MPO<double>::~MPO();
template MPO< complex<double> >::~MPO();

template int MPO<double>::gD() const;
template int MPO< complex<double> >::gD() const;

template double MPO<double>::dot(const MPO<double> &bra) const;
template  complex<double>  MPO< complex<double> >::dot(const MPO< complex<double> > &bra) const;
